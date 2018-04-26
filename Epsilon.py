
# coding: utf-8

# In[1]:

import numpy as nmp
from numpy import pi as PI
import math
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import re
import time
import sys

help_msg = "The script takes a PDB ID, a corresponding PQR file and a cube-format epsilon map from Delphi to write out the average epsilon as a function of distance normalized by the Radius of gyration as calculated from the input structure (PDB/PQR) file. However, that will only happen if do_avgEps_per_normDist ==  TRUE. If it is false, then the avergae eps around a user-specied 3D point with a radius of 1.0 A (default) will be returned."

syntax   = "Arguments: \n\t<PDB_ID>  \n\t<PQR file> \n\t<CUBE format EPS file> \n\t< do_avgEps_per_normDist? TRUE | FALSE > \n\t< 3D coord of a point around which avergae epsilon is sought with white space separation (DEFAULt: 0 0 0) > \n\t< radius of the shell (Ang) around that point (DEFAULT = 1.0) > \n\t<LOG FILE NAME> "

if len(sys.argv) != 10:
    print("number of Arguments = ",len(sys.argv))
    print("\n\n",help_msg)
    print(syntax,"\n\n")
    sys.exit()

#output tags
infotag = "INFO>"
epstag  = "AVGEPS>"
timetag = "TIME>"
warntag = "WARNING>"

#inputs
pid = sys.argv[1]
fpqr = sys.argv[2]
feps = sys.argv[3]
bAvgEps_nDist = sys.argv[4].capitalize().strip()
sb_center = sys.argv[5] + " " + sys.argv[6] + " " + sys.argv[7]
dR = float(sys.argv[8])
outname = sys.argv[9]

pid = pid.upper()

bAvgEps_nDist =  True if bAvgEps_nDist == "True" else  False

# In[2]:

################################ HOUSE operations ########################################
# cosines, phis = nmp.meshgrid(cosines, phis)
v=nmp.linspace(0,2*nmp.pi,12)
u=nmp.linspace(-1,1, 12)
u,v=nmp.meshgrid(u,v)

u = u.flatten()
v = v.flatten()

points2D = nmp.vstack([u,v]).T
triangulation = Delaunay(points2D)
ntri = triangulation.simplices.shape[0]
print("{:<8s} No. of triangles = {}".format(infotag, ntri))


# In[3]:

mass_table = {"C" : 12.01, 
             "H"  : 1.008,
             "CL" : 35.45,
             "NA" : 22.99,
             "S" : 32.06,
             "N" : 14.01,
             "O" : 16.0}


# In[4]:

def get_Rgyr(fpqr):
    #assuming max number of heavy atoms is 5000 for input structure
    maxAtoms = 5000
    coords = nmp.zeros( (maxAtoms,4)) # x,y,z, mass (look up table)
    numHatoms = 0
    
    with open(fpqr,'r') as f_in:
        for line in f_in:
            line = line.strip()
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atomname = line[12:16]
                atomtype = re.sub(r'^\d',"",atomname)[0:1]
                
                if not atomname.startswith("H"):
                    coords[numHatoms, :] = float(line[30:38]), float(line[38:46]), float(line[46:54]),mass_table[atomtype]
                    numHatoms += 1
    
    if numHatoms < maxAtoms:
        print("{:<8s} Will delete extra memory reserved for the coordiantes".format(warntag))
        print("{:<8s} Found {}  in the structure".format(infotag, numHatoms))
        coords = nmp.resize( coords, (numHatoms, 4) )
        
    #measure center
    center = [ nmp.mean(coords[:,0]), nmp.mean(coords[:,1]), nmp.mean(coords[:,2])]
    mass_tot = nmp.sum(coords[:,3])
    print("{:<8s} Center is : {:.2f} {:.2f} {:.2f}".format(infortag, *center))
    print("{:<8s} Max coords are: {:.2f} {:.2f} {:.2f}".format(infotag, nmp.max(coords[:,0]), nmp.max(coords[:,1]), nmp.max(coords[:,2])))
    
    # redefine coords and normalize mass
    for i in range(0,3):
        coords[:,i] -= center[i]
    
    coords[:,3] /= mass_tot   
    
    #calculate Rgyr
    R_sqrd = nmp.zeros( (len(coords[:,3]), 2) ) # only stores R2 and m
    R_sqrd[:,0], R_sqrd[:,1] = nmp.power(coords[:,0],2) + nmp.power(coords[:,1],2) + nmp.power(coords[:,2],2), coords[:,3]
    nmp.resize(coords, (0,0) )
    
    Rg = nmp.sum(R_sqrd[:,1] * R_sqrd[:,0])
    Rg = nmp.sqrt(Rg)  
    print("{:<8s} Radius of gyration = {:.3f}".format(infotag, Rg))
    
    maxR = nmp.sqrt(nmp.max(R_sqrd[:,0]))
    print("{:<8s} Max R = {:.3f}".format(infotag, maxR))
    return Rg, maxR

################################ HOUSE operations ########################################



# In[5]:

# feps = "/common/compbio/ARGO/Gaussian_multiDielectric/allPDB_noW_emSD5000_PQR/2NLS_noW_emSD5000.pqr_epsMap.cube"
# fpqr = "/common/compbio/ARGO/Gaussian_multiDielectric/allPDB_noW_emSD5000_PQR/2NLS_noW_emSD5000.pqr"


# In[6]:
if bAvgEps_nDist:
    Rgyr, Rmax = get_Rgyr(fpqr)


# In[7]:

start1 = time.time()

RBohr = 0.529177
NGrids = 0
dGrids = 0.0
origin, mid = [], []
triangles = []
lines_read = 0
points_read = 0
ix, iy, iz = 0, 0, 0

with open(feps) as epscube:
    for line in epscube:
        lines_read += 1
        if lines_read == 3:
            origin = [ float(oc)*RBohr for oc in line.split()[1:4] ]
            
        elif lines_read == 4:
            NGrids, dGrids = int(line.split()[0]), float(line.split()[1])*RBohr
            print("{:<8s} {} grid points in each direction".format( infotag, NGrids))
            epsmap = nmp.zeros( (NGrids*NGrids*NGrids,3) )   # eps,r2,dTringle_idx
            mid = [ oc + ((NGrids - 1)*dGrids/2) for oc in origin]
            maxCorner = [oc + (NGrids-1)*dGrids for oc in origin]

            
        elif lines_read > 7:
            for eps in line.split():
                ez = origin[2] + (dGrids * iz)
                ey = origin[1] + (dGrids * iy)
                ex = origin[0] + (dGrids * ix)
                
                r2 = (ex - mid[0])**2 + (ey-mid[1])**2 + (ez-mid[2])**2
                if r2 != 0 : 
                    cth = (ez-mid[2])/math.sqrt(r2)
                    phi = math.atan2((ey-mid[1]),(ex-mid[0]))
                    if ey < mid[1]:
                        phi += (2 * PI)
                    
                    dTriangle_idx = int(triangulation.find_simplex( [cth, phi] ))
                    epsmap[points_read, : ] = [float(eps), r2, dTriangle_idx ]
                    
                    
                points_read += 1
                # if points_read % 50000 == 0: print("{} points read. DT = {}".format(points_read,dTriangle_idx))
                if points_read % 50000 == 0: print("{:<8s} {:>10d} of {} points read".format(infotag, points_read,NGrids*NGrids*NGrids))
                iz += 1
                if iz == (NGrids):
                    iy += 1
                    iz = 0
                if iy == (NGrids):
                    ix += 1
                    iy = 0
                    
# print(ix, iy, iz)
end1 = time.time()

print("{:<8s} Takes {:.3f} mins to read the eps cube file.".format(timetag, (end1-start1)/60))
                


# In[8]:

if not bAvgEps_nDist:
        
        fout = open(outname,"w")

        print("{:<8s} Will only calculate the avg eps in the vicinity of {}".format(infotag, sb_center))
        fout.write("{:<8s} Will only calculate the avg eps in the vicinity of {}\n".format(infotag, sb_center))
        
        inx, iny, inz = [float(cc) for cc in sb_center.split()]
        
        inRange = True
        if inx >= maxCorner[0] or iny >= maxCorner[1] or inz >= maxCorner[2] or inx <= origin[0] or iny <= origin[1] or inz <= origin[2]:
            print("{:<8s} Bad Input value for the Coordinates. Point outside of cube range.".format(warntag))
            inRange = False

        if inRange:
            inr = math.sqrt((inx - mid[0])**2 + (iny - mid[1])**2 + (inz - mid[2])**2)

            incth = (inz-mid[2])/inr
            inphi = math.atan2(iny - mid[1],inx - mid[0])
            if iny < mid[1]: inphi += 2*PI

            target_id = int(triangulation.find_simplex( [incth, inphi] ))
            print("{:<8s} Target matched in simplex number {}".format(infotag,target_id))
            fout.write("{:<8s} Target matched in simplex number {}\n".format(infotag,target_id))

            r = []  #distance I guess
            e = []  #epsilon

            # er2d stands for "eps, r2, dtriangle_idx"
            for row in filter(lambda er2d: er2d[2] == target_id and math.sqrt(er2d[1]) <= (inr + dR) and math.sqrt(er2d[1]) >= (inr - dR), epsmap):
                r.append(math.sqrt(row[1]))
                e.append(row[0])

            print("{:<8s} Minimium and maximum epsilons are {} and {}".format(epstag, min(e),max(e)))
            fout.write("{:<8s} Minimium and maximum epsilons are {} and {}\n".format(epstag, min(e),max(e)))
            print("{:<8s} Avergae epsilon in the vicinity: {:.3f}\n".format(epstag, nmp.mean(e)))
            fout.write("{:<8s} Avergae epsilon in the vicinity: {:.3f}\n".format(epstag, nmp.mean(e)))

        fout.close()

# In[10]:

if bAvgEps_nDist:
        print("{:<8s} Will calculate the avg eps from r = 1 to r = {:.3f}".format(infotag, math.floor(Rmax)))
        start2 = time.time()
        avg_eps = []
        norm_dist = []
        for inr in range(1,math.floor(Rmax + 1)):
            sum_eps = 0
            num_eps = 0
            dr = 1.0
            num_scans = 0
            for incth in nmp.linspace(-1,1, 12):
                for inphi in nmp.linspace(0,2*nmp.pi,12):
                    target_id = int(triangulation.find_simplex( [incth, inphi] ))
                    # print("Target matched in simplex number {}".format(target_id))
                    
                    # er2d stands for "eps, r2, dtriangle_idx"
                    for row in filter(lambda er2d: er2d[2] == target_id and math.sqrt(er2d[1]) <= (inr + dr) and math.sqrt(er2d[1]) >= (inr - dr), epsmap):
                        sum_eps += row[0]
                        num_eps += 1
                    
                    num_scans += 1
                    if num_scans % 24 == 0:
                        print("{:<8s} {} combinations scanned".format(infotag, num_scans))
                    # print("Minimium and maximum epsilons are {} and {}".format(min(e),max(e)))
                    
            avg_eps.append(sum_eps/num_eps)
            norm_dist.append(inr/Rgyr)
                    
        end2 = time.time()

        print("{:<8s} Takes {} mins to calculate avg_eps per distance.".format(timetag, (end2-start2)/60))


        #print output
        with open(outname, 'w') as fout:
            fout.write("#NORMDIST\tAVGEPS\n")
            for data in zip(norm_dist, avg_eps):
                fout.write("{:.1f}\t{:.3f}\n".format(*data))



        # In[ ]:

        # plt.plot(norm_dist,avg_eps,'o')
        # plt.xlabel('Norm Distance')
        # plt.ylabel('Epsilon')
        # plt.show()


# In[ ]:



