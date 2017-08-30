import os

os.chdir("/Users/arghyachakravorty/Box Sync/Gaussian_multiDielectric/allPDB_MD23_20ns_EDR")

run = 2
with open("../common/pdb_id.list","r") as fPDB:
    pids = [pid.strip() for pid in fPDB]

print("A total of {} ids read".format(len(pids)))

# pid = "1aho"
# 

for pid in pids:
    pid = pid.upper()
    run_num = str(run)
    
    xvg = pid + "_md"+ run_num + "_energy.xvg"
    
    out = xvg.replace(".xvg","_formatted.dat")
    fout = open(out, 'w')
    
    with open(xvg, 'r') as fxvg:
        for line in fxvg:
            if not (line.startswith("#") or line.startswith("@")):
                time, pe = line.split()
                fout.write("{}\t{}\t{}\n".format(time, run, pe))
                
    fout.close()
    print("Formatting XVG for {} is done".format(pid))

            
            
                