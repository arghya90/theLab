import os, sys

help_m = "This script takes a XVG file as an input and outputs a txt file whise name should also be supplied as an input."
syntax = "python readXVG_data.py <XVG file> <OUTPUT file name with extension>"

if len(sys.argv) != 3:
	print(help_m)
	print(syntax)
	sys.exit()

#input arguments
xvgn = sys.argv[1]
outn = sys.argv[2]

fout = open(outn, 'w')
fout.write("# Converted from " + xvgn + "\n")
with open(xvgn, 'r') as fxvg:
    for line in fxvg:
        if not (line.startswith("#") or line.startswith("@")):
            element_array = line.split()
            for ele in element_array:
            	fout.write("{}\t".format(ele))
            fout.write("\n")
            
fout.close()
print("Successfully completed converting " + xvgn + " --> " + outn)

            
            
                
