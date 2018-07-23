import os
paths = ["./TC_ILP/","./RC_ILP/"]
for path in paths:
	outfile = open(path + "sol", "w")
	for file in os.listdir(path):
		if (file[-3:] == ".lp"):
			outfile.write("read "+file+"\n")
			outfile.write("opt"+"\n")
			outfile.write("write "+file[:-3]+".sol"+"\n")
			outfile.write("\n")