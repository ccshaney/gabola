import sys
import os
wd = sys.argv[1]
gapCnt = int(sys.argv[2])
Length = int(sys.argv[3])
Coverage = float(sys.argv[4])

def getContig(filename, idx):
	skip = 0
	with open(filename, "r") as f:
		for l in f:
			y = l.strip()
			if '>' == y[0]:
				name = y.split("_")
				if float(name[5]) >= Coverage and int(name[3]) >= Length:
					skip = 0
					print y+"_"+str(idx)
					continue
                                else:
					skip = 1

			if not skip:
				print y

for i in range(1, gapCnt+1):
	fname = wd+"/scaffolds_gapID"+str(i)+".fasta"
	if not os.path.exists(fname): continue
	getContig(fname, i)

