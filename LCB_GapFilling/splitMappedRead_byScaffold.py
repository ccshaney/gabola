import sys
import os

RenameFile = sys.argv[1]
SamFile = sys.argv[2]
Num = sys.argv[3]
wd = os.path.dirname(RenameFile)
#RenameFile format
# >ref_scafname >scaffoldX|sizeXXX


def countCIGAR(cigar):
	cnt = ''
	mappedPos = 0
	for i in range(len(cigar)):
		if cigar[i].isalpha():
			if(cigar[i] == 'M'):
					mappedPos = mappedPos + int(cnt)
			if(cigar[i] == 'D'):
					mappedPos = mappedPos + int(cnt)
			cnt = ''
		else:
			cnt = cnt+cigar[i]

	return mappedPos

scafName = {}
with open(RenameFile, "r") as f:
    for line in f:
        l = line.split()
        scafName[l[0][1:]] = l[1].split("|")[0][1:]


preFile = ''
preContent = ''
with open(SamFile, "r") as sfile:
	for line in sfile:
                if line[0] == "@": continue

		readInfo = line.split()
                idx = int(scafName[readInfo[2]][8:])
                if idx > int(Num): continue

		C = countCIGAR(readInfo[5])
		endPos = int(readInfo[3]) + C - 1
                BX = readInfo[0].split("_")[1]
                outFile = wd + "/scaffold" + str(idx) + "_C70M60ReadPos.tsv"

                if outFile != preFile and preFile != '':
                    fd = open(preFile, "a")
		    fd.write(preContent)
                    fd.close()
                    preContent = ''
                
                preContent += BX + "\t" + readInfo[3] + "\t" + str(endPos) + "\n"
                preFile = outFile


if preContent != '':
    fd = open(preFile, "a")
    fd.write(preContent)
    fd.close()
