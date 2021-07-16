import sys

BX = [sys.argv[1], sys.argv[2]]

bxList = []
with open(BX[0], 'r') as f1:
	for l1 in f1:
		l = l1.split()
		bxList.append([l[0], int(l[1])])
f1.close()

idx = 0
Cnt = len(bxList)
ReadPair = 0
with open(BX[1], 'r') as f:
	for line in f:
		bxInfo = line.split()		
		while idx < Cnt and bxList[idx][0] < bxInfo[0]:
			idx += 1
		if idx >= Cnt:
			break
		if bxList[idx][0] == bxInfo[0]:
			ReadPair += bxList[idx][1] * int(bxInfo[1])
			idx += 1


print BX[0].split("/")[-2], ReadPair
	
