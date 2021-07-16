import sys

wd = sys.argv[1]
Rp = sys.argv[2]
minBXListSize = int(sys.argv[3])
minFlnkingRp = int(sys.argv[4])
scaf = wd.replace('/', '').split('scaf')[-1]
gapPOSFILE = wd + '/gap_pos_scaf' + scaf + '.txt' 
MappedReadFILE = wd + '/scaffold' + scaf + '_rp' + Rp + 'BX_C70M60ReadPos.tsv'
outRangeFILE = wd + '/gapRange_BXs.txt'
#Reads' end position is sorted by decreasing order 

def CollectBX():
	if not gapInfo: 
		print "scaf"+scaf, "has no gap."
		return
        cnt = 0
        totalCnt = len(gapInfo)
        with open(MappedReadFILE, "r") as bfile:
                for line in bfile:
                        x = line.split()
                        barcode = x[0]
                        pos1 = int(x[1])
			pos2 = int(x[2])
                        if pos1 > gapInfo[cnt][3]:
                                continue
                        while cnt < totalCnt and pos2 < gapInfo[cnt][2]:
                                cnt += 1

                        if cnt >= totalCnt:
                                break

                        if pos1 <= gapInfo[cnt][3] and pos2 >= gapInfo[cnt][2]:
                                for next_cnt in range(cnt, totalCnt):
                                        if pos1 <= gapInfo[next_cnt][3] and pos2 >= gapInfo[next_cnt][2]:
                                                if len(gapInfo[next_cnt]) == 4:
                                                        gapInfo[next_cnt].append({barcode: 1})
                                                elif barcode in gapInfo[next_cnt][4]:
                                                        gapInfo[next_cnt][4][barcode] += 1
                                                else:
                                                        gapInfo[next_cnt][4][barcode] = 1
                                        elif pos1 > gapInfo[next_cnt][3] :
                                                break




def GapRange(gapStart, gapEnd):
        gapLength = gapEnd - gapStart + 1
        careRange = [5000, 10000, 15000, 20000]
        for r in careRange:
                if gapLength < r:
                        break

        return [1, gapEnd+r] if gapStart-r < 0 else [gapStart-r, gapEnd+r]


gapInfo = []
with open(gapPOSFILE, "r") as gfile:
        for line in gfile:
                g = line.split()
                gapStart = int(g[0])
                gapEnd = int(g[1])
                r = GapRange(gapStart, gapEnd)
                gapInfo.append([gapStart, gapEnd, r[0], r[1]])

#range sort by decreasing order
gapInfo = sorted(gapInfo, key=lambda x: (x[3], -x[2]), reverse=True)
rm = []
for i in range(len(gapInfo)-1):
	if i in rm: continue
	for j in range(i+1, len(gapInfo)):
		if j in rm: continue
		if gapInfo[i][2] <= gapInfo[j][2]:
			rm.append(j)
		elif gapInfo[i][2] >= gapInfo[j][3]:
			break
rm = sorted(rm, reverse=True)
for r in rm:
	del gapInfo[r]
CollectBX()

outFILE = open(outRangeFILE, "w")
idx = 0
for info in gapInfo:

        if len(info) > 4:
                bx_set = [bx for bx, value in info[4].items() if value >= 2 * minFlnkingRp]
		
		if len(bx_set) >= minBXListSize:
			idx += 1
			outFILE.write(str(idx)+" "+str(info[0])+" "+str(info[1])+" "+str(info[2])+" "+str(info[3])+" "+','.join(bx_set)+"\n")
        else:
                bx_set = []
                print "gap:", info[0], info[1], "has no barcode."

