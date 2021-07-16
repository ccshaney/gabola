import sys
strand = ["+", "-"]

def countM(cigar, MappedRange, strand, node_len):
	cnt = ''
	M_length = 0
	node_pos = 0
        tempRange = [MappedRange[0], MappedRange[1]]
	for i in range(len(cigar)):
		if cigar[i].isalpha():
			if cigar[i]=='M':
				if strand == 0:
                                    tempRange[0] = min(MappedRange[0], node_pos + 1, tempRange[0])
                                    tempRange[1] = max(MappedRange[1], node_pos + int(cnt), tempRange[1])
				else:
                                    tempRange[0] = min(MappedRange[0], node_len - (node_pos + int(cnt)) + 1, tempRange[0])
                                    tempRange[1] = max(MappedRange[1], node_len - node_pos, tempRange[1])

				M_length += int(cnt)
				node_pos += int(cnt)
			elif cigar[i] == 'D':
				M_length += int(cnt)
			elif cigar[i] in ['S', 'H', 'I']:
				node_pos += int(cnt)
			cnt = ''
		else:
			cnt = cnt + cigar[i]
	return [M_length, tempRange]


samname = [sys.argv[1], sys.argv[3]]
part = [sys.argv[2], sys.argv[4]]
ContigMappedLen = int(sys.argv[5])
ContigCIGAR = float(sys.argv[6])
status = {}
EndContig = []
for z in range(2):
	Seen = []
	with open(samname[z], "r") as file:
		for line in file:
			collect = 0
			x = line.split()
			if "@SQ" == x[0]:
				scafName = x[1][3:]
				scafLen = int(x[2][3:])
				continue
			if "@PG" == x[0] or x[2] == "*": continue

			if x[0] not in Seen:
				Seen.append(x[0])
			else:
				continue

			contigLen = int(x[0].split("size")[1])
			#contigCov = float(x[0].split("_")[5])
			if contigLen < 1000 : continue
			s = (int(x[1])&16)/16

			match = [0] * 2
			matchRange = [[float("inf"), float("-inf")], [float("inf"), float("-inf")]]
			temp = countM(x[5], matchRange[s], s, contigLen)
			match[s] = temp[0]
			matchRange[s] = temp[1]
			end = int(x[3])+ match[s] - 1

			SA = ['', '']

			if len(x) < 16 or x[15][:5] != 'SA:Z:':
				if match[s] > 500:
					C = float(match[s])/(end-int(x[3])+1)
					if part[z] == "H" and int(x[3]) < 20000 :
						if match[s] >= ContigMappedLen and C >= ContigCIGAR/100:
							co = [x[0], strand[s], x[3], end, SA[s], int(x[3]) - 1, scafName + ' H', matchRange[s]]
							collect = 1
					elif part[z] == "T" and abs(end - scafLen) < 20000:
							if match[s] >= ContigMappedLen and C >= ContigCIGAR/100:
								co = [x[0], strand[s], x[3], end, SA[s], abs(end - scafLen), scafName + ' T', matchRange[s]]
								collect = 1
					if collect:
						if not x[0] in status:
							status[x[0]] = [0, 0]
						status[x[0]][z] = 1
						EndContig.append(co)

				continue

			contigRange = [[float("inf"), float("-inf")], [float("inf"), float("-inf")]]
			contigRange[s][0] = int(x[3])
			contigRange[s][1] = contigRange[s][0] + match[s]

			aln = x[15][5:-1].split(";")
			for i, info in enumerate(aln):
				read = info.split(",")
				aln_s = strand.index(read[2])
				temp = countM(read[3], matchRange[aln_s], aln_s, contigLen)
				M = temp[0]
				if M < 300: continue
				start = int(read[1])
				end = int(read[1]) + M - 1

				if start < contigRange[aln_s][0]:
					if abs(end - contigRange[aln_s][0]) < contigLen:
						contigRange[aln_s][0] = start
                                                matchRange[aln_s] = temp[1]
                                                match[aln_s] += M
					elif M > (contigRange[aln_s][1] - contigRange[aln_s][0]):
						contigRange[aln_s][0] = start
						contigRange[aln_s][1] = end
						SA[aln_s] = info+","+str(end)+";"
                                                temp = countM(read[3], [float("inf"), float("-inf")], aln_s, contigLen)
                                                matchRange[aln_s] = temp[1]
                                                match[aln_s] = temp[0]
						continue

				if end > contigRange[aln_s][1]:
					if abs(start - contigRange[aln_s][1])  < contigLen:
						contigRange[aln_s][1] = end
                                                matchRange[aln_s] = temp[1]
                                                match[aln_s] += M
					elif M > (contigRange[aln_s][1] - contigRange[aln_s][0]):
						contigRange[aln_s][0] = start
						contigRange[aln_s][1] = end
						SA[aln_s] = info+","+str(end)+";"
                                                temp = countM(read[3], [float("inf"), float("-inf")], aln_s, contigLen)
                                                matchRange[aln_s] = temp[1]
                                                match[aln_s] = temp[0]
						continue

				if start >= contigRange[aln_s][0] and end <= contigRange[aln_s][1]:
					SA[aln_s] += info+","+str(end)+";"



			if match[0] > match[1]:

				if match[0] > 500:
					C = float(match[0])/(contigRange[0][1]-contigRange[0][0]+1)
					if part[z] == "H" and contigRange[0][0] < 20000:
						if contigRange[0][1]-contigRange[0][0] >= ContigMappedLen and C >= ContigCIGAR/100:
							co = [x[0], "+", contigRange[0][0], contigRange[0][1], SA[0], contigRange[0][0]-1, scafName + ' H', matchRange[0]]
							collect = 1
					elif part[z] == "T" and abs(contigRange[0][1] - scafLen) < 20000:
						if contigRange[0][1]-contigRange[0][0] >= ContigMappedLen and C >= ContigCIGAR/100:
							co = [x[0], "+", contigRange[0][0], contigRange[0][1], SA[0], abs(contigRange[0][1] - scafLen), scafName + ' T', matchRange[0]]
							collect = 1
			else:

				if match[1] > 500:
					C = float(match[1])/(contigRange[1][1]-contigRange[1][0]+1)
					if part[z] == "H" and contigRange[1][0] < 20000 :
						if contigRange[1][1]-contigRange[1][0] >= ContigMappedLen and C >= ContigCIGAR/100:
							co = [x[0], "-", contigRange[1][0], contigRange[1][1], SA[1], contigRange[1][0]-1, scafName + ' H', matchRange[1]]
							collect = 1
					elif part[z] == "T" and abs(contigRange[1][1] - scafLen) < 20000:
						if contigRange[1][1]-contigRange[1][0] >= ContigMappedLen and C >= ContigCIGAR/100:
							co = [x[0], "-", contigRange[1][0], contigRange[1][1], SA[1], abs(contigRange[1][1] - scafLen), scafName + ' T', matchRange[1]]
							collect = 1
			if collect:
				if not x[0]	in status:
					status[x[0]] = [0, 0]
				status[x[0]][z] = 1
				EndContig.append(co)


for contig, value in status.items():
	if value[0] and value[1]:
		print samname, "Cross!"
		for s in EndContig:
			if s[0] == contig:
				print s[6]
				print s[0], s[1], s[2], s[3], s[7], str(s[5])+"_"+s[6][-1]
													
		print "----------------"


