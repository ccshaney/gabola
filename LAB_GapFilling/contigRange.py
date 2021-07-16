import sys
fname=sys.argv[1]
#readName flag pos MAPQ CIGAR SAtag/none 
strand = ["+", "-"]

def countM(cigar):
	cnt = ''
	M_length = 0
	for i in range(len(cigar)):
		if cigar[i].isalpha():
				if(cigar[i]=='M' or cigar[i]=='D'):
						M_length = M_length + int(cnt)
				cnt = ''
		else:
				cnt = cnt+cigar[i]
	return M_length


with open(fname, "r") as file:
	for line in file:
		x = line.split()
		if "@SQ" in x[0] or "@PG" in x[0]: continue
		s = (int(x[1])&16)/16
		match = [0] * 2
		match[s] = countM(x[4])
		end = int(x[2])+ match[s] - 1

		SA = ['', '']
		SA[s] = ','.join(["scaffold", x[2], strand[s], x[4], x[3], "0", str(end)+";"])

		if len(x) < 6 or x[5][:5] != 'SA:Z:':
			#M = countM(x[4])
			if match[s] > 100:
				print x[0], strand[s], x[2], end, SA[s]
			continue

		contigRange = [[float("inf"), float("-inf")], [float("inf"), float("-inf")]]
		contigRange[s][0] = int(x[2])
		contigRange[s][1] = contigRange[s][0] + match[s]

		contigLen = int(x[0].split("_")[3])
		aln = x[5][5:-1].split(";")
		for i, info in enumerate(aln):
			read = info.split(",")
			aln_s = strand.index(read[2])
			M = countM(read[3])
			match[aln_s] += M
			start = int(read[1])
			end = int(read[1]) + M - 1

			if start < contigRange[aln_s][0]:
				if abs(end - contigRange[aln_s][0]) < contigLen:
					contigRange[aln_s][0] = start
				elif M > (contigRange[aln_s][1] - contigRange[aln_s][0]):
					contigRange[aln_s][0] = start
					contigRange[aln_s][1] = end
					SA[aln_s] = info+","+str(end)+";"
					continue

			if end > contigRange[aln_s][1]:
				if abs(start - contigRange[aln_s][1])  < contigLen:
					contigRange[aln_s][1] = end
				elif M > (contigRange[aln_s][1] - contigRange[aln_s][0]):
					contigRange[aln_s][0] = start
					contigRange[aln_s][1] = end
					SA[aln_s] = info+","+str(end)+";"
					continue

			if start >= contigRange[aln_s][0] and end <= contigRange[aln_s][1]:
				SA[aln_s] += info+","+str(end)+";"



		if match[0] > match[1]:
			if contigRange[0][1] - contigRange[0][0] > 100:
				print x[0], "+", contigRange[0][0], contigRange[0][1], SA[0]
		else:
			if contigRange[1][1] - contigRange[1][0] > 100:
				print x[0], "-", contigRange[1][0], contigRange[1][1], SA[1]



