from Bio.Seq import Seq
import sys
import os

#gap_coveredContig: gapStart gapEnd "compliete:" completelyCoveredContig "partialL:" LeftSideCoveredContig "partialR:" RightSideCoveredContig
gapFILE = sys.argv[1]

#contigRange: ContigName strand mappedStart mappedEnd SAtag
contigFILE = sys.argv[2]
wd = os.path.dirname(contigFILE)

#scaffoldXX.fa
scafFILE = sys.argv[3]

#local assembly fasta (scaffolds.fasta by SPAdes)
seqFILE = sys.argv[4]

#Threshold of mapped contigs
AlignedLen = int(sys.argv[5])
AlignedIdn = float(sys.argv[6])/100

#Record gap-filling info of each gap
out_fname = wd + '/gapStatus_Record.txt'

#Output fasta
filledScaf_fname = scafFILE.split(".fa")[0] + "_filled.fa"


filling = []
def getContigInfo(contigName):
	with open(contigFILE, "r") as f:
		for l in f:
			if contigName in l:
				return l.strip()

def countCigar(cigar):
	cnt = ''
	c_length = 0
	for i in range(len(cigar)):
		if cigar[i].isalpha():
				if(cigar[i] != 'D'):
						c_length = c_length + int(cnt)
				cnt = ''
		else:
				cnt += cigar[i]
	return c_length

def checkIdentity(cigar, threshold):
	cnt = ''
	matchCount = 0
	totalLength = 0
	for i in range(len(cigar)):
		if cigar[i].isalpha():
			if cigar[i] == 'M':
				matchCount += int(cnt)
			if cigar[i] != 'S' and cigar[i] != 'H' and cigar[i] != 'I':
				totalLength += int(cnt)
			cnt = ''
		else:
			cnt = cnt + cigar[i]

	if float(matchCount)/totalLength >= threshold:
		return 1
	else:
		return 0

def getContigSeq(contigName, strand):
	seq = ''
	with open(seqFILE) as fp:
		flag = 0
		for i, contig_line in enumerate(fp):
			if contigName in contig_line:
				flag = 1
				continue
			if flag:
				if ">" in contig_line:
					break
				seq += contig_line.strip()
	if seq == '':
		print "cant find " + contigName + " in fasta file"
		return seq

	if strand == "+":
		return seq
	else:
		return Seq(seq).reverse_complement()


pre_contig = ''
contig_seq = ''
with open(gapFILE, "r") as file:
	for line in file:
		status = "noContig"
		filled = 0
		gapInfo = line.split()
		gapStart = int(gapInfo[0])
		gapEnd = int(gapInfo[1])
                
                #Fully covered contig
		if gapInfo[3] != "{}":
			ccc = gapInfo[3][1:-1].split(";")[:-1]
                        contig = sorted(ccc, key=lambda x: (int(x.split("_")[3])), reverse=True)
			for contigName in contig:
				alnLeft = [0] * 7
				alnRight = [0] * 7
				inScaf = 0
				alignedScafPos1 = 0
				contigInfo = getContigInfo(contigName).split()
				alignment = contigInfo[4][:-1].split(";")
				strand = contigInfo[1]
				for x in alignment:
					aln = x.split(",")
					aln_start = int(aln[1])
					aln_end = int(aln[6])
					C = aln[3]
					scafCur = aln_start
					readCur = 0
                			#Gap locates in single alignment
					if gapStart >= aln_start and gapEnd <= aln_end: #inScaf
						inScaf = 1
						start = 0
						cnt = ''
						for i in range(len(C)):
							if C[i].isalpha():
								if C[i] != 'D': readCur += int(cnt)
								if C[i]=='M' or C[i]=='D':
									scafCur = scafCur + int(cnt)
									#get contig pos when gap starts
									if not start and scafCur - 1 >= gapStart:
										start = 1
										if C[i] == 'D': 
											if scafCur > gapEnd: #all gaps in deletion
												if (gapEnd - gapStart + 1) == int(cnt):
													filling.append([gapStart, gapStart, gapEnd, '', "Delete", 0, contigName])
													filled = 1
													print "Delete gaps:", str(gapStart), "Same length as deletion"
													alignedScafPos1 =  0
													break

												print "Gaps starting from", str(gapStart), "are inside deletion"
												status = "inDeletion"
												break

											alignedScafPos1 =  scafCur - int(cnt) 
											alignReadPos1 = readCur + 1
										else:
											alignedScafPos1 =  gapStart
											alignReadPos1 = readCur - (scafCur - 1 - gapStart)
									if scafCur > gapEnd:
										if C[i] == 'D':
												alignedScafPos2 = scafCur - 1
												alignReadPos2 = readCur
										else:
												alignedScafPos2 = gapEnd
												alignReadPos2 = readCur - (scafCur - 1 - gapEnd)
										
										filled = 1
										break

								cnt = ''
							else:
								cnt = cnt + C[i]


						break
                                        #Check if alignment is close enough to the gaps
					if abs(gapStart - aln_end) <= 50:
						if aln_end - aln_start > int(alnLeft[6]) - int(alnLeft[1]):
							alnLeft = aln
					elif abs(gapEnd - aln_start) <= 50:
						if aln_end - aln_start > int(alnRight[6]) - int(alnRight[1]):
							alnRight = aln

				if filled: break

                        	#Gap have contigs aligned on both ends
				if alnLeft[1] != 0 and alnRight[1] != 0:
					contigLen = int(contigName.split('_')[3])
					CLeft = countCigar(alnLeft[3].split("M")[-1])
					CRight = countCigar(alnRight[3].split("M")[0])
					if not checkIdentity(alnLeft[3], AlignedIdn) or not checkIdentity(alnRight[3], AlignedIdn):
						status = 'LowIdentity'
						continue
					
					if int(alnLeft[6]) + 1 <= gapStart:
						alignedScafPos1 = int(alnLeft[6]) + 1
						alignReadPos1 = contigLen - CLeft + 1
					else: #mismatch with N
						alignedScafPos1 = gapStart
						alignReadPos1 = contigLen - CLeft - (int(alnLeft[6]) - gapStart)

					if int(alnRight[1]) > gapEnd:
						alignedScafPos2 = int(alnRight[1]) - 1
						alignReadPos2 = CRight - 1
					else: #mismatch with N
						alignedScafPos2 = gapEnd
						alignReadPos2 = CRight + (gapEnd - int(alnRight[1])) + 1

					if alignReadPos2 - alignReadPos1 + 1 > 0:
						filled = 1
						break
					else:
						print "negative length! check out Gap position:", gapStart, "contig:", contigName
						status = "Conflict"
						continue
				elif not inScaf:
					print "Inconsistent contig:", contigName, "Gap position:", gapStart
					status = "Inconsistent"

			
			if alignedScafPos1 and filled:
				if pre_contig != contigName:
					contig_seq = getContigSeq(contigName, strand)
				pre_contig = contigName
				filling.append([gapStart, alignedScafPos1, alignedScafPos2, contig_seq[alignReadPos1-1: alignReadPos2], "Complete", gapEnd-gapStart+1, contigName])
				
				continue
		L = 0
		R = 0
                #Paritally covered contig on left side of gap
		if not filled and gapInfo[5] != "{}":
			ccc = gapInfo[5][1:-1].split(";")[:-1]
			contig = sorted(ccc, key=lambda x: (int(x.split("_")[3])), reverse=True)
			for contiggg in contig:
				contigInfo = getContigInfo(contiggg).split()
				contigLen = int(contiggg.split('_')[3])
				MappedLen = int(contigInfo[3]) - int(contigInfo[2]) + 1
				alignment = contigInfo[4][:-1].split(";")
				rightMost = 0
				tmep_max = float("-inf")
				for x in alignment:
					aln = x.split(",")
					if int(aln[6]) > rightMost:
						rightMost = int(aln[6])
						rightAln = aln
						
				C = countCigar(rightAln[3].split("M")[-1])

				#Nothing to fill
				L = 0 if C == 0 else 1

				#Contig at least 300 Match
				if MappedLen < AlignedLen:
					L = 0
					status = "PartialFailed_M"

			
				if int(rightAln[6]) < gapStart:
					alignReadPos1_L = contigLen - C + 1
					alignReadPos2_L = contigLen
					alignedScafPos1_L = int(rightAln[6]) + 1
					alignedScafPos2_L = int(rightAln[6]) + C
				else:
					alignReadPos1_L = contigLen - C - (int(rightAln[6]) - gapStart)
					alignReadPos2_L = contigLen
					alignedScafPos1_L = gapStart
					alignedScafPos2_L = int(rightAln[6]) + C

				#Alignment identity
				if L and not checkIdentity(rightAln[3], AlignedIdn):
					L = 0
					status = "PartialFailed_I"
					print "PartialL covered gap:", gapStart, "Contig Idnetity <", AlignedIdn, ":", contiggg

				#Can't reach gap region
				if alignedScafPos2_L < gapStart:
					L = 0

				if L:
					if pre_contig != contiggg:
						contig_seq = getContigSeq(contiggg, contigInfo[1])
						pre_contig = contiggg
					contig_seqL = contig_seq[alignReadPos1_L - 1:alignReadPos2_L]
					contigNameL = contiggg
					break
		#Paritally covered contig on right side of gap		
		if not filled and gapInfo[7] != "{}":
			ccc = gapInfo[7][1:-1].split(";")[:-1]
			contig = sorted(ccc, key=lambda x: (int(x.split("_")[3])), reverse=True)
			for contiggg in contig:
				contigInfo = getContigInfo(contiggg).split()
				contigLen = int(contiggg.split('_')[3])
				MappedLen = int(contigInfo[3]) - int(contigInfo[2]) + 1
				alignment = contigInfo[4][:-1].split(";")
				leftMost = int(contigInfo[3])
				tmep_min = float("inf")
				for x in alignment:
					aln = x.split(",")
					if int(aln[1]) < leftMost:
						leftMost = int(aln[1])
						leftAln = aln
						
				C = countCigar(leftAln[3].split("M")[0])
				R = 0 if C == 0 else 1

				#Contig with sufficient numbers of Match
				if MappedLen < AlignedLen:
					R = 0
					status = "PartialFailed_M"


				if int(leftAln[1]) > gapEnd:
					alignReadPos1_R = 1
					alignReadPos2_R = C
					alignedScafPos1_R = int(leftAln[1]) - C
					alignedScafPos2_R = int(leftAln[1]) - 1
				else:
					alignReadPos1_R = 1
					alignReadPos2_R = C + (gapEnd - int(leftAln[1])) + 1
					alignedScafPos1_R = int(leftAln[1]) - C
					alignedScafPos2_R = gapEnd
					
				#Alignment identity
				if R and not checkIdentity(leftAln[3], AlignedIdn):
					R = 0
					status = "PartialFailed_I"
					print "PartialR covered gap:", gapStart, "Contig Idnetity <", AlignedIdn, ":", contiggg
				
				#Can't reach gap region
				if alignedScafPos1_R > gapEnd:
					R = 0

				if R:
					if pre_contig != contiggg:
						contig_seq = getContigSeq(contiggg, contigInfo[1])
						pre_contig = contiggg
					contig_seqR = contig_seq[alignReadPos1_R - 1:alignReadPos2_R]
					contigNameR = contiggg
					break
				
                #Both sides are fillable       
		if L and R:
                        #ContigL and ContigR intersect, add Ns to represent the gap
			if alignedScafPos2_L >= alignedScafPos1_R:
				filling.append([gapStart, alignedScafPos1_L, alignedScafPos2_R, contig_seqL+"N"*30+contig_seqR, "L+30N+R", gapEnd-gapStart+1, len(contig_seqL), len(contig_seqR), alignedScafPos2_L, alignedScafPos1_R, contigNameL+"_"+contigNameR])

			else:
				filling.append([gapStart, alignedScafPos1_L, alignedScafPos2_L, contig_seqL, "L+N+R", gapEnd-gapStart+1, contigNameL])
				filling.append([gapStart, alignedScafPos1_R, alignedScafPos2_R, contig_seqR, "L+N+R", gapEnd-gapStart+1, contigNameR])
                #Only left side is fillable
		elif L:
			m = "L"
			if alignedScafPos2_L > gapEnd:
				alignedScafPos2_L = gapEnd
				contig_seqL += "N"*30 
				m = "L+30N"
			filling.append([gapStart, alignedScafPos1_L, alignedScafPos2_L, contig_seqL, m, gapEnd-gapStart+1, contigNameL])
                #Only right side is fillable
		elif R:
			m = "R"
			if alignedScafPos1_R < gapStart:
				alignedScafPos1_R = gapStart
				contig_seqR = "N"*30 + contig_seqR
				m = "30N+R"

			filling.append([gapStart, alignedScafPos1_R, alignedScafPos2_R, contig_seqR, m, gapEnd-gapStart+1, contigNameR])

		elif not filled:
			filling.append([gapStart, gapStart, gapEnd, status, "X", gapEnd-gapStart+1, "X"])

gap_record = []
scaf=''
with open(scafFILE, "r") as scaf_file:
	for line in scaf_file:
		if ">" in line:
			scaf_name = line.strip().split("|")[0]
			scaf = ''
			continue

		scaf += line.strip()

#Fill the gaps on scaffold
for i, info in enumerate(reversed(filling)):
	p1 = len(info[3])
	p2 = len(scaf[info[2]:])
	if info[4] == "L+30N+R":
		gap_record.append([info[0], info[5], p2+int(info[7])-1, p2, info[4], "Filled", info[10], info[9], info[2]])
		gap_record.append([info[0], info[5], p1+p2-1, p1+p2-int(info[6]), info[4], "Filled", info[10], info[1], info[8]])
	elif info[4] == "X":
		gap_record.append([info[0], info[5], info[5]+p2-1, p2, info[4], info[3], info[6], 0, 0])
		continue
	else:
		gap_record.append([info[0], info[5], p1+p2-1, p2, info[4], "Filled", info[6], info[1], info[2]])
	scaf = scaf[:info[1]-1] + info[3].lower() + scaf[info[2]:]

#Output gap-filled scaffold
scafFilled_file = open(filledScaf_fname, "w")
cnt = 0
scaf_len = len(scaf)
scafFilled_file.write(scaf_name + "|size" + str(scaf_len) + "\n")
while cnt < scaf_len:
	end = scaf_len if cnt + 100 > scaf_len else cnt + 100
	scafFilled_file.write(str(scaf[cnt:end]) + "\n")
	cnt += 100

#Gap-filling's record
out_file=open(out_fname, "w")
for gap_info in gap_record:
	out_file.write(str(gap_info[0])+" "+str(gap_info[1])+" "+gap_info[4]+" "+str(scaf_len - gap_info[2])+" "+str(scaf_len - gap_info[3])+" "+gap_info[5]+" "+gap_info[6]+" "+str(gap_info[7])+" "+str(gap_info[8])+"\n")


