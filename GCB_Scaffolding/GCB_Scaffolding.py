from Bio.Seq import Seq
import os
import sys

cnt = 0
connect = 0
ConnectedScaf = []
s2R = False
ScafStrand = {"R": "F", "F": "R"}
ConnectInfoLog = sys.argv[1]
wd = os.path.dirname(ConnectInfoLog) 
fastaFILE = sys.argv[2]

def getScaffold(fastaname, scafName, Reverse):
	scaf = ''
	skip = 1
	with open(fastaname, "r") as fafile:
		for l in fafile:
			x = l.split()
			if x[0][0] == ">":
				if scaf != '': break
				if x[0][1:] == scafName:
					skip = 0
					continue
				else:
					skip = 1
			if not skip:
				scaf += x[0]

	if Reverse:
		return Seq(scaf).reverse_complement()
	else:
		return scaf

def getContig(folder, contigName, Reverse, p1, p2):
	skip = 1
	c = ''
	with open(sys.argv[3], 'r') as f:
		for l in f:
			if l[0] == ">":
				name = l.strip()[1:]
				if name == contigName:
					skip = 0
					continue
				else:
					skip = 1
				if c != '': break
			if not skip:
				c += l.strip()
	seq = c[p1-1:p2]
	if Reverse:
		return Seq(seq).reverse_complement()
	else:
		return seq

def WriteFastq(scaf, fname):
	cnt = 0
	scaf_len = len(scaf)
	scafFilled_file = open(fname, "w")
	scafFilled_file.write(">" + fname.split(".fa")[0] + "\n")
	while cnt < scaf_len:
		end = scaf_len if cnt + 100 > scaf_len else cnt + 100
		scafFilled_file.write(str(scaf[cnt:end]) + "\n")
		cnt += 100

def scaffolding(s1_info, s2_info, s2R, folder, NewScafIDX):
	S = "R" if s2R else "F"
	s1R = False
	s1_scafListInfo = [[s1_info[0]+"F"], "temp.fa"]
	s2_scafListInfo = [[s2_info[0]+S], "temp.fa"]
	for scafListInfo in ConnectedScaf:
		scafList = scafListInfo[0]
		l = len(scafList)
		if s1_info[0]+"F" in scafList or s1_info[0]+"R" in scafList:
			if s1_info[0] == scafList[0][:-1]:
				#Check if is the same end
				if scafList[0][-1] == "F" and s1_info[1] != "H": return False
				if scafList[0][-1] == "R" and s1_info[1] != "T": return False
				s1_scafListInfo = scafListInfo

			elif s1_info[0] == scafList[l-1][:-1]:
				if scafList[l-1][-1] == "F" and s1_info[1] != "T": return False
				if scafList[l-1][-1] == "R" and s1_info[1] != "H": return False
				s1_scafListInfo = scafListInfo
			else:
				print s1_info[0], "is already connected by other scaffolds."
				return False

		if s2_info[0]+"F" in scafList or s2_info[0]+"R" in scafList:
			if s2_info[0] == scafList[0][:-1]:
				if scafList[0][-1] == "F":
					if s2R and s2_info[1] != "T": return False
					if not s2R and s2_info[1] != "H": return False
				else:
					if s2R and s2_info[1] != "H": return False
					if not s2R and s2_info[1] != "T": return False
				s2_scafListInfo = scafListInfo

			elif s2_info[0] == scafList[l-1][:-1]:
				if scafList[0][-1] == "F":
					if s2R and s2_info[1] != "H": return False
					if not s2R and s2_info[1] != "T": return False
				else:
					if s2R and s2_info[1] != "T": return False
					if not s2R and s2_info[1] != "H": return False
				s2_scafListInfo = scafListInfo
			else:
				print s2_info[0], "is already connected by other scaffolds."
				return False

	if s1_scafListInfo == s2_scafListInfo:
		print "Creat a cycle:", s1_info[0], s2_info[0]
		return False


	s1_scaffold = getScaffold(fastaFILE, s1_info[0], False)
	s1_scafLength = len(s1_scaffold)
	s2_scaffold = getScaffold(fastaFILE, s2_info[0], s2R)
	s2_scafLength = len(s2_scaffold)

	if s1_scafListInfo in ConnectedScaf:
		ConnectedScaf.remove(s1_scafListInfo)
		s1_scaffold = getScaffold(s1_scafListInfo[1], s1_scafListInfo[1].split(".fa")[0], False)
		#delte old fa
		os.remove(s1_scafListInfo[1])
	if s2_scafListInfo in ConnectedScaf:
		ConnectedScaf.remove(s2_scafListInfo)
		s2_scaffold = getScaffold(s2_scafListInfo[1], s2_scafListInfo[1].split(".fa")[0], False)
		os.remove(s2_scafListInfo[1])

	s1_scafList = s1_scafListInfo[0]
	s2_scafList = s2_scafListInfo[0]
	s1_contigpos = [int(s1_info[2][4][1:-1]), int(s1_info[2][5][:-1])]
	s2_contigpos = [int(s2_info[2][4][1:-1]), int(s2_info[2][5][:-1])]
	contigName = s1_info[2][0]
	if s1_info[1] > s2_info[1]:
		#s1_T---s2_H
		if s1_scafList[len(s1_scafList)-1][:-1] != s1_info[0]:
			s1_tempList = []
			for i in range(len(s1_scafList)-1, -1, -1):
				s1_tempList.append(s1_scafList[i][:-1]+ScafStrand[s1_scafList[i][-1]])
			s1_scafList = s1_tempList
			s1_scaffold = Seq(s1_scaffold).reverse_complement()
		if s2_scafList[0][:-1] != s2_info[0]:
			s2_tempList = []
			for i in range(len(s2_scafList)-1, -1, -1):
				s2_tempList.append(s2_scafList[i][:-1]+ScafStrand[s2_scafList[i][-1]])
			s2_scafList = s2_tempList
			s2_scaffold = Seq(s2_scaffold).reverse_complement()

		#concat two scaffold
		if s1_info[2][1] == '+':
			if s1_contigpos[0] > s2_contigpos[1]: return False
			contig = getContig(folder, contigName, False, s1_contigpos[0], s2_contigpos[1])
			
		else:
			if s2_contigpos[0] > s1_contigpos[1]: return False
			contig = getContig(folder, contigName, True, s2_contigpos[0], s1_contigpos[1])

		newScaf = s1_scaffold[:int(s1_info[2][2]) - 1 - s1_scafLength] + contig + s2_scaffold[int(s2_info[2][3]):]

		#write file
		outFile = wd + "/NewScaf"+str(NewScafIDX)+".fa"
		WriteFastq(newScaf, outFile)
		s1_scafList.extend(s2_scafList)
		ConnectedScaf.append([s1_scafList, outFile])
	else:
		#s2_T---s1_H
		if s2_scafList[len(s2_scafList)-1][:-1] != s2_info[0]:
			s2_tempList = []
			for i in range(len(s2_scafList)-1, -1, -1):
				s2_tempList.append(s2_scafList[i][:-1]+ScafStrand[s2_scafList[i][-1]])
			s2_scafList = s2_tempList
			s2_scaffold = Seq(s2_scaffold).reverse_complement()
		if s1_scafList[0][:-1] != s1_info[0]:
			s1_tempList = []
			for i in range(len(s1_scafList)-1, -1, -1):
				s1_tempList.append(s1_scafList[i][:-1]+ScafStrand[s1_scafList[i][-1]])
			s1_scafList = s1_tempList
			s1_scaffold = Seq(s1_scaffold).reverse_complement()

		
		if s1_info[2][1] == '+':
			if s2_contigpos[0] > s1_contigpos[1]: return False
			contig = getContig(folder, contigName, False, s2_contigpos[0], s1_contigpos[1])
		else:
			if s1_contigpos[0] > s2_contigpos[1]: return False
			contig = getContig(folder, contigName, True, s1_contigpos[0], s2_contigpos[1])

		newScaf = s2_scaffold[:int(s2_info[2][2]) - 1 - s2_scafLength] + contig + s1_scaffold[int(s1_info[2][3]):]

		outFile = wd + "/NewScaf"+str(NewScafIDX)+".fa"
		WriteFastq(newScaf, outFile)
		s2_scafList.extend(s1_scafList)
		ConnectedScaf.append([s2_scafList, outFile])

	return True



NewScafIDX = 0
if os.path.exists(wd + "/NewScaf.log"):
	with open (wd + "/NewScaf.log", 'r') as ccfile:
		w = 0
		for line in ccfile:
			x = line.split()
			l = []
			for i in range(1,len(x)):
				l.append(x[i])
			ConnectedScaf.append([l, x[0] + ".fa"])
			NewScafIDX = max(NewScafIDX, int(x[0][7:]))

	NewScafIDX += 1


with open(ConnectInfoLog, "r") as file:
	for line in file:
		x = line.split()

		if cnt == 0:
			folder = x[0]
			if len(x) > 1: s2R = True
		elif cnt % 6 == 2:
			s1 = x[0]
			p1 = x[1]
		elif cnt % 6 == 3:
			mappedInfo1 = x
		elif cnt % 6 == 4:
			s2 = x[0]
			p2 = x[1]
		elif cnt % 6 == 5:
			mappedInfo2 = x
		elif cnt % 6 == 0:
			if not connect:
				contigLen = int(mappedInfo1[0].split("size")[1])
				if int(mappedInfo1[6].split("_")[0]) < contigLen/5 and int(mappedInfo2[6].split("_")[0]) < contigLen/5:
					if mappedInfo1[1] == mappedInfo2[1]:
						if scaffolding([s1, p1, mappedInfo1], [s2, p2, mappedInfo2], s2R, wd + "/" + folder, NewScafIDX):
							connect = 1
							NewScafIDX += 1

					#else:
					#	print "Conflict direction of contig."
		elif x[0][0] == '=':
			cnt = 0
			connect = 0
			s2R = False
			continue

		cnt += 1


outFile = open(wd + "/NewScaf.log", "w")
p = ''
for slist in ConnectedScaf:
	p = slist[1].split(".fa")[0]+"\t"
	for i in slist[0]:
		p = p + i + "\t"
	outFile.write(p[:-1] + "\n")
outFile.close()	


