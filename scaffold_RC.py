from Bio.Seq import Seq
import sys

fname = sys.argv[1]
scaf = ''
with open(fname, "r") as f:
	for l in f:
		if ">" in l:
			contigName = l.strip()
		else:
			scaf += l.strip()

cnt = 0
scaf = Seq(scaf).reverse_complement()
scaf_len = len(scaf)
print contigName
while cnt < scaf_len:
	end = scaf_len if cnt + 100 > scaf_len else cnt + 100
	print scaf[cnt:end]
	cnt += 100


