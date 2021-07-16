import sys

#Get gap position for a single scaffold
FAfile = sys.argv[1]
with open(FAfile, "r") as f:
	pos = 1
	flag = 0

	for line in f:
		if ">" in line:
			continue

		l = line.strip().lower()

		for i in range(len(l)):
			if not flag and l[i] == 'n':
				flag = 1
				pos1 = pos

			if flag and l[i] != 'n':
				print pos1, pos-1, pos-pos1
				flag = 0

			pos += 1

if flag:
	print pos1, pos-1, pos-pos1
