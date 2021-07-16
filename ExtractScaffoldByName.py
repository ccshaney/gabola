import sys

fa = sys.argv[1]
scafName = sys.argv[2]
skip = 1
with open(fa, "r") as file:
    for line in file:
        l = line.split()
        if l[0][0] == ">":
            if l[0] == ">" + scafName:
                skip = 0
            else:
                if not skip: break
                skip = 1

        if not skip:
            print line.strip()
