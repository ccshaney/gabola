import argparse
import re
import os

# #########
# Fuction:
# #########

# 1. check CIGAR values:
def is_cigar_largethanN(CIAGR,seq_len, N):
    # input: CIGAR string, seq_len INT, N INT
    # output: True / False
    m_sum = 0
    if CIAGR != '*':
        # For M: "M=X"皆算，但此 dataset 沒有 = & X
        m_list = re.findall(r"\d+[M=X]", CIAGR)
        for m in m_list:
            m_num = int(re.split(r'[M=X]', m)[0])
            m_sum += m_num
        # If ratio_m >= N%, Return True:
        if (m_sum/seq_len) >= round(N/100, 2):
            return True
        else:
            return False
    else:
        return False

# #########
# Program:
# #########

version_number = '0.0.1'
parser = argparse.ArgumentParser(description='''Filter sam file by CIGAR match percentage and MAPQ score.''',
                                 epilog='''Any questions or comments, please contact Fu Po-Ying <spashleyfu@gmail.com> 
filter_CIGAR_MAPQ.py version: {}'''.format(version_number))
parser.add_argument('-s', '--sam', dest="sam", required=True, type=str,
                    help='SAM file path. (BWA MEM output format)')
parser.add_argument('-o', '--out', dest="out", required=True, type=str,
                    help='Output sam file path.')
parser.add_argument('-c', '--cigar', dest="cigar_percentage", default=70, type=int,
                    help='''CIGAR match percentage: count([M=X])/len(seq)*100. 
                    From 0 to 100. [default: 70]''')
parser.add_argument('-m', '--mapq', dest="mapq", type=int, default=60,
                    help='''MAPQ score. From 0 to 60. [default: 60]''')

options = parser.parse_args()
input_samfile = options.sam
out_sam_path = options.out
cigar_in = int(options.cigar_percentage)
mapQ_in = int(options.mapq)

# sam file path:
input_samfile = os.path.abspath(input_samfile)

# generate OUTPUT sam file path:
out_sam_path = os.path.abspath(out_sam_path)
print("out_sam_path = {}".format(out_sam_path))

########## v1:

write_out_list = []
i = 0
readID_befor = ''
line_befor = ''

with open(input_samfile, 'r') as f:
    for line in iter(f.readline, ''):
        line_list = line.split('	')
        readID_next = line_list[0]
        MAPQ = int(line_list[4])
        CIAGR = line_list[5]
        seq_len = len(line_list[9])
        # First, check MAPQ score:
        if MAPQ >= mapQ_in:
            # Second, check CIGAR match%:
            if is_cigar_largethanN(CIAGR, seq_len, cigar_in):
                # Third, make sure read pair pass both conditions:
                if readID_befor == readID_next:
                    write_out_list.append(line_befor)
                    write_out_list.append(line)
                    line_befor = ''
                    readID_befor = ''
                    i += 1
                else:
                    readID_befor = readID_next
                    line_befor = line

        if (i%10000 == 0):
            f = open(out_sam_path, 'a')
            for x in write_out_list:
                f.write(x)
            f.close()
            write_out_list = []
    f = open(out_sam_path, 'a')
    for x in write_out_list:
        f.write(x)
    f.close()
    print("[CIGAR{} MAPQ{} TOTAL READ PAIRS]\t{}".format(cigar_in, mapQ_in, i))

##########
