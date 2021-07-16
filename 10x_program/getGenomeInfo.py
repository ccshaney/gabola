import argparse
import pandas as pd
import gzip
from Bio import SeqIO
import os
import re

# #############
# ERROR class:
# #############

class InputFileNameError(Exception):
    pass

# ##########
# Function:
# ##########

# 1. For genome fa file:
#    use recursive to find next position, which is not gaps, so set ori_start, ori_end separately.
def genomeFa_to_genomeLen_POS(seq, length, head_tail, ori_start, ori_end, head_tail_range):
    import re
    N_start = 0
    N_end = 0
    minus_len = 0 # N_length
    sikp = re.compile("[Nn]+")
    final_len = 0
    # Head/Tail range over half of scaffold:
    if (head_tail == 'head') and (ori_end > (length//2)):
        print("ori_end > (length//2)")
        print("ori_start = {}, ori_end = {}".format(ori_start, ori_end))
        return (ori_start, ori_end)
    elif (head_tail == 'tail') and (ori_start < (length//2)):
        print("ori_start < (length//2)")
        print("ori_start = {}, ori_end = {}".format(ori_start, ori_end))
        return (ori_start, ori_end)
    # 1) if no gaps in head_tail_range! (if first loop is not True...it will always False.)
    if (sikp.search(seq[ori_start:ori_end])==None):
        final_len = head_tail_range - minus_len
        if head_tail == 'head':
            ori_end = head_tail_range
        # print("No any gaps(N) in {} (range:{})!".format(head_tail, head_tail_range))
    # 2) walk through gaps:
    else:
        if (head_tail=='head'):
            ### For Head: add larger range to find position without gaps (go through RIGHT direct --> ),
            # ie, ori_end - len(N_length) == head_tail_range (i.e. 20k)
            for pos in re.finditer("[Nn]+", seq[ori_start:ori_end]):
                # find all gaps and count the length of all gaps in this range
                N_start = int(pos.start())
                N_end = int(pos.end())
                minus_len += (N_end - N_start)
            ori_end += (minus_len-1) # if range of length - gaps length < head_tail_range (i.e. 20k), then iterate this function again.
        elif (head_tail=='tail'):
            ### For Tail: add larger range to find position without gaps (go through LEFT direct <-- ),
            # ie, ori_start - len(N_length) == head_tail_range (i.e. 20k)
            for pos in re.finditer("[Nn]+", seq[ori_start:ori_end]):
                # find all gaps and count the length of all gaps in this range
                N_start = int(pos.start())
                N_end = int(pos.end())
                minus_len += (N_end - N_start)
            ori_start -= minus_len # if range of length - gaps length < head_tail_range (i.e. 20k), then iterate this function again.
        final_len = len(seq[ori_start:ori_end]) - minus_len
    if (final_len == head_tail_range):
        return (ori_start, ori_end)
    else:
        ### recursive to find next position
        genomeFa_to_genomeLen_POS(seq, length, head_tail, ori_start, ori_end, head_tail_range)


def splitFileSuffixAndRemoveLatest(InputfileName):
    # InputstrReturn: only file name.
    # Return: The file name that remove original suffix (ie, .fa, .fq, etc).
    strReturn = ''
    for part in InputfileName.split(".")[:-1]:
        strReturn += part
    return strReturn

# ###############
# argparse part:
# ###############

version_number = '0.0.1'
parser = argparse.ArgumentParser(
    description='''Process genome fasta file to get Head/Tail position for each scaffold.''',
    epilog='''Any questions or comments, please contact Fu Po-Ying <spashleyfu@gmail.com> 
getGenomeInfo.py version: {}'''.format(version_number))

parser.add_argument('-f', '--fasta', dest="fasta", required=True, type=str,
                    help='Reference genome. (fasta format, named as .fa(.gz) or .fasta(.gz))')
parser.add_argument('-r', '--range', dest="range", type=int, default=20000,
                    help='Range length of each scaffold Head/Tail. [default:20000]')
parser.add_argument('-g', '--withgap', dest="withgap", type=bool, default=True,
                    help='''Boolean value to present the Head/Tail range is with or without gap.
                    True is mean the Head/Tail range with gap. 
                    If you input [-r 20000], then the Head/Tail range will calculate from both end (Head/Tail) to 20000 
                    base, and in that range may exist some gaps (N base).
                    False is mean the Head/Tail range will count without gap (N base). (default: True)''')
parser.add_argument('-q', '--quiet', dest="quiet", type=str, default=False,
                    help='Process without any output. [default:False]')

options = parser.parse_args()
input_genomeFa_path = options.fasta
head_tail_range = options.range
with_or_without = options.withgap
quiet = options.quiet

# ######################################################################################
# Process Genome FASTA then Get each Scaffold's Head/Tail Range: (with and without gap)
# ######################################################################################

if not quiet:
    print("\n### START TO PROCESS GENOME FASTA FILE...")

# Check Genome FASTA file name and Get output TSV file name:
# Get FASTA absolute path:
input_genomeFa_path = os.path.abspath(input_genomeFa_path)
if not quiet:
    print("[INPUT GENOME FASTA]\t{}".format(input_genomeFa_path))

if re.search(r"\.((fa)|(fasta))\.gz$", input_genomeFa_path):
    input_genomeFa_path_GZ = input_genomeFa_path
    input_genomeFa_path = splitFileSuffixAndRemoveLatest(input_genomeFa_path)
elif re.search(r"\.((fa)|(fasta))$", input_genomeFa_path):
    input_genomeFa_path_GZ = ''
else:
    raise InputFileNameError("Input fasta file must suffix in .fa(.gz) or .fasta(.gz), Please Check your INPUT FASTA")

outpath_FaInfo, n = re.subn(r"\.((fa)|(fasta))$",
                            "_range{}_scafLen.tsv".format(head_tail_range),
                            input_genomeFa_path)

# process output genome info path: turn to absolute path!
outpath_FaInfo = os.path.abspath(outpath_FaInfo)

# read FASTA file: (with or without gzip)
if not quiet:
    print("Processing...")
if input_genomeFa_path_GZ == '':
    genomeFA_dict = SeqIO.to_dict(SeqIO.parse(input_genomeFa_path, "fasta"))
else:
    with gzip.open(input_genomeFa_path_GZ, "rt") as handle:
        genomeFA_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

Scaf_list = []
while len(genomeFA_dict) > 0:
    lessThanRange = False
    key1, values1 = genomeFA_dict.popitem()
    scafID = str(key1)
    seq = str(values1.seq)
    scafLen = len(seq)
    if scafLen < (head_tail_range*2): # len(seq) < head_tail_range*2 (ie, 20k*2)
        lessThanRange = True
        Scaf_list.append([scafID, scafLen, 0, scafLen//2, scafLen//2,
                          (scafLen//2 + 1), (scafLen//2 + 1), scafLen, lessThanRange])
    else:
        ### Find Head/Tail position for each scaffold:
        ### 1. Head:
        head_start_end = genomeFa_to_genomeLen_POS(seq, scafLen, "head", 0, head_tail_range+1, head_tail_range)
        ### 2. Tail:
        tail_ori_start = len(seq) - head_tail_range
        tail_ori_end = len(seq)
        tail_start_end = genomeFa_to_genomeLen_POS(seq, scafLen, "tail", tail_ori_start, tail_ori_end, head_tail_range)
        Scaf_list.append([scafID, scafLen, head_start_end[0], head_tail_range, head_start_end[1],
                          tail_ori_start, tail_start_end[0], tail_start_end[1], lessThanRange])

Scaf_df = pd.DataFrame(Scaf_list, columns=['ScafID', 'ScafLen', 'Head_start', 'Head_end_withgap', 'Head_end_withoutgap',
                                           'Tail_start_withgap', 'Tail_start_withoutgap', 'Tail_end', 'lessThanRange'])
Scaf_df = Scaf_df.sort_values(['ScafLen'], ascending=False).reset_index(drop=True)

### Save FASTA info to file:
Scaf_df.to_csv(outpath_FaInfo, index=False, sep='\t')
if not quiet:
    print("SUCCESSFULLY PROCESSED.\n[SAVE] GENOME Head/Tail INFO at... {}".format(outpath_FaInfo))
