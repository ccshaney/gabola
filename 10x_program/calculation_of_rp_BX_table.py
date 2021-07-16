from src.gapfill_lib.calculate_combinations import rp1toN_BX1toN
from src.gapfill_lib.process_path import newPathName
import argparse
import os
import pandas as pd

# ###ARGS:
version_number = '0.0.1'
parser = argparse.ArgumentParser(description='''Process genome fasta file to get Head/Tail position for each scaffold.
And use this scaffold Head/Tail info with input sam file to process each reads' position information.''',
                                 epilog='''Any questions or comments, please contact Fu Po-Ying <spashleyfu@gmail.com> 
process_sam_to_HeadTailInfo.py version: {}'''.format(version_number))

# required:
parser.add_argument('-f', '--tsv_file', dest="tsv_file_path", required=True, type=str,
                    help='''Path of tsv file which has information about [BX, scaffoldEnd, read pair number], 
                    ie. [PROJ_NAME]_ScafHeadTail_BX_pairSum.tsv.''')
# Set parameter - for calculate combination:
parser.add_argument('--min_rp', dest="threshold_of_BX_rp_per_end", type=int, default=0,
                    help='''Threshold of BX read pair/per end. The certified barcode control value. 
                    If value is 5 means barcode mapping to genome must support by 6 read pair. [default:0]''')
parser.add_argument('--rpN', dest="establish_combination_of_read_pair", type=int, default=10,
                    help='''Establish Number of barcode read pair at each end of scaffold end.
                    If value is 10 means this program will generate combination from (1,1) to (10*,10*). [default:10]''')
parser.add_argument('--BXN', dest="establish_number_of_bx", type=int, default=20,
                    help='''Establish Number of barcode support for pair of scaffold end.
                    If value is 20 means this program will generate BX1 to BX20*. [default:20]''')

options = parser.parse_args()
input_tsvfile_path = options.tsv_file_path
ScafEndCnt_limit = options.scafendcnt
dropOut = options.threshold_of_BX_rp_per_end
rpN = options.establish_combination_of_read_pair
BXNumb = options.establish_number_of_bx

assert os.path.exists(input_tsvfile_path), "Your Input tsv file path file is not exists! Please Check it."
assert dropOut >= 0, "scafendcnt must large than zero."
assert rpN >= 1, "rpN must large than one."
assert BXNumb >= 1, "BXNumb must large than one"

ScafAB_rpTuple_detail_Fpath = newPathName(input_tsvfile_path, '_ScafHT_BXrpList_rpTupleBXCnt.tsv')
ScafHT_BXSupportTableForScafPair_Fpath = newPathName(input_tsvfile_path, '_ScafHT_BXSupportTableForScafPair_byAffinity.tsv')
ScafHeadTail_BX_pairSum_df = pd.read_csv(input_tsvfile_path, sep='\t')
print(rp1toN_BX1toN(ScafHeadTail_BX_pairSum_df, dropOut, rpN, BXNumb, ScafAB_rpTuple_detail_Fpath, ScafHT_BXSupportTableForScafPair_Fpath))
print(["\n[PROCESS SUCCESS]"])
