import pandas as pd
import numpy as np
import argparse
from datetime import datetime
import re
import os
from src.gapfill_lib.process_path import newPathName

# ##########
# Function:
# ##########

def remove_scafHT_by_nameList(df, removeScafHT_nameList_filePath, new_df_path):
    # Input: df [Scaf_HT, BX, pairSum] & file path of removeEnds.txt
    # Output: removed df (alread save as file)
    import pandas as pd
    print("Ori file - Number of data: {}".format(len(df)))
    scaf_list_to_remove = pd.read_csv(removeScafHT_nameList_filePath, sep='\t', header=None, names=['Scaf_HT'])
    print("Number of scaffold have to remove: {}".format(len(scaf_list_to_remove)))
    new_df = df[~df.Scaf_HT.isin(scaf_list_to_remove.Scaf_HT)]
    print("Total number of data been removed: {}".format(len(BX_scafHT_pairSum) - len(BX_scafHT_pairSum_removed)))
    new_df.to_csv(new_df_path, sep='\t', index=False)
    print("[SAVE] new data without scaffolds which have been removed... at {}".format(new_df_path))
    return new_df

def BXScafHTrpSum_to_BXScafHTCntrpSum(df):
    BX_scafHTCnt_pairSum = df.groupby(['BX']).agg({'Scaf_HT':np.size, 'pairSum':np.sum}).reset_index()
    BX_scafHTCnt_pairSum.columns = ['BX', 'ScafHTCnt', 'pairSum']
    return BX_scafHTCnt_pairSum

### Output as: BX, Scaf_A, Scaf_B (combinations, only (a,b), no (b,a))
def list_scaf_to_network_df_withBX(df, column_name):
    from itertools import combinations
    ans_list = []
    for index, row in df.iterrows():
        row_BX = row['BX']
        combinations_result = list(combinations(row[column_name], 2))
        for list_n in combinations_result:
            ans_list.append([row_BX, list_n[0], list_n[1]])
    return ans_list

# ########
# main():
# ########

version_number = '0.0.1'
parser = argparse.ArgumentParser(description='''Input original ScafEnd-BX-PairSum tsv, then, 
Output new list of ScafEnd pair (sort by BX support). If input -l is not na, program will remove scaffold end
by removeEnds.txt.''',
                                 epilog='''Any questions or comments, please contact Fu Po-Ying <spashleyfu@gmail.com> 
process_sam_to_HeadTailInfo.py version: {}'''.format(version_number))
parser.add_argument('-f', '--ori_tsv', dest="ori_tsv", required=True, type=str,
                    help='path of original ScafEnd-BX-PairSum tsv file. (tsv format)')
parser.add_argument('-l', '--removeList', dest="removeEndList", default='na', type=str,
                    help='''If you want to remove some scaffold end, you need to input path of removeEnds.txt file.
                    (This file is output From Inter-Scaffold Program) [default: "na"]''')
parser.add_argument('-n', '--scafCntN', dest="scafCnt_number", type=int, default=2,
                    help='''Setting the number of scaffold end attached by one barcode. [default: 2]''')


options = parser.parse_args()
removeEndList_path = options.removeEndList
ori_tsv_path = options.ori_tsv
scafCnt_number = options.scafCnt_number

tStart = datetime.now()

if removeEndList_path != 'na':
    assert os.path.exists(removeEndList_path), "Your Input file ({}) is not exists! Please use absolute path.".format(removeEndList_path)
assert os.path.exists(ori_tsv_path), "Your Input file ({}) is not exists! Please use absolute path.".format(ori_tsv_path)

# Ori tsv file:
BX_scafHT_pairSum = pd.read_csv(ori_tsv_path, sep='\t')
print("Ori file - Number of data: {}".format(len(BX_scafHT_pairSum)))

# ScafHT list have to removed:
if removeEndList_path != 'na':
    scaf_list_to_remove = pd.read_csv(removeEndList_path, sep='\t', header=None, names=['Scaf_HT'])
    print("Number of scaffold have to remove: {}".format(len(scaf_list_to_remove)))
    # A[~A.key.isin(B.key)]: Keep Scaf_HT ONLY in A, delete Scaf_HT present in B
    BX_scafHT_pairSum_removed = BX_scafHT_pairSum[~BX_scafHT_pairSum.Scaf_HT.isin(scaf_list_to_remove.Scaf_HT)]
    minusVal = len(BX_scafHT_pairSum) - len(BX_scafHT_pairSum_removed)
    minusPercent = round((minusVal/len(BX_scafHT_pairSum))*100, 2)
    print("Total number of data been removed: {} ({}%)".format(minusVal, minusPercent))
    # Calculate ScafCnt: BXScafHTrpSum_to_BXScafHTCntrpSum()
    BX_scafHTCnt_pairSum = BXScafHTrpSum_to_BXScafHTCntrpSum(BX_scafHT_pairSum_removed)
    # Get ScafCnt == 2: get_BXList_which_scafCnt_is_two()
    ScafCnt_is2_BX_df = BX_scafHTCnt_pairSum.loc[(BX_scafHTCnt_pairSum['ScafHTCnt'] == 2)]['BX']
    ScafCnt_2_len = len(ScafCnt_is2_BX_df)
    ScafCnt2_percent_of_removed = round((ScafCnt_2_len / len(BX_scafHT_pairSum_removed)) * 100, 2)
    print("\nNumber of BX (ScafCnt==2): {} ({}%)".format(ScafCnt_2_len, ScafCnt2_percent_of_removed))
    # select data by ScafCnt_is2_BX_df:
    BX_scafHT_pairSum_removed_ScafHT_is2 = BX_scafHT_pairSum_removed[
        BX_scafHT_pairSum_removed.BX.isin(ScafCnt_is2_BX_df)]
    print(BX_scafHT_pairSum_removed_ScafHT_is2.head())
else:
    # Calculate ScafCnt: BXScafHTrpSum_to_BXScafHTCntrpSum()
    BX_scafHTCnt_pairSum = BXScafHTrpSum_to_BXScafHTCntrpSum(BX_scafHT_pairSum)
    # Get ScafCnt == 2: get_BXList_which_scafCnt_is_two()
    ScafCnt_is2_BX_df = BX_scafHTCnt_pairSum.loc[(BX_scafHTCnt_pairSum['ScafHTCnt'] == 2)]['BX']
    ScafCnt_2_len = len(ScafCnt_is2_BX_df)
    ScafCnt2_percent_of_removed = round((ScafCnt_2_len / len(BX_scafHT_pairSum)) * 100, 2)
    print("\nNumber of BX (ScafCnt==2): {} ({}%)".format(ScafCnt_2_len, ScafCnt2_percent_of_removed))
    # select data by ScafCnt_is2_BX_df:
    BX_scafHT_pairSum_removed_ScafHT_is2 = BX_scafHT_pairSum[BX_scafHT_pairSum.BX.isin(ScafCnt_is2_BX_df)]
    print(BX_scafHT_pairSum_removed_ScafHT_is2.head())

# get BX list of ScafHT:
ScafCnt2_GroupbyBX_ListScafEnd = BX_scafHT_pairSum_removed_ScafHT_is2.groupby(['BX']).agg({'Scaf_HT': list,
                                                                                           'pairSum': np.size}).reset_index()
ScafCnt2_GroupbyBX_ListScafEnd.columns = ['BX', 'ScafHT_LIST', 'ScafEnd_cnt']

# combination:
combination_df_of_ScafCnt2 = pd.DataFrame(list_scaf_to_network_df_withBX(ScafCnt2_GroupbyBX_ListScafEnd, 'ScafHT_LIST'),
                                          columns=['BX', 'ScafA', 'ScafB'])

# If ScafA's Scaf is equal to ScafB's Scaf? --> same scaffold's head/tail (removed!)
combination_df_of_ScafCnt2['SameScaf'] = combination_df_of_ScafCnt2.apply(lambda x: 1 if re.sub(r"_((Tail)|(Head))$", '', x['ScafA']) == re.sub(r"_((Tail)|(Head))$", '', x['ScafB']) else 0, axis=1)
allScafCnt2_Scaf_Len = len(combination_df_of_ScafCnt2)

### Only keep SameScaf == 0:
print("\nDelete some data which ScafA and ScafB at same scaffold...")
combination_df_of_ScafCnt2_notSameScaf = combination_df_of_ScafCnt2[combination_df_of_ScafCnt2['SameScaf']==0].reset_index(drop=True)
notSameScafLen = len(combination_df_of_ScafCnt2_notSameScaf) # 297225
print("""Number of pair count (ScafEnd==2): {}
Number of pair count (ScafEnd==2 & not at same scaffold): {}""".format(allScafCnt2_Scaf_Len, notSameScafLen))

# Output inter-scaffold list: (sort by BX support...)
ScafAB_BXlist_BXcnt = combination_df_of_ScafCnt2_notSameScaf.groupby(['ScafA', 'ScafB']).agg({'BX': list, 'SameScaf': np.size}).reset_index()
ScafAB_BXlist_BXcnt.columns = ['ScafA', 'ScafB', 'BXList', 'BXCnt']
ScafAB_BXlist_BXcnt = ScafAB_BXlist_BXcnt.sort_values(['BXCnt'], ascending=False).reset_index(drop=True)

# Save:
if removeEndList_path != 'na':
    savePath = newPathName(ori_tsv_path, '_ScafCntIs2_deleteEnds.tsv')
    ScafAB_BXlist_BXcnt.to_csv(savePath, sep='\t', index=False)
else:
    savePath = newPathName(ori_tsv_path, '_ScafCntIs2.tsv')
    ScafAB_BXlist_BXcnt.to_csv(savePath, sep='\t', index=False)
print("\n[SAVE] New list of ScafA-ScafB pair for inter-scaffold... save at {}".format(savePath))
assert os.path.exists(savePath), "Your output tsv file ({}) is not write out success! Please Check it.".format(savePath)

# print info:
combinationCnt = len(ScafAB_BXlist_BXcnt)
combinationCnt_largeThanONE = len(ScafAB_BXlist_BXcnt.loc[ScafAB_BXlist_BXcnt['BXCnt']>1])
print("""\nAfter remove ScafEnd of inter-scaffold filled...

Total combination of ScafCnt==2 is {},
Move out BX support==1, number of combination still have {} pair,
Top 5 combination is:\n{}""".format(combinationCnt, combinationCnt_largeThanONE, ScafAB_BXlist_BXcnt.head()))

print("\n[DONE] This program spent {}".format((datetime.now()-tStart)))