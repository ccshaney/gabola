##################
## MAIN Function:
##################
def rp1to3_X_BX1to3(df, savePrefixF):
    import pandas as pd
    df_modifyColumn = addPairSum_atScafHT(df)
    # Generate [ BX - ScafA_withPair - ScafB_withPair ]:
    df_groupbyBX = df_modifyColumn.groupby(['BX']).agg({'ScafHT_withPairSum': list}).reset_index()
    df_groupbyBX.columns = ['BX','ScafHT_withPairSum_list']
    df_BX_ScafAB = pd.DataFrame(list_scaf_to_network_df_withBX(df_groupbyBX, 'ScafHT_withPairSum_list'),
                                columns=['BX', 'ScafA__pair', 'ScafB__pair'])
    # Generate [ ScafA - ScafB - BXCnt ]:
    df_BX_ScafAB = ScafHT_withPairSum_splitOff(df_BX_ScafAB)
    df_BX_ScafAB_rpN_dict = select_pairSum(df_BX_ScafAB)
    result_df = sep_each_rp_by_BX1to3(df_BX_ScafAB_rpN_dict, savePrefixF)
    return result_df

def addPairSum_atScafHT(df):
    df['ScafHT_withPairSum'] = df.apply(lambda x: x['Scaf_HT'] + '__' + str(x['pairSum']), axis=1)
    return df[['ScafHT_withPairSum', 'BX']]

### To BX, Scaf_A, Scaf_B (combinations, no repeat: just a,b not b,a)
def list_scaf_to_network_df_withBX(df, column_name):
    from itertools import combinations
    ans_list = []
    for index, row in df.iterrows():
        row_BX = row['BX']
        combinations_result = list(combinations(row[column_name], 2))
        for list_n in combinations_result:
            ans_list.append([row_BX,list_n[0],list_n[1]])
    return ans_list

def ScafHT_withPairSum_splitOff(df):
    df['ScafA'] = df['ScafA__pair'].apply(lambda x: x.split('__')[0])
    df['ScafA_pairSum'] = df['ScafA__pair'].apply(lambda x: int(x.split('__')[1]))
    df['ScafB'] = df['ScafB__pair'].apply(lambda x: x.split('__')[0])
    df['ScafB_pairSum'] = df['ScafB__pair'].apply(lambda x: int(x.split('__')[1]))
    return df[['BX', 'ScafA', 'ScafA_pairSum', 'ScafB', 'ScafB_pairSum']]

def select_pairSum(df):
    print("len(df) = {}".format(len(df)))
    # 1. Sep df into ScafA - rp1~3
    df_ScafA_rp1 = df.loc[df['ScafA_pairSum'] == 1]
    df_ScafA_rp2 = df.loc[df['ScafA_pairSum'] == 2]
    df_ScafA_rp3 = df.loc[df['ScafA_pairSum'] >= 3]
    print("len(df_ScafA_rp1) = {}".format(len(df_ScafA_rp1)))
    print("len(df_ScafA_rp2) = {}".format(len(df_ScafA_rp2)))
    print("len(df_ScafA_rp3) = {}".format(len(df_ScafA_rp3)))
    # 2. Sep ScafA_rp1~3 into 3 df by ScafB rp1~3:
    # rp1:
    df_ScafA_rp1_rp1 = df_ScafA_rp1.loc[df_ScafA_rp1['ScafB_pairSum'] == 1]
    df_ScafA_rp1_rp2 = df_ScafA_rp1.loc[df_ScafA_rp1['ScafB_pairSum'] == 2]
    df_ScafA_rp1_rp3 = df_ScafA_rp1.loc[df_ScafA_rp1['ScafB_pairSum'] >= 3]
    # rp2:
    df_ScafA_rp2_rp1 = df_ScafA_rp2.loc[df_ScafA_rp2['ScafB_pairSum'] == 1]
    df_ScafA_rp2_rp2 = df_ScafA_rp2.loc[df_ScafA_rp2['ScafB_pairSum'] == 2]
    df_ScafA_rp2_rp3 = df_ScafA_rp2.loc[df_ScafA_rp2['ScafB_pairSum'] >= 3]
    # rp3:
    df_ScafA_rp3_rp1 = df_ScafA_rp3.loc[df_ScafA_rp3['ScafB_pairSum'] == 1]
    df_ScafA_rp3_rp2 = df_ScafA_rp3.loc[df_ScafA_rp3['ScafB_pairSum'] == 2]
    df_ScafA_rp3_rp3 = df_ScafA_rp3.loc[df_ScafA_rp3['ScafB_pairSum'] >= 3]
    # 3. Calculate rp1_1, rp1_2, rp1_3, rp2_2, rp2_3, rp3_3
    rp1_1 = len(df_ScafA_rp1_rp1)
    rp1_2 = len(df_ScafA_rp1_rp2) + len(df_ScafA_rp2_rp1)
    rp1_3 = len(df_ScafA_rp1_rp3) + len(df_ScafA_rp3_rp1)
    rp2_2 = len(df_ScafA_rp2_rp2)
    rp2_3 = len(df_ScafA_rp2_rp3) + len(df_ScafA_rp3_rp2)
    rp3_3 = len(df_ScafA_rp3_rp3)
    print("rp1_1 = {}".format(rp1_1))
    print("rp1_2 = {}".format(rp1_2))
    print("rp1_3 = {}".format(rp1_3))
    print("rp2_2 = {}".format(rp2_2))
    print("rp2_3 = {}".format(rp2_3))
    print("rp3_3 = {}".format(rp3_3))
    # 4. return 6 DF in dict:
    return {'df_rp1_1': df_ScafA_rp1_rp1,
            'df_rp1_2': df_ScafA_rp1_rp2.append(df_ScafA_rp2_rp1, ignore_index=True, sort=False),
            'df_rp1_3': df_ScafA_rp1_rp3.append(df_ScafA_rp3_rp1, ignore_index=True, sort=False),
            'df_rp2_2': df_ScafA_rp2_rp2,
            'df_rp2_3': df_ScafA_rp2_rp3.append(df_ScafA_rp3_rp2, ignore_index=True, sort=False),
            'df_rp3_3': df_ScafA_rp3_rp3}

def sep_each_rp_by_BX1to3(dict_rp_df, savePrefixF):
    import pandas as pd
    result_list = []
    # rp1_1:
    rp1_1_DF = dict_rp_df['df_rp1_1']
    rp1_1_BX = generate_BX1to3(rp1_1_DF, savePrefixF, '1_1')
    result_list.append(rp1_1_BX)
    # rp1_2:
    rp1_2_DF = dict_rp_df['df_rp1_2']
    rp1_2_BX = generate_BX1to3(rp1_2_DF, savePrefixF, '1_2')
    result_list.append(rp1_2_BX)
    # rp1_3:
    rp1_3_DF = dict_rp_df['df_rp1_3']
    rp1_3_BX = generate_BX1to3(rp1_3_DF, savePrefixF, '1_3')
    result_list.append(rp1_3_BX)
    # rp2_2:
    rp2_2_DF = dict_rp_df['df_rp2_2']
    rp2_2_BX = generate_BX1to3(rp2_2_DF, savePrefixF, '2_2')
    result_list.append(rp2_2_BX)
    # rp2_3:
    rp2_3_DF = dict_rp_df['df_rp2_3']
    rp2_3_BX = generate_BX1to3(rp2_3_DF, savePrefixF, '2_3')
    result_list.append(rp2_3_BX)
    # rp3_3:
    rp3_3_DF = dict_rp_df['df_rp3_3']
    rp3_3_BX = generate_BX1to3(rp3_3_DF, savePrefixF, '3_3')
    result_list.append(rp3_3_BX)
    # index:
    index_list = ['rp1_1','rp1_2','rp1_3','rp2_2','rp2_3','rp3_3']
    return pd.DataFrame(result_list, index=index_list, columns=['BX1','BX2','BX3'])

def generate_BX1to3(df, savePrefixF, rpN):
    import numpy as np
    df_ScafAB_BXCnt = df.groupby(['ScafA', 'ScafB']).agg({'BX':np.size}).reset_index()
    df_ScafAB_BXCnt.columns = ['ScafA', 'ScafB', 'BXCnt']
    df_ScafAB_BXCnt = df_ScafAB_BXCnt.sort_values(['BXCnt','ScafA'], ascending=False).reset_index(drop=True)
    # BX1:
    save_RPxBX_tsv(df_ScafAB_BXCnt, savePrefixF, rpN, 1)
    BX1 = len(df_ScafAB_BXCnt)
    # BX2:
    BX2_df = df_ScafAB_BXCnt.loc[df_ScafAB_BXCnt['BXCnt'] == 2]
    save_RPxBX_tsv(BX2_df, savePrefixF, rpN, 2)
    BX2 = len(BX2_df)
    # BX3:
    BX3_df = df_ScafAB_BXCnt.loc[df_ScafAB_BXCnt['BXCnt'] >= 3]
    save_RPxBX_tsv(BX3_df, savePrefixF, rpN, 3)
    BX3 = len(BX3_df)
    return [BX1, BX2, BX3]

def save_RPxBX_tsv(df, savePrefixF, rpN, bxN):
    save_path = newPathName(savePrefixF, '_rp{}_BX{}_ScafHTAB_BXCnt.tsv'.format(rpN, bxN))
    df.to_csv(save_path, index=False, sep='\t')
    print("[SAVED] rp{} BX{}: ScafA-ScafB-BXCnt at...{}".format(rpN, bxN, save_path))

def newPathName(path, suffix):
    import re
    returnPath = ''
    strList = path.split(".")
    if len(strList) < 2:
        # is mean no suffix:
        returnPath = strList[0] + suffix
    else:
        # cut off last . part:
        for strN in strList[:-1]:
            returnPath += (strN + ".")
        returnPath = re.sub("\.$","",returnPath)
        returnPath += suffix
    return returnPath

### test:
# df = rp1to3_X_BX1to3(ScafHeadTail_BX_pairSum_df, '/home/pyfu/eel_data/bbbTest.sam')
# print(df)
