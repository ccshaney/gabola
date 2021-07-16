def addPairSum_atScafHT(df):
    df['ScafHT_withPairSum'] = df.apply(lambda x: x['Scaf_HT'] + '__' + str(x['pairSum']), axis=1)
    return df[['ScafHT_withPairSum', 'BX']]


def ScafHT_withPairSum_splitOff(df):
    df['ScafA'] = df['ScafA__pair'].apply(lambda x: x.split('__')[0])
    df['ScafA_pairSum'] = df['ScafA__pair'].apply(lambda x: int(x.split('__')[1]))
    df['ScafB'] = df['ScafB__pair'].apply(lambda x: x.split('__')[0])
    df['ScafB_pairSum'] = df['ScafB__pair'].apply(lambda x: int(x.split('__')[1]))
    return df[['BX', 'ScafA', 'ScafA_pairSum', 'ScafB', 'ScafB_pairSum']]


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


# 去除最小值小於 limit 的： ie, (S1-S2) = {(1,3), (2,2), (2,3), ...} 去除 rp1 剩下 {(2,2), (2,3), ...}
# Input: all rp list & limit valuse (if limit=1, then keep rp number > 1)
# Output: all rp tuple large than limit (not equal)
def keepAbove(tupleList, limit):
    returnTupleList = []
    for tuple_i in tupleList:
        if tuple_i[0] > limit:
            returnTupleList.append(tuple_i)
    return returnTupleList


# Input: all rp tuple in keep list & all combination list (input one-by-one) & max rp number
# Output: bxCnt
def countBX(rpTupleList, rpTupleAll_i, rpN):
    bxCnt = 0
    # For rpN*-rpN*:
    if rpTupleAll_i[0] == rpN & rpTupleAll_i[1] == rpN:
        for tupleRP in rpTupleList:
            if tupleRP[0] >= rpTupleAll_i[0]:
                if tupleRP[1] >= rpTupleAll_i[1]:
                    bxCnt += 1
        return bxCnt
    # For rpN-rpN*:
    elif rpTupleAll_i[0] != rpN & rpTupleAll_i[1] == rpN:
        for tupleRP in rpTupleList:
            if tupleRP[0] == rpTupleAll_i[0]:
                if tupleRP[1] >= rpTupleAll_i[1]:
                    bxCnt += 1
        return bxCnt
    # Others:
    if rpTupleAll_i in rpTupleList:
        for tupleRP in rpTupleList:
            if tupleRP[0] >= rpTupleAll_i[0]:
                if tupleRP[1] >= rpTupleAll_i[1]:
                    bxCnt += 1
        return bxCnt
    else:
        return bxCnt


# Generate columns name:
def formatColumnName(tupleRP, rpN):
    if tupleRP[0] == rpN:
        return '({}*, {}*)'.format(rpN, rpN)
    else:
        if tupleRP[1] == rpN:
            return str(tupleRP).replace(')', '*)')
        else:
            return str(tupleRP)


def rp1toN_BX1toN(df, dropOut, rpN, BXN, ScafAB_rpTuple_detail_Fpath, ScafHT_BXSupportTableForScafPair_Fpath):
    import itertools
    import pandas as pd
    # Inupt df: 'BX', 'Scaf_HT', 'pairSum'
    # Output df: 'BX', 'ScafHT_withPairSum'
    df_modifyColumn = addPairSum_atScafHT(df)
    # Generate [ BX - ScafA_withPair - ScafB_withPair ]:
    df_groupbyBX = df_modifyColumn.groupby(['BX']).agg({'ScafHT_withPairSum': list}).reset_index()
    df_groupbyBX.columns = ['BX', 'ScafHT_withPairSum_list']
    # Generate combinations:
    df_BX_ScafAB = pd.DataFrame(list_scaf_to_network_df_withBX(df_groupbyBX, 'ScafHT_withPairSum_list'),
                                columns=['BX', 'ScafA__pair', 'ScafB__pair'])
    # Generate [ ScafA - ScafB - BXCnt ]:
    df_BX_ScafAB = ScafHT_withPairSum_splitOff(df_BX_ScafAB)
    # NEW 算法：
    df_BX_ScafAB['rp_tuple'] = df_BX_ScafAB.apply(
        lambda x: (x['ScafA_pairSum'], x['ScafB_pairSum']) if x['ScafA_pairSum'] < x['ScafB_pairSum'] else (
        x['ScafB_pairSum'], x['ScafA_pairSum']), axis=1)
    ScafAB_rpTuple = df_BX_ScafAB.groupby(['ScafA', 'ScafB']).agg({'rp_tuple': list}).reset_index()
    # Remove tuple by drop limit:
    columnName = 'keep_above_{}'.format(dropOut)
    ScafAB_rpTuple[columnName] = ScafAB_rpTuple['rp_tuple'].apply(lambda x: keepAbove(x, dropOut))
    ScafAB_rpTuple['keepLen'] = ScafAB_rpTuple[columnName].apply(lambda x: len(x))
    ScafAB_rpTuple_keep = ScafAB_rpTuple.loc[ScafAB_rpTuple['keepLen'] > 0].reset_index(drop=True)
    colNameList = []
    for rpTupleAll_i in list(itertools.combinations_with_replacement(range(dropOut + 1, rpN + 1), 2)):
        # 一個一個 (1,1), (1,2), 依序算，一個 raw 算一次！
        colName = formatColumnName(rpTupleAll_i, rpN)
        colNameList.append(colName)
        print("Match {}".format(colName))
        ScafAB_rpTuple_keep[colName] = ScafAB_rpTuple_keep[columnName].apply(lambda x: countBX(x, rpTupleAll_i, rpN))
    ScafAB_rpTuple_keep.to_csv(ScafAB_rpTuple_detail_Fpath, sep='\t', index=False)
    print('[SAVE] ScafAB Read Pair Tuple detail at...{}'.format(ScafAB_rpTuple_detail_Fpath))
    # 計算 BX count:
    BXCnt_dict = {}
    for colName in colNameList:
        print('Groupby {}'.format(colName))
        eachRpN_BXCnt = ScafAB_rpTuple_keep.groupby([colName]).size()
        print(eachRpN_BXCnt)
        BXCnt_index_list = list(eachRpN_BXCnt.index)
        ansList = list(itertools.repeat(0, BXN))
        #         print("BXCnt_index_list = ", BXCnt_index_list)
        for BXCnt_i in BXCnt_index_list:
            if BXCnt_i >= BXN:
                newValue = ansList[BXN - 1] + eachRpN_BXCnt[BXCnt_i]
                ansList[BXN - 1] = newValue
            elif BXCnt_i == 0:
                pass
            else:
                ansList[BXCnt_i - 1] = eachRpN_BXCnt[BXCnt_i]
        BXCnt_dict.update({colName: ansList})
    # New df column name list:
    columnList = []
    for i in range(1, BXN + 1):
        if i == BXN:
            BXName = 'BX' + str(i) + '*'
            columnList.append(BXName)
        else:
            BXName = 'BX' + str(i)
            columnList.append(BXName)
    result_df = pd.DataFrame.from_dict(BXCnt_dict, orient='index', columns=columnList)
    print("\n{}".format(result_df))
    result_df.to_csv(ScafHT_BXSupportTableForScafPair_Fpath, sep='\t')
    print('[SAVE] ScafHT BX Support Table For ScafPair at...{}'.format(ScafHT_BXSupportTableForScafPair_Fpath))
    return result_df


### test:
# import pandas as pd
# ScafAB_rpTuple_detail_Fpath = newPathName(input_samfile_path, '_ScafHT_BXrpList_rpTupleBXCnt.tsv')
# ScafHT_BXSupportTableForScafPair_Fpath = newPathName(input_samfile_path, '_ScafHT_BXSupportTableForScafPair_byAffinity.tsv')
# remove_ScafEndAbove20 = '/mnt/biogalaxy/pyfu_bk/eelv42_nonDup_10X/v1_out/eelv42_nonDup_C70M60_ScafHeadTail_BX_pairSum_remove3BX.tsv'
# ScafHeadTail_BX_pairSum_df = pd.read_csv(remove_ScafEndAbove20, sep='\t')
# rp1toN_BX1toN(ScafHeadTail_BX_pairSum_df, 1, 10, 20, ScafAB_rpTuple_detail_Fpath, ScafHT_BXSupportTableForScafPair_Fpath)
