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

def generate_xmax_PlusRange(maxNumb, xy):
    digitsNumb = digits_count(maxNumb, xy)
    if digitsNumb <= 2:
        return 0
    else:
        return digitsNumb - 2


def digits_count(numb, xy):
    assert numb > 0, "max number of {} axis must large than 0!".format(xy)
    i = 1
    while round(numb / 10, 2) >= 1:
        i += 1
        numb = numb // 10
    return i


def plot_with_hist(x_Series, y_Series, bins, figsize_tuple, xLabel, yLabel, saveTF, figPath):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    ### data:
    x = x_Series  # BX_ScafHTCnt_pairSum['ScafHTCnt']
    y = y_Series  # BX_ScafHTCnt_pairSum['pairSum']

    ### bins:
    xbins = bins  # 100
    ybins = bins  # 100

    # no labels:
    nullfmt = NullFormatter()

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.03

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=figsize_tuple)

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, c='c', s=5)
    axScatter.grid(True)

    # now determine nice limits by hand:
    binwidth = 1
    xmax_PlusRange = generate_xmax_PlusRange(x.max(), 'X')
    ymax_PlusRange = generate_xmax_PlusRange(y.max(), 'Y')
    xmax = np.max(np.fabs(x)) + 10 ** xmax_PlusRange  # np.fabs: 可以求整數與浮點數的絕對值
    ymax = np.max(np.fabs(y)) + 10 ** ymax_PlusRange
    xlim = (int(xmax / binwidth) + 1) * binwidth
    ylim = (int(ymax / binwidth) + 1) * binwidth

    axScatter.set_xlim((-10, xlim))
    axScatter.set_ylim((-10, ylim))

    axHistx.hist(x, bins=xbins, color='r')
    axHisty.hist(y, bins=ybins, color='b', orientation='horizontal')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    axScatter.set_xlabel(xLabel, fontsize=14)
    axScatter.set_ylabel(yLabel, fontsize=14)

    if saveTF:
        plt.savefig(figPath)
    plt.show()


def draw_BX_distribution_of_ScafEndCnt_EndPairSum(df, PROJ_table_path, PROJ_fig_Path):
    import numpy as np
    # from src.gapfill_lib.process_path import newPathName
    df_2d = df.groupby(['BX']).agg({'Scaf_HT': np.size, 'pairSum': np.sum}).reset_index()
    df_2d.columns = ['BX', 'ScafEnd_Count', 'ScafEnd_PairSum']
    # save df as Table:
    df_2d.to_csv(PROJ_table_path, index=False, sep='\t')
    print("[SAVE BX-ScafEndCnt-EndPairSum TABLE]... at {}".format(PROJ_table_path))
    # draw BX info: plot and Hist
    print("[SAVE BX INFO FIG]... at {}".format(PROJ_fig_Path))
    plot_with_hist(df_2d['ScafEnd_Count'], df_2d['ScafEnd_PairSum'], 100, (16, 16),
                   'number of Scaffold End relative with BX',
                   'sum of read pair on Scaffold End', True, PROJ_fig_Path)
