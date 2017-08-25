# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 21:48:46 2017

@author: k
"""

def value2textNotes(value):
    '''
    0 to 9, return '0' to '9'
    10 to 'A', 11 to 'B'
    '''
    s = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    v = round(value, 0)
    v = int(v)
    if v>=36:
        print('value is greater than 36, the maximum value allowed!')
        return None
    return s[v]

def fpkmHierarchicalCluster(fileCSV,outfileHead='',col_cluster=False, row_cluster=True):
    '''
    given a csv file, plot one pdf figure with values in cells, one pdf with no values in cells,
    return a csv file with values for annotation
    fileCSV is a string or dataframe
    '''
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    if type(fileCSV) is str:
        data = pd.read_csv(fileCSV,index_col = 0)
    elif type(fileCSV) is pd.core.frame.DataFrame:
        data = fileCSV
    else:
        print('input is not filename, not DataFrame')
        return None
    
    dataPlot = np.log2(data+1)
    
#    plt.rcParams['figure.figsize'] = (11, 8)
    plt.rcParams['font.size'] = 6
#    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['pdf.fonttype'] = 42
    
    fig = sns.clustermap(dataPlot, cmap = 'jet', col_cluster=col_cluster, row_cluster = row_cluster, vmin = 0, vmax = 10, square = True, figsize = (12,8+0.08*data.shape[0]))
    plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation = 0, fontsize = 6)
    fig.ax_heatmap.xaxis.tick_top()
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), rotation = 90, fontsize = 6)
    fig.savefig(outfileHead+'HierarchicalClusterFig.pdf')

    dataFill = dataPlot.copy()
    if row_cluster:
        dataFill = dataFill.iloc[fig.dendrogram_row.reordered_ind, :]
    dataFill = dataFill.applymap(value2textNotes)
    plt.close()
    fig = sns.clustermap(dataPlot, cmap = 'jet', col_cluster=col_cluster, row_cluster = row_cluster, vmin = 0, vmax = 10, square = True, annot = dataFill, fmt = 's', annot_kws ={'weight':'bold'}, figsize = (12,5+0.08*data.shape[0]))
    plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation = 0, fontsize = 6)
    fig.ax_heatmap.xaxis.tick_top()
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), rotation = 90, fontsize = 6)
    hm=fig.ax_heatmap.get_position()
    fig.ax_row_dendrogram.set_position([hm.x0-0.05, hm.y0, 0.05, hm.height])
    fig.savefig(outfileHead+'HierarchicalClusterAnnotation.pdf')
    plt.close()
    dataFill.to_csv(outfileHead+'HierarchicalClusterAnnotation.csv')
    print('Done!')

