# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 13:16:46 2016
practise hierarchical Cluster
@author: k
"""
def hiearchicalClusteringStudy():
    import scipy
    import pylab
    import scipy.cluster.hierarchy as sch
    
    # Generate random features and distance matrix.
    x = scipy.rand(40)
    D = scipy.zeros([40,40])
    for i in range(40):
        for j in range(40):
            D[i,j] = abs(x[i] - x[j])
    
    import numpy as np
    D = np.random.rand(40,40)
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='centroid')
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    
    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.show()

def fillInContentStudy():
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams['pdf.fonttype'] = 42
    data = np.random.rand(5, 4)
    heatmap = plt.pcolor(data)
    
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            plt.text(x + 0.5, y + 0.5, '%.4f' % data[y, x],
                     horizontalalignment='center',
                     verticalalignment='center',fontname='Arial',color = 'w'
                     )
    
    plt.colorbar(heatmap)
    plt.savefig('hiearchical.pdf')
    plt.show()

def hiearchicalClusteringStudyTestRSEMresult20170305():
    '''
    try to make heatmap with csv files using python
    '''
    #read in csv file
    import pandas as pd
    import numpy as np
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170305AgSPSPH_RSEM_FPKM_expression340.csv",header = 1,index_col = 0)
    df2 = df.iloc[:,1:]
    
    dfo = df2.copy()
    
    
    df3 = np.log2(df2+1)
    df4 = df3.copy()
    df4[df4>10] = 10
    
    
    def value2textNotes(value):
        '''
        0 to 9, return '0' to '9'
        10 to 'A', 11 to 'B'
        '''
        s = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        v = round(value, 0)
        v = int(v)
        return s[v]
    df5 = df3.copy()
    df5 = df5.applymap(value2textNotes)
    
    
    
    
    #try to plot FPKM values directly
    import matplotlib.pyplot as plt
    plt.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots()
    fig.set_size_inches(10.5,34)
    heatmap = ax.pcolormesh(df4,cmap = 'jet')
    plt.colorbar(heatmap, orientation='horizontal',pad=0.05, shrink = 0.5)
    ax.set_frame_on(False)
    ax.set_yticks(np.arange(df2.shape[0])+0.5, minor = False)
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.set_xticks(np.arange(df2.shape[1])+0.5, minor = False)
    ax.set_yticklabels(df2.index,fontsize=4)
    ax.set_xticklabels(df2.columns,fontsize=4)
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    plt.xticks(rotation=90)
    
    for y in range(df5.shape[0]):
        for x in range(df5.shape[1]):
            plt.text(x + 0.5, y + 0.4, df5.iloc[y,x],
                     horizontalalignment='center',
                     verticalalignment='center',fontname='Arial',color = 'w',fontsize = 5, weight = 'bold'
                     )
    
    plt.savefig('heatmap.pdf')
    plt.close()

def plotWithSeabornStudy20170306():
    #read in some data
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170305AgSPSPH_RSEM_FPKM_expression340.csv",header = 1,index_col = 0)
    df2 = df.iloc[:,1:]
    
    data = df2.iloc[:30,:]
    dataPlot = np.log2(data+1)
    def value2textNotes(value):
        '''
        0 to 9, return '0' to '9'
        10 to 'A', 11 to 'B'
        '''
        s = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        v = round(value, 0)
        v = int(v)
        return s[v]
    
#    plt.rcParams['figure.figsize'] = (11, 8)
    plt.rcParams['font.size'] = 6
#    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['pdf.fonttype'] = 42
    
    fig = sns.clustermap(dataPlot, cmap = 'jet', col_cluster=False, vmin = 0, vmax = 10, square = True, annot = True, fmt = '.0f')
    
    dataFill = dataPlot.copy()
    dataFill = dataFill.iloc[fig.dendrogram_row.reordered_ind, :]
    dataFill = dataFill.applymap(value2textNotes)
    
    fig = sns.clustermap(dataPlot, cmap = 'jet', col_cluster=False, vmin = 0, vmax = 10, square = True, annot = dataFill, fmt = 's', annot_kws ={'weight':'bold'}, figsize = (12,5+0.07*dataFill.shape[0]))
    plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation = 0, fontsize = 6)
    fig.ax_heatmap.xaxis.tick_top()
    plt.setp(fig.ax_heatmap.xaxis.get_majorticklabels(), rotation = 90, fontsize = 6)
    fig.ax_row_dendrogram.autoscale()
    hm=fig.ax_heatmap.get_position()
    fig.ax_row_dendrogram.set_position([hm.x0-0.05, hm.y0, 0.05, hm.height])
    fig.savefig('cluster.pdf')
    plt.close()