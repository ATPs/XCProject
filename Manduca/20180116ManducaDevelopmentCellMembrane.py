# -*- coding: utf-8 -*-
"""
Created on Wed May 23 19:30:09 2018

@author: ATPs
"""

# 20180523 check iBAQ and intensity
def check_iBAQ_intensity_relation20180523():
    # read in files
    import pandas as pd
    file = r"F:\Insects\ManducaSexta\1520\txt_group\proteinGroups_group_develop.txt"
    df = pd.read_csv(file,sep='\t')
    # extract iBAQ and intensity
    df_iBAQ = df.iloc[:,100:113]
    df_intensity = df.iloc[:,87:100]
    df_ratio = df_iBAQ.copy()
    for i in range(13):
        df_ratio.iloc[:,i] = df_iBAQ.iloc[:,i] / df_intensity.iloc[:,i]
    df_ratioGood = df_ratio.dropna()
    df_ratioGood.mean(axis=1)
    df_ratioGood.corr()
    

def annotate_identified_proteins20180528():
    # clean the data of contamination and rev_proteins
    import pandas as pd
    fDevGroup = r"F:\Insects\ManducaSexta\1520\txt_group_20times\proteinGroups.txt"
    dfDevGroup = pd.read_csv(fDevGroup, sep = "\t", low_memory=False)
    f_filter = lambda x: x[0][:5] != "CON__" and x[0][:5] != "REV__"
    dfDevGroup_filter1 = dfDevGroup.iloc[:,:(dfDevGroup.shape[1] - 13)]
    dfDevGroup_filter2 = dfDevGroup_filter1[dfDevGroup_filter1.apply(f_filter, axis=1)]
    dfDevGroup_filter2.to_csv(r"C:\Users\ATPs\OneDrive\Lab\works\2017DmSPSPH\20180166ManducaMassSpectDevelopMembrane\20180118MsDevGroup20Adjusted_filter2.csv")
    
    # get all major protein ids to a list and try to annotate them
    proteins = [e.split(';') for e in list(dfDevGroup_filter2['Majority protein IDs'])]
    proteins_all = [e2 for e1 in proteins for e2 in e1]
    
    
    filename = r"C:\Users\ATPs\OneDrive\Lab\works\2017DmSPSPH\20180166ManducaMassSpectDevelopMembrane\20180118MsDevGroup_filter2.csv"
    df = pd.read_csv(filename,index_col=0)
    proteins = list(df.iloc[:,0])
    #only keep the first protein id for annotation
    proteins1 = [e.split(';')[0] for e in proteins]
    df['gene']
    
def summarize_identified_proteins20180528():
    pass

