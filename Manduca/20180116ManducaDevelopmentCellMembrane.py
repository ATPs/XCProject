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

import pandas as pd
import os


folder = r'C:\Users\ATPs\OneDrive\Lab\works\2017DmSPSPH\20180166ManducaMassSpectDevelopMembrane'
fDevBand = os.path.join(folder,'20190209MsDevBand_ori.xlsx')
fDevGroup = os.path.join(folder,'20190209MsDevGroup_ori.xlsx')
fMemBand = os.path.join(folder,'20190209MsMemBand_ori.xlsx')
fMemGroup = os.path.join(folder,'20190209MsMemGroup_ori.xlsx')
df_MemBand = pd.read_excel(fMemBand,sheet_name='Sheet1')
df_MemGroup = pd.read_excel(fMemGroup,sheet_name='Sheet1')
df_DevBand = pd.read_excel(fDevBand,sheet_name='Sheet1')
df_DevGroup = pd.read_excel(fDevGroup,sheet_name='Sheet1')

column_naming = ['Protein', 'Annotated', 'description', 'Protein IDs']

# iBAQ for membrane
column_MemBand = [e for e in df_MemBand.columns if 'iBAQ' in e]
column_MemGroup = [e for e in df_MemGroup.columns if 'iBAQ' in e]
df_MemBand_iBAQ = df_MemBand.loc[:,column_naming + column_MemBand].copy()
df_MemBand_iBAQ = df_MemBand_iBAQ[(df_MemBand_iBAQ['iBAQ'] != 0)].copy()
df_MemGroup_iBAQ = df_MemGroup.loc[:,column_naming + column_MemGroup].copy()
df_MemGroup_iBAQ = df_MemGroup_iBAQ[(df_MemGroup_iBAQ['iBAQ'] != 0)].copy()

#Normalize so that sum of each sample have the sampe iBAQ value
samples = ['FP','FS','HP','HS']
samples_intensity = [df_MemGroup_iBAQ['iBAQ '+ e].sum() for e in samples]
average_intensity = sum(samples_intensity)/len(samples_intensity)
dc_normalize_factor = {k:average_intensity/v for k,v in zip(samples,samples_intensity)}

for s in samples:
    column = 'iBAQ '+s
    df_MemGroup_iBAQ[column] = df_MemGroup_iBAQ[column] * dc_normalize_factor[s]
df_MemGroup_iBAQ['iBAQ'] = df_MemGroup_iBAQ[column_MemGroup[1:]].sum(axis='columns')

for s in samples:
    for b in range(1,9):
        column = 'iBAQ '+s+str(b)
        df_MemBand_iBAQ[column] = df_MemBand_iBAQ[column] * dc_normalize_factor[s]
df_MemBand_iBAQ['iBAQ'] = df_MemBand_iBAQ[column_MemBand[1:]].sum(axis='columns')

df_MemBand_iBAQ.to_excel(os.path.join(folder, '20190219MsMemBand_iBAQ.xlsx'))
df_MemGroup_iBAQ.to_excel(os.path.join(folder, '20190219MsMemGroup_iBAQ.xlsx'))

# iBAQ for development
column_DevBand = [e for e in df_DevBand.columns if 'iBAQ' in e]
column_DevGroup = [e for e in df_DevGroup.columns if 'iBAQ' in e]
df_DevBand_iBAQ = df_DevBand.loc[:,column_naming + column_DevBand].copy()
df_DevBand_iBAQ = df_DevBand_iBAQ[(df_DevBand_iBAQ['iBAQ'] != 0)]
df_DevGroup_iBAQ = df_DevGroup.loc[:,column_naming + column_DevGroup].copy()
df_DevGroup_iBAQ = df_DevGroup_iBAQ[(df_DevGroup_iBAQ['iBAQ'] != 0)]

# calibrate dilution of bands
for column in column_DevBand:
    if column.endswith('3'):
        df_DevBand_iBAQ[column] = df_DevBand_iBAQ[column] * 20
        print(column)
    if 'Adult' in column:
        if column.endswith('1') or column.endswith('2') or column.endswith('7') or column.endswith('8'):
            df_DevBand_iBAQ[column] = df_DevBand_iBAQ[column] * 20
            print(column)
df_DevBand_iBAQ['iBAQ'] = df_DevBand_iBAQ[column_DevBand[1:]].sum(axis=1)

# calculate iBAQ in groups from bands
for column in column_DevGroup[1:]:
    df_DevGroup_iBAQ[column] = df_DevBand_iBAQ.loc[:,[column+str(n+1) for n in range(8)]].sum(axis=1)
    print(column)
df_DevGroup_iBAQ['iBAQ'] = df_DevGroup_iBAQ[column_DevGroup[1:]].sum(axis=1)

# Normalize so that sum of each sample have the sampe iBAQ value
samples = column_DevGroup[1:]
samples_intensity = [df_DevGroup_iBAQ[e].sum() for e in samples]
average_intensity = sum(samples_intensity)/len(samples_intensity)
dc_normalize_factor = {k:average_intensity/v for k,v in zip(samples,samples_intensity)}

for s in samples:
    column = s
    df_DevGroup_iBAQ[column] = df_DevGroup_iBAQ[column] * dc_normalize_factor[s]
df_DevGroup_iBAQ['iBAQ'] = df_DevGroup_iBAQ[column_DevGroup[1:]].sum(axis='columns')

for s in samples:
    for b in range(1,9):
        column = s+str(b)
        df_DevBand_iBAQ[column] = df_DevBand_iBAQ[column] * dc_normalize_factor[s]
df_DevBand_iBAQ['iBAQ'] = df_DevBand_iBAQ[column_DevBand[1:]].sum(axis='columns')

df_DevBand_iBAQ.to_excel(os.path.join(folder, '20190219MsDevBand_iBAQ.xlsx'))
df_DevGroup_iBAQ.to_excel(os.path.join(folder, '20190219MsDevGroup_iBAQ.xlsx'))