# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:17:25 2017
Ag proteome 2017 all data
@author: ATPs
"""

def cleanTheProteinGroupsFiles():
    '''
    proteinGroups.txt, clean the file.
    first, for unknown reason, the "protein ID" column is really ugly. clean "protein ID"
    only keep LFQ columns and iBAQ columns, currently.
    remove matched reverse sequences
    '''
    folder = 'X:\\Insects\\Anopheles_gambiae\\2017AgAllProteomeDataAnalysis\\'
    import pandas as pd
    df_group = pd.read_csv(folder+'proteinGroups_group.txt',sep='\t')
    df_band = pd.read_csv(folder+'proteinGroups_band.txt',sep='\t')
    
    #remove reverse matched lines
    df_group1 = df_group[- df_group.loc[:,'Protein IDs'].str.startswith('REV__')]
    df_band1 = df_band[- df_band.loc[:,'Protein IDs'].str.startswith('REV__')]
    print('in group', df_group.shape, df_group1.shape)
    print('in band', df_band.shape, df_band1.shape)
    
    #keep only LFQ and iBAQ columns
    df_group2 = df_group1.loc[:,['Protein IDs', 'Majority protein IDs', 'Peptide counts (all)',
       'Peptide counts (razor+unique)', 'Peptide counts (unique)',
       'Fasta headers', 'Number of proteins', 'Peptides',
       'Razor + unique peptides', 'Unique peptides']+[e for e in df_group1.columns[10:] if 'Unique peptides' in e or 'iBAQ' in e or 'LFQ' in e]]
    df_band2 = df_band1.loc[:,['Protein IDs', 'Majority protein IDs', 'Peptide counts (all)',
       'Peptide counts (razor+unique)', 'Peptide counts (unique)',
       'Fasta headers', 'Number of proteins', 'Peptides',
       'Razor + unique peptides', 'Unique peptides']+[e for e in df_band1.columns[10:] if 'Unique peptides' in e or 'iBAQ' in e or 'LFQ' in e]]
    print('in group', df_group1.shape, df_group2.shape)
    print('in band', df_band1.shape, df_band2.shape)
    
    #clean names
    ##define function
    def cleanName(s):
        '''
        s is a string in 'Protein IDs' or 'Majority protein IDs' column. return a new string
        '''
        names = s.split(';')
        newnames = ''
        for name in names:
            if len(name)>13:
                if name[:3] == 'MCD':
                    newnames = newnames + name.split('Len:')[0]+';'
                if name[:4] == 'AGAP':
                    newnames = newnames + name[:13] + ';'
            else:
                newnames = newnames + name+';'
        if newnames[-1] == ';':
            newnames = newnames[:-1]
        return newnames
    ##calculate new columns
    df_group2.loc[:,'Protein IDs'] = df_group2.loc[:,'Protein IDs'].apply(cleanName)
    df_band2.loc[:,'Protein IDs'] = df_band2.loc[:,'Protein IDs'].apply(cleanName)
    df_group2.loc[:,'Majority protein IDs'] = df_group2.loc[:,'Majority protein IDs'].apply(cleanName)
    df_band2.loc[:,'Majority protein IDs'] = df_band2.loc[:,'Majority protein IDs'].apply(cleanName)
    
    #save files
    df_group2.to_excel(folder +'20170916AgProteomeAll_group.xlsx')
    df_band2.to_excel(folder+'20170916AgProteomeAll_band.xlsx')
    
    #extract SPSPH and PPO data from df_group2
    names = df_group2.loc[:,'Protein IDs']
    newnames = names.apply(lambda x:';'.join(e for e in x.split(';') if len(e) <8))
    df_group2.loc[:,'annotation'] = newnames
    df_group2.to_csv(folder+'20170916AgProteomeAll_group.csv')
    