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
    
def annotateMCOT2forProteome20171018():
    '''
    find MCOT2 sequences with no match with AGAP proteins
    '''
    filename = r"X:\Insects\Anopheles_gambiae\MCOT2\MCDpeptides.fasta"
    fout = open(r"X:\Insects\Anopheles_gambiae\MCOT2\MCDpeptidesToAnnotate20171018.fasta",'w')
    import re
    from Bio import SeqIO
    for s in SeqIO.parse(open(filename),'fasta'):
        if not re.findall('AGAP\d\d\d\d\d\d', s.description):
            fout.write('>'+s.description+'\n'+str(s.seq)+'\n')
    fout.close()

def addNamesForResult20171028():
    '''
    simplify the naming column for the two xlsx files
    '''
    import pandas as pd
    df = pd.read_excel(r"C:\Users\ATPs\OneDrive\Lab\works\2017DmSPSPH\20170823AgMassSpecDataAnalysisAll\20170916AgProteomeAll_group.xlsx")
    def preferenceFunc(x):
        if "AGAP" in x:
            return 2
        if "MCD" in x:
            return 3
        return 1
    for _row in df.index:
        _r = df.loc[_row,:]
        _proteins = _r['Protein IDs'].split(';')
        _counts = [int(e) for e in _r['Peptide counts (all)'].split(';')]
        _protein_keep = _proteins[:_counts.count(_counts[0])]
        _protein_keep.sort(key = lambda x: preferenceFunc(x))
        _proteinID = _protein_keep[0]
        df.loc[_row,'Protein Best'] = _proteinID
    
    df.to_csv('20170916AgProteomeAll_group.csv')
    
    #get AGAP naming
    from Bio import SeqIO
    ls_fa = list(SeqIO.parse(open(r"X:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.7.fa"),'fasta'))
    ls_name = []
    for _s in ls_fa:
        _id = _s.id
        _name = _s.description.split(' ',1)[1].split('|')[0]
        ls_name.append((_id,_name))
    fout = open('list.txt','w')
    for _id, _name in ls_name:
        fout.write(_id+'\t'+_name+'\n')
    fout.close()
        
            