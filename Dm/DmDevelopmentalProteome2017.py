# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 15:37:04 2017

@author: k
"""

def dmProteomeDevelopment20170417():
    '''
    deal with the result of maxquant search
    '''
    import pandas as pd
    #read in proteinGroups
    f_proteinGroups1 = r"D:\Insects\Drosophila_melanogaster\2017DmProteomeDevelopmental\txt\proteinGroups.txt"
    df_pr1 = pd.read_csv(f_proteinGroups1,sep = '\t')
    f_proteinGroups2 = r"D:\Insects\Drosophila_melanogaster\2017DmProteomeDevelopmental\txt2\proteinGroups.txt"
    df_pr2 = pd.read_csv(f_proteinGroups2,sep = '\t')
    
    #remove last 9 columns 
    df_pr1 = df_pr1.iloc[:,:-9]
    df_pr2 = df_pr2.iloc[:,:-9]
    
    #remove contaminant proteins
    df_pr1 = df_pr1[df_pr1.loc[:,"Potential contaminant"] != '+']
    df_pr2 = df_pr2[df_pr2.loc[:,"Potential contaminant"] != '+']
    
    #keep only those normal Dm Proteins
    df_pr1 = df_pr1[df_pr1.iloc[:,0].apply(lambda x:x[:4] == 'FBtr')]
    df_pr2 = df_pr2[df_pr2.iloc[:,0].apply(lambda x:x[:4] == 'FBtr')]
    
    #save csv file
    df_pr1.to_csv(r"D:\Insects\Drosophila_melanogaster\2017DmProteomeDevelopmental\proteinGroups20170417Em.csv")
    df_pr2.to_csv(r"D:\Insects\Drosophila_melanogaster\2017DmProteomeDevelopmental\proteinGroups20170417De.csv")
    
    #read the excel file of developmental stages, add column of protein ID and gene id
    df_dev = pd.read_excel(r"D:\Insects\Drosophila_melanogaster\2017DmProteomeDevelopmental\proteinGroups20170417DeLFQ.xlsx")
    ##read in Dm gene, transcript, protein file
    dc_t2p = {}
    dc_t2n = {}
    for line in open(r"D:\Insects\Drosophila_melanogaster\fbgn_fbtr_fbpp_fb_2017_01.tsv"):
        if len(line) >10 and line[:4] == 'FBgn':
            eles = line.split()
            if len(eles) == 2:
                dc_t2n[eles[1]] = eles[0]
                dc_t2p[eles[1]] = 'NA'
            else:
                dc_t2n[eles[1]] = eles[0]
                dc_t2p[eles[1]] = eles[2]
    
    for index in df_dev.index:
        fts = df_dev.loc[index,'Protein IDs'].split(';')
        fps = ';'.join([dc_t2p[ft] for ft in fts if ft in dc_t2p])
        fns = ';'.join(set([dc_t2n[ft] for ft in fts if ft in dc_t2n]))
        df_dev.loc[index,'proteinID'] = fps
        df_dev.loc[index,'geneID'] = fns
#        print(fps,fns)
#        break
    df_dev.to_csv(r"D:\Insects\Drosophila_melanogaster\2017DmProteomeDevelopmental\proteinGroups20170417DeLFQ.csv")
    
    #get LFQs for SPSPHs
    ls = set(open('list.txt').read().split())
    dcGene2index = {}
    for index in df_dev.index:
        geneIDs = df_dev.loc[index,'geneID'].split(';')
        for geneID in geneIDs:
            if geneID not in dcGene2index:
                dcGene2index[geneID] = []
            dcGene2index[geneID].append(index)
    dcSP2index = {}
    SPnoExpression = []
    for geneID in ls:
        if geneID not in dcGene2index:
            SPnoExpression.append(geneID)
        else:
            dcSP2index[geneID] = dcGene2index[geneID][0]
    df_SP = df_dev.loc[dcSP2index.values(),:]
    df_SP.to_csv('DmExpression.csv')

