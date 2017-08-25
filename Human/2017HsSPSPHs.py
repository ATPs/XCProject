# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:22:57 2017

@author: k
"""

def getHumanSPSPHs20170622():
    '''
    with domain predicted by SMART and PFAM, annotate human SP genes
    '''
    #filter SP domain containing genes
    #SM00020, Serine proteases, trypsin domain
    #PF00089, Serine proteases, trypsin domain
    lsDomainLines = open(r"D:\Insects\human\2017RefSeq\interim_GRCh38.p10.RefSeqProts.pfam.tab").readlines()+open(r"D:\Insects\human\2017RefSeq\interim_GRCh38.p10.RefSeqProtsUniNonwithin.fastaSmartTab2").readlines()
    lsDomainLines = [e for e in lsDomainLines if e[0]!='#']
    lsDomainLinesSP = [e for e in lsDomainLines if 'SM00020' in e.split()[3] or 'PF00089' in e.split()[3]]
    SPids = set([e.split()[0] for e in lsDomainLinesSP])
    
    #get sequences in SPids
    import largeFastaFile
    dcHs = largeFastaFile.open_fasta_to_dic(r"D:\Insects\human\2017RefSeq\interim_GRCh38.p10.RefSeqProtsUniNonwithin.fasta")
    lsSPs = [dcHs[e] for e in dcHs if e in SPids]
    largeFastaFile.saveFastaListToFile(lsSPs,r"D:\Insects\human\2017HsSPSPH\20170622HumanSPs.txt")
    
    #get human accession number gene relationship
    fout = open(r"D:\Insects\human\2017RefSeq\gene2accession_human",'w')
    fo = open(r"D:\Insects\human\2017RefSeq\gene2accession")
    line = fo.readline()
    fout.write(line)
    for line in fo:
        if line.split()[0] == '9606':
            fout.write(line)
    fout.close()
    fo.close()
    
    #get SPids and their gene id
    import pandas as pd
    df = pd.read_csv(r"D:\Insects\human\2017RefSeq\gene2accession_human",sep='\t')
    df_SP = pd.DataFrame()
    df_SP.loc[:,'accession'] = list(SPids)
    df_SP = df_SP.merge(df,how='left',left_on = 'accession',right_on = 'protein_accession.version')
    df_SP = df_SP.iloc[:,[0,2,4]]
    df_SP.to_csv('list.csv')
    
    #the gene id in the gene2accession is not so good. only 266 sequences belongs to 123 genes, seems unreasonable
    #use the gff file of refseq to get accession-gene relationship
    fo = open(r"D:\Insects\human\2017RefSeq\GCF_000001405.36_GRCh38.p10_genomic.gff")
    dcSPgene = {}
    for line in fo:
        if line[0]!= '#':
            if 'protein_id' in line:
#                break
                des = line.split('\t')[-1]
                proteinid = des.split('=')[-1].replace('\n','')
                geneid = des.split('GeneID:')[1].split(',')[0]
                if proteinid in SPids:
                    dcSPgene[proteinid] = geneid
    print(len(set(dcSPgene.values())))
    #still 123 genes. Which means, there is only 123 SPSPHs in human genome.
    
    