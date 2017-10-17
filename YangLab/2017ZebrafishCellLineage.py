# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:46:06 2017

@author: ATPs
"""

def downloadConvertSRAFiles20170921():
    '''
    download all sra files from the paper
    Whole organism lineage tracing by combinatorial and cumulative genome editing
    '''
    #readin srrs
    import pandas as pd
    df = pd.read_excel(r"C:\Users\ATPs\OneDrive\Lab\YangLab\2017CellLineage\zebraFishPaper\2017zebrafish_SraRunInfo.xlsx")
    _srrs = list(df.loc[:,'Run'])
    
    #download files
    fout = open('20170921downloadZebrafishSRA.txt','w')
    fout.write('cd /mnt/data/home/ATPs/W/zebrafish/SRA/\n')
    for n in range(len(_srrs)):
        _srr = _srrs[n]
        fout.write('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'+_srr[:3]+'/'+_srr[:6]+'/'+_srr+'/'+_srr+'.sra')
        if (n+1)%6 !=0:
            fout.write(' &\n')
        else:
            fout.write(' \n')
            fout.write('wait $(jobs -p)\n')
    fout.close()
    
    #convert to fastq format
    fout = open('20170921convertZebrafishSRA.txt','w')
    for _srr in _srrs:
        fout.write('''/mnt/data/Programs/sratoolkit.2.8.2/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri'  --split-files /mnt/data/home/ATPs/W/zebrafish/SRA/%s.sra -O /mnt/data/home/ATPs/W/zebrafish/Reads/ --gzip &\n'''%_srr)
    fout.close()
    
    