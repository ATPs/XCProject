# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 11:09:18 2016

@author: k

do Trinity assembling for insect with RNA-seq data
"""

def lepidoptera_TSA20161215():
    """
    deal with lepidoptera species.
    with less than 10 runs
    """
    
    folder = 'D:\\mine\\OneDrive\\Lab\\works\\2016GeneModelingInsectsTSA\\Lepidoptera_TSA\\'
    f_sraInfo = '20161215SRA_Lepidoptera_RNAseq.txt'
#    import numpy as np
    import pandas as pd
    df = pd.read_csv(folder + f_sraInfo, sep='\t') #store SRX information
    
    dc_speciesSRX = {} # dictionary use species as key, list of SRX as value
    for n in range(len(df)):
        _s = df.iloc[n,2]
        _x = df.iloc[n,0]
        if _s not in dc_speciesSRX:
            dc_speciesSRX[_s] = []
        dc_speciesSRX[_s].append(_x)
    
    f_srrInfo = 'SraRunInfo.csv'
    dfR = pd.read_csv(folder+f_srrInfo) #store SRR information
    dc_SRX_SRR = {} # dictionary use SRX as key, list of SRR as value
    ls_pairedSRR = [] #store paired SRRs
    for n in range(len(dfR)):
        _r = dfR.iloc[n,0]#SRR
        _x = dfR.iloc[n,10]#SRX
        _ps = dfR.iloc[n,15]#paired or single
        if _ps == 'PAIRED':
            ls_pairedSRR.append(_r)
        if _x not in dc_SRX_SRR:
            dc_SRX_SRR[_x] = []
        dc_SRX_SRR[_x].append(_r)
    
    #save file, store species and SRR
    #species name: SRRs
    dc_speciesSRR = {} #species as key, list of SRR as value
    f_out = open(folder+'species2SRR.txt','w')
    for _s in dc_speciesSRX:
        _x = dc_speciesSRX[_s]
        _rs = []
        for _xx in _x:
            _r = dc_SRX_SRR[_xx]
            _rs = _rs + _r
        if len(_rs) == 0:
            print(_s)
        else:
            f_out.write(_s + ':'+'\t'.join(_rs)+'\n')
            dc_speciesSRR[_s] = _rs
    f_out.close()
    
    #writing scripts for qsub
    upperlimitSRR = 10
    for _s in dc_speciesSRR:
        _sp =  _s.replace(' ','_')
        _filename = '20161215_qsub_' + _sp +'.txt'
        _srrs = dc_speciesSRR[_s]
        if len(dc_speciesSRR[_s]) <= upperlimitSRR:
            fout = open(folder+'scriptsForLessThan10SRR\\'+_filename,'w')
            fout.write('#!/bin/bash\n')
            fout.write('#PBS -q batch\n')
            fout.write('#PBS -N %s\n'%_sp)
            fout.write('#PBS -l nodes=1:ppn=12\n')
            fout.write('#PBS -l walltime=120:00:00\n')
            fout.write('#PBS -m abe -M atps@outlook.com\n')
            fout.write('#PBS -j oe\n')
            fout.write('module load python\n')
            fout.write('module load perl\n')
            fout.write('module load gcc-4.9.2\n')
            fout.write('module load jdk/1.8.0_45\n')            
            fout.write('module load bowtie2\n')
            fout.write('\n')
            fout.write('cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA\n')
            fout.write('mkdir %s\n'%_sp)
            fout.write('cd %s\n'%_sp)
            fout.write('mkdir RNA\n')
            fout.write('cd RNA')
            
            fout.write('\n')
            for _srr in _srrs:
                fout.write('''/panfs/panfs.cluster/home/ks2073/p/2016/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri'  --split-files %s\n'''%_srr)
            
            fout.write('cat *.fastq >> ../all.fastq\n')
            fout.write('rm *.fastq\n') #combine fastq files
            
            fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --trimmomatic ')
            fout.write(' --quality_trimming_params "ILLUMINACLIP:/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50" ')
            fout.write(' --seqType fq --max_memory 30G --single ../all.fastq ')
            fout.write(' --CPU 12 --no_normalize_reads --min_kmer_cov 2   --no_distributed_trinity_exec  ')
            fout.write(' --output ../Trinity\n')
            fout.close()
    
    #generate qsub scripts for species with more than 10 libraries
    #download only
    upperlimitSRR = 10
    for _s in dc_speciesSRR:
        _sp =  _s.replace(' ','_')
        _filename = '20161215_qsub_' + _sp +'_download.txt'
        _srrs = dc_speciesSRR[_s]
        if len(dc_speciesSRR[_s]) > upperlimitSRR:
            fout = open(folder+'scriptsForGreaterThan10SRR\\'+_filename,'w')
            fout.write('#!/bin/bash\n')
            fout.write('#PBS -q batch\n')
            fout.write('#PBS -N %s\n'%_sp)
            fout.write('#PBS -l nodes=1:ppn=12\n')
            fout.write('#PBS -l walltime=96:00:00\n')
            fout.write('#PBS -m abe -M atps@outlook.com\n')
            fout.write('#PBS -j oe\n')
            fout.write('module load python\n')
            fout.write('module load perl\n')
            fout.write('module load gcc-4.9.2\n')
            fout.write('module load jdk/1.8.0_45\n')            
            fout.write('module load bowtie2\n')
            fout.write('\n')
            fout.write('cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA2\n')
            fout.write('mkdir %s\n'%_sp)
            fout.write('cd %s\n'%_sp)
            fout.write('mkdir RNA\n')
            fout.write('cd RNA')
            
            fout.write('\n')
            for n in range(len(_srrs)):
                _srr = _srrs[n]
                fout.write('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'+_srr[:3]+'/'+_srr[:6]+'/'+_srr+'/'+_srr+'.sra')
                if (n+1)%4 !=0:
                    fout.write(' &\n')
                else:
                    fout.write(' \n')
            fout.write('wait $(jobs -p)\n')
            for _srr in _srrs:
                fout.write('''/panfs/panfs.cluster/home/ks2073/p/2016/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri'  --split-files ./%s.sra  '''%_srr)
                fout.write(''' && rm ./%s.sra  &\n'''%_srr)
            fout.write('\n')
            fout.write('wait $(jobs -p)\n')
            for n in range(len(_srrs)):
                _srr = _srrs[n]
                for ele in [_srr+'_1.fastq',_srr+'_2.fastq']:
                    fout.write('java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic/trimmomatic.jar ')
                    fout.write(' SE -threads 2 -phred33 ')
                    fout.write(ele +' ' + ele+'.trim ')
                    fout.write('ILLUMINACLIP:/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:25 && ')
                    fout.write('/panfs/panfs.cluster/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i ' + ele +'.trim -o ' + ele + '.fa -n -Q33 &\n' )
                if (n+1)%3 == 0:
                    fout.write('wait $(jobs -p)\n')
            fout.write('wait $(jobs -p)\n')
            fout.write('cat *.fa >>../single.fa && ')
            fout.write('rm *\n')
            fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity ')
            fout.write(' --seqType fa --max_memory 30G --single ../single.fa ')
            fout.write(' --CPU 12 --no_normalize_reads --min_kmer_cov 1   --no_distributed_trinity_exec  ')
            fout.write(' --output ../Trinity\n')
            fout.close()
    
    #trim reads, and combine to fasta format 20170102
    #list.txt stores filenames of fastq files
    fo = open('list.txt')
    ls = fo.read().split()
    fo.close()
    for ele in ls:
        _e = ele.split('/')[-1]
        fout = open('20170202Trim'+_e+'.txt','w')
        fout.write('#!/bin/bash\n')
        fout.write('#PBS -q batch\n')
        fout.write('#PBS -N %s\n'%_e)
        fout.write('#PBS -l nodes=1:ppn=12\n')
        fout.write('#PBS -l walltime=12:00:00\n')
        fout.write('#PBS -m abe -M atps@outlook.com\n')
        fout.write('#PBS -j oe\n')
        fout.write('module load python\n')
        fout.write('module load perl\n')
        fout.write('module load gcc-4.9.2\n')
        fout.write('module load jdk/1.8.0_45\n')
        fout.write('java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic/trimmomatic.jar ')
        fout.write(' SE -threads 12 -phred33 ')
        fout.write(ele +' ' + ele+'.trim ')
        fout.write('ILLUMINACLIP:/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50\n')
        fout.write('/panfs/panfs.cluster/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i ' + ele +'.trim -o ' + ele + '.fa -n -Q33 \n' )
        fout.close()

    

def speedUpTrinitySplitRun20161216():
    """
    
    """
    import argparse
    parser = argparse.ArgumentParser(description = 'splitFile')
    parser.add_argument('f', type = str, help = 'folder')
    arg = parser.parse_args()
    f = arg.f
    
    import os
    os.chdir(f)
    files = os.listdir('./')
    
    s1 = set(open('recursive_trinity.cmds').readlines())
    if 'recursive_trinity.cmds.completed' in files:
        s2 = set(open('recursive_trinity.cmds.completed').readlines())
    else:
        s2 = set()
    
    s3 = s1-s2
    
    outfolder = r'/panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA/0cmds/cmdjobs'
    
    dc = {}
    for i in range(40):
        dc[i] = open(outfolder+str(i),'w')
    
    s3 = list(s3)
    for n in range(len(s3)):
        i = n % 40
        dc[i].write(s3[n])
    
    for i in range(40):
        dc[i].close()            
    
    
    #generate qsub files
    for i in range(40):
        fout = open('20161216TrinityCursive.txt%d'%i,'w')
        fout.write('#!/bin/bash\n')
        fout.write('#PBS -q batch\n')
        fout.write('#PBS -N cur%d\n'%i)
        fout.write('#PBS -l nodes=1:ppn=12\n')
        fout.write('#PBS -l walltime=120:00:00\n')
        fout.write('#PBS -m abe -M atps@outlook.com\n')
        fout.write('#PBS -j oe\n')
        fout.write('module load python\n')
        fout.write('module load perl\n')
        fout.write('module load gcc-4.9.2\n')
        fout.write('module load jdk/1.8.0_45\n')            
        fout.write('module load bowtie2\n')
        fout.write('\n')
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/parafly-code/bin/ParaFly -c /panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA/0cmds/cmdjobs%d -CPU 12 -v \n'%i)
        fout.close()
        
def dealwithTrinityAssemblyResult20170102():
    '''
    gene names of Trinity looks like 'TRINITY_DN87_c1_g1_i1 len=243 path=[221:0-242] [-1, 221, -2]'
    simplify to T87_c1_g1_i1
    '''
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA/result/'
    import os
    files = os.listdir(folder)
    from Bio import SeqIO
    outfolder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA/DNA/'
    for file in files:
        fout = open(outfolder + file,'w')
        for ele in SeqIO.parse(open(folder + file),'fasta'):
            fout.write('>'+ele.id.replace('TRINITY_DN','T') +'\n' + str(ele.seq) +'\n')
        fout.close()

def translateTSAResult20170102():
    '''
    translate transcripts to proteins
    '''
    files = open('list.txt').read().split()
    for file in files:
        fn = file.split('/')[-1]
        fout = open('20170102translation'+fn,'w')
        fout.write('#!/bin/bash\n')
        fout.write('#PBS -q batch\n')
        fout.write('#PBS -N %s\n'%fn)
        fout.write('#PBS -l nodes=1:ppn=6\n')
        fout.write('#PBS -l walltime=12:00:00\n')
        fout.write('#PBS -m abe -M atps@outlook.com\n')
        fout.write('#PBS -j oe\n')
        fout.write('module load cd-hit\n')
        fout.write('module load perl\n')
        fout.write('module load gcc-4.9.2\n')
        fout.write('module load jdk/1.8.0_45\n')            
        fout.write('\n')
        fout.write('cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA/protein/\n')
        fout.write('mkdir ' + fn +'\n')
        fout.write('cd ' +fn +'\n')
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/TransDecoder-3.0.0/TransDecoder.LongOrfs -m 60 -t '+file +' \n')
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/TransDecoder-3.0.0/TransDecoder.Predict -t '+file + ' --cpu 12 \n')
        fout.write('')
        fout.close()

def simplifyTransdecoderResultName20170102():
    '''
    change name from
    >T100074_c0_g2::T100074_c0_g2_i2::g.42404::m.42404 T100074_c0_g2::T100074_c0_g2_i2::g.42404  ORF type:complete len:67 (-) T100074_c0_g2_i2:97-297(-)
    to
    >m.42404 T100074_c0_g2_i2 97-297(-) complete
    '''
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/TSA/result/'
    import os
    files = os.listdir(folder)
    from Bio import SeqIO
    for file in files:
        file = folder + file
        fn = file.replace('result','protein')
        fn = fn.replace('.fasta.','.')
        fout = open(fn,'w')
        for ele in SeqIO.parse(open(file),'fasta'):
            _ls = ele.description.split(' ')
            part1 = _ls[0].split(':')[-1]
            part2 = _ls[0].split(':')[2]
            part3 = _ls[-1].split(':')[-1]
            part4 = _ls[4].split(':')[-1]
            name = part1 + ' ' + part2 + ' ' + part3 + ' ' + part4
            fout.write('>' + name +'\n' + str(ele.seq)+'\n')
        fout.close()



