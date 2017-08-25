# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 14:05:58 2016

@author: k
"""

def getIndex20161130():
    """
    get indexes used for 72 files.
    reads
    @J00138:89:HG7NLBBXX:6:1101:2676:1068 1:N:0:NGACCA
    GNCGATCACACAGACACGTCCAGCGCTGCAGTGAACCGTCAGACACGGATCTTCCATGTTTAGTTCATTGTCAGTA
    +
    -#AAF<A--FJFFJFFJ<AFJJAFAFJ<JJ<--FJ<F-<JJ7FJJJFFJ-F<FAFJ<<-A-<-<7<FAA7AA<<--
    
    index is NGACCA
    
    run as python scripts directly
    """
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/'
    fout = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHindex.txt",'w')
    import os
    files = os.listdir(folder)
    for file in files:
        fo = open(folder + file)
        line = fo.readline()
        filename = file.split(".")[0]
        index = line.split(":")[-1]
        fout.write(filename + '\t'+index)
        fo.close()
    fout.close()
    
def fastqc20161130():
    """
    get fastQC information for 72 files
    generates 72 qsub scripts
    """
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/'
    files = open('list.txt').read().split()
    
    for filename in files:
        fo = open(filename,'w')
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+filename+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd " + folder +'\n')
        fo.write("/home/ks2073/p/FastQC0.11.3/fastqc -f fastq -o /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/QC ")
        fo.write(filename)
        fo.close()

def fastqTrim20161130():
    """
    get fastQC information for 72 files
    generates 72 qsub scripts
    """
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/'
    files = open('list.txt').read().split()
    
    for filename in files:
        fo = open(filename,'w')
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+filename+"\n")
        fo.write("#PBS -l nodes=1:ppn=4\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd " + folder +'\n')
        fo.write("module load jdk\n")
        fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic.jar SE -phred33 -threads 4 " +\
        filename+" ../CIFHG_Trim/" + filename +\
        " ILLUMINACLIP:/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/SE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50 ")

        fo.close()


def fastaDealwith20161201():
    """
    find some way to compress fastq seq
    test
    conclusion: will not try it any more. Take too long time, and only save a little spaces in disk
    """
    fo = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/4_R1_001.fastq.unique")
    fout = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/4_R1_001.fastq.unique_moreThan1",'w')
    from Bio import SeqIO
    for read in SeqIO.parse(fo,'fasta'):
        num = read.id.split('-')[-1][:-1]
        if num == '1':
            break
        else:
            fout.write(str(read.seq)+'\n')
    fout.close()
    fo.close()
    
    fo = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/4_R1_001.fastq.unique_moreThan1",'r')
    myset = set(fo.read().split())
    fo.close()
    fo = "/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/4_R1_001.fastq"
    fout = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/CIFHG/4_R1_001.fastq.order",'w')
    
    def fq_reader(fqfile):
        fo = open(fqfile)
        
        for i, line in enumerate(fo):
            j = i%4
            if j == 1:
                read_seq = line.rstrip('\n')
            if j == 3:
                read_qual = line.rstrip('\n')
                yield(read_seq, read_qual)
            
    fqs = fq_reader(fo)
    for read_seq in myset:
        tempout = open('/dev/shm/'+read_seq,'w')
        tempout.close()
    for read_seq, read_qual in fqs:
        if read_seq not in myset:
            fout.write("@1\n%s\n+\n%s\n"%(read_seq,read_qual))
        else:
            tempout = open('/dev/shm/'+read_seq,'r+')
            tempout.write(read_qual+'\n')
            tempout.close()
    files = myset
    for read_seq in files:
        for read_qual in open('/dev/shm/'+read_seq):
            fout.write("@1\n%s\n+\n%s\n"%(read_seq,read_qual))
    fout.close()
        
def rsemWithAnnotatedseq():
    """
    684 immune related genes. check their express specifically
    """
    ogsfile = r"E:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2_2016\ms_ogs_transcripts.fa"
    toremovefile = r"F:\Insects\ManducaSexta\AnnotatedSequences\RSEM_specific\idToRemove.txt"
    fout = open("20161201OGS2_new.txt",'w')
    toremove = set(open(toremovefile).read().split())
    from Bio import SeqIO
    for fa in SeqIO.parse(open(ogsfile),'fasta'):
        if fa.id not in toremove:
            fout.write(">"+fa.id+'\n'+str(fa.seq)+'\n')
    fout.close()
    
    toaddfile = r"F:\Insects\ManducaSexta\AnnotatedSequences\all.fa"
    fout = open("20161201OGS2_new2.txt",'w')
    nameOrder = open("20161201NameOrder.txt",'w')
    count = 0
    for fa in SeqIO.parse(open(toaddfile),'fasta'):
        count+=1
        fout.write(">zzz%d"%count+'\n'+str(fa.seq)+'\n')
        nameOrder.write(">zzz%d\t"%count + fa.id + '\n')
    fout.close()
    nameOrder.close()

def readsInforExtractionTrimmomatic20161202():
    folder = 'F:\\Insects\\ManducaSexta\\CIFGH\cont\\'
    import os
    files = os.listdir(folder)
    infos = []
    for file in files:
        fo = open(folder + file)
        lib = file.split('_')[0]
        librf = file.split('_')[1]
        templs = fo.readlines()
        templn = templs[7]
        readsn = templn.split()[2]
        readssur = templn.split()[4]
        infos.append([lib, librf, readsn, readssur])
        fo.close()
    fout = open('libReadsInfo.txt','w')
    for ele in infos:
        fout.write('\t'.join(ele) +'\n')
    fout.close()

def rsem_RunningInfor_Extraction20161205():
    """
    mapping rate of trimmed reads to cufflinks gene models, OGS2, MCOT
    """
    def rsemExtraction(folder, outname):
        """
        mapping rate of trimmed reads gene models
        """
        folder = folder
        import os
        files = os.listdir(folder)
        infos = []
        for file in files:
            fo = open(folder + file)
            templs = fo.readlines()
            fo.close()
            libn = file.split('.')[0][4:]
            unmap = templs[3].split()[0]
            map1 = templs[4].split()[0]
            mapmorethan1 = templs[5].split()[0]
            infos.append([libn, unmap, map1, mapmorethan1])
        fout = open(outname,'w')
        for ele in infos:
            fout.write('\t'.join(ele) +'\n')
        fout.close()
    
    folder = "F:\\Insects\\ManducaSexta\\CIFGH\RSEM\\RSEM_MCOT\\log\\"
    outname = "RSEMinfo_MCOT.txt"
    rsemExtraction(folder, outname)
    
    folder = "F:\\Insects\\ManducaSexta\\CIFGH\RSEM\\RSEMogs2\\log\\"
    outname = "RSEMinfo_ogs2.txt"
    rsemExtraction(folder, outname)
    

def rsem_result_extraction20161205():
    """
    combine running result for 36 libraries
    """
    import pandas as pd
    import numpy as np
    import os
    folder = 'F:\\Insects\\ManducaSexta\\CIFGH\RSEM\\RSEM_cufflinks\\'
    files = os.listdir(folder)
    files.sort(key = lambda x:int(x.split(".")[0]))
    
    #save FPKM
    fpkm = pd.read_csv(folder+files[0], sep = '\t').iloc[:,(0,2)]
    for file in files:
        df = pd.read_csv(folder+file, sep = '\t')
        lib = file.split('.')[0]
        fpkm[lib] = df.iloc[:,6]
    fpkm.to_csv('20161214CufflinksFPKM_CIFGH.csv')
    
    #save TPM
    fpkm = pd.read_csv(folder+files[0], sep = '\t').iloc[:,(0,2)]
    for file in files:
        df = pd.read_csv(folder+file, sep = '\t')
        lib = file.split('.')[0]
        fpkm[lib] = df.iloc[:,5]
    fpkm.to_csv('20161214Cufflinks_TPM_CIFGH.csv')
    
    #save expected_count
    fpkm = pd.read_csv(folder+files[0], sep = '\t').iloc[:,(0,2)]
    for file in files:
        df = pd.read_csv(folder+file, sep = '\t')
        lib = file.split('.')[0]
        fpkm[lib] = df.iloc[:,4]
    fpkm.to_csv('20161214Cufflinks_Count_CIFGH.csv')
    
    def saveCount_TPM_FPKM(folder, database = ''):
        import pandas as pd
        import numpy as np
        import glob
        files = glob.glob(folder+'\\*\\*.isoforms.results')
        files.sort(key = lambda x:int(x.split('\\')[-1].split(".")[0]))
        
        fpkm = pd.read_csv(files[0], sep = '\t').iloc[:,(0,2)]
        tpm = pd.read_csv(files[0], sep = '\t').iloc[:,(0,2)]
        rcount = pd.read_csv(files[0], sep = '\t').iloc[:,(0,2)]
        for file in files:
            df = pd.read_csv(file, sep = '\t')
            lib = file.split('\\')[-1].split(".")[0]
            fpkm[lib] = df.iloc[:,6]
            tpm[lib] = df.iloc[:,5]
            rcount[lib] = df.iloc[:,4]
        fpkm.to_csv('20161214%s_FPKM_CIFGH.csv'%database)
        tpm.to_csv('20161214%s_TPM_CIFGH.csv'%database)
        rcount.to_csv('20161214%s_count_CIFGH.csv'%database)
        
    
    folder = r'F:\Insects\ManducaSexta\CIFGH\RSEM\RSEM_MCOT'
    saveCount_TPM_FPKM(folder,'MCOT')
    
    folder = r'F:\Insects\ManducaSexta\CIFGH\RSEM\RSEMogs2'
    saveCount_TPM_FPKM(folder,'OGS2')
    
    folder = r'F:\Insects\ManducaSexta\CIFGH\RSEM\RSEM_annotated'
    saveCount_TPM_FPKM(folder,'Annotated')
    
    folder = r'F:\Insects\ManducaSexta\CIFGH\RSEM\RSEM_annotated'
    #keep target gene infor only. gene name begin with zzz
    import glob
    files = glob.glob(folder+'\\*')
    for file in files:
        fout = open(file +'.csv','w')
        fo = open(file)
        line = fo.readline()
        fout.write(line)
        for line in fo:
            if line.split(',')[1][:3] == 'zzz':
                fout.write(line)
        fo.close()
        fout.close()









    
        