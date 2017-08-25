# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 21:44:36 2017

@author: k
"""

def ApTSA20170612():
    import pandas as pd
    f_srrInfo = r"C:\Users\k\OneDrive\Lab\works\2016GeneModelingInsectsTSA\ApTSA\SraRunInfo.csv"
    dfR = pd.read_csv(f_srrInfo)
    SRRs = list(dfR.iloc[:,0])
    for _srr in SRRs:
        fout = open('2017ApSRR_%s.txt'%(_srr),'w')
        fout.write('#!/bin/bash\n')
        fout.write('#PBS -q batch\n')
        fout.write('#PBS -N%s\n'%(_srr))
        fout.write('#PBS -l nodes=1:ppn=6\n')
        fout.write('#PBS -l walltime=120:00:00\n')
        fout.write('#PBS -m abe -M atps@outlook.com\n')
        fout.write('#PBS -j oe\n')
        fout.write('module load python\n')
        fout.write('module load perl\n')
        fout.write('module load gcc-4.9.2\n')
        fout.write('module load jdk/1.8.0_45\n')            
        fout.write('module load bowtie2\n')
        fout.write('\n')
        fout.write('cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/SRR\n')
        fout.write('mkdir '+_srr+'\ncd ' + _srr +'\n')
        fout.write('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'+_srr[:3]+'/'+_srr[:6]+'/'+_srr+'/'+_srr+'.sra &&\n')
        fout.write('''/panfs/panfs.cluster/home/ks2073/p/2016/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri'  --split-files ./%s.sra  '''%_srr)
        fout.write(''' && rm ./%s.sra  \n'''%_srr)
        for ele in [_srr+'_1.fastq',_srr+'_2.fastq']:
            fout.write('java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic/trimmomatic.jar ')
            fout.write(' SE -threads 6 -phred33 ')
            fout.write(ele +'.gz ' + ele+'.trim ')
            fout.write('ILLUMINACLIP:/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50 \n ')
            fout.write('/panfs/panfs.cluster/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i ' + ele +'.trim -o ' + ele + '.fa -r -Q33 \n' )
            fout.write('rm ' + ele +'.gz\n')
            fout.write('rm ' + ele +'.trim\n')
        fout.write('cat ' + _srr+'_2.fastq.fa >> ' + _srr+'_1.fastq.fa && rm ' + _srr+'_2.fastq.fa\n')
#        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/util/insilico_read_normalization.pl --seqType fa --JM 30G --max_cov 250 --output ./Trinity --CPU 12 --single '  + _srr+'_1.fastq.fa\n')
#        fout.write('mv ./Trinity/'+ _srr+'_1.fastq.fa.normalized_K25_C250_pctSD200.fa ../../RNA/'+_srr+'.fa\n')
#        fout.write('rm -rf ../'+_srr+'\n')
        fout.close()

def TcTSA20170612():
    import pandas as pd
    f_srrInfo = r"C:\Users\k\OneDrive\Lab\works\2016GeneModelingInsectsTSA\TcTSA\SraRunInfo.csv"
    dfR = pd.read_csv(f_srrInfo)
    SRRs = list(dfR.iloc[:,0])
    for _srr in SRRs:
        fout = open('2017TcSRR_%s.txt'%(_srr),'w')
        fout.write('#!/bin/bash\n')
        fout.write('#PBS -q batch\n')
        fout.write('#PBS -N%s\n'%(_srr))
        fout.write('#PBS -l nodes=1:ppn=6\n')
        fout.write('#PBS -l walltime=120:00:00\n')
        fout.write('#PBS -m abe -M atps@outlook.com\n')
        fout.write('#PBS -j oe\n')
        fout.write('module load python\n')
        fout.write('module load perl\n')
        fout.write('module load gcc-4.9.2\n')
        fout.write('module load jdk/1.8.0_45\n')            
        fout.write('module load bowtie2\n')
        fout.write('\n')
        fout.write('cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Tc/SRR\n')
        fout.write('mkdir '+_srr+'\ncd ' + _srr +'\n')
        fout.write('wget -q ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/'+_srr[:3]+'/'+_srr[:6]+'/'+_srr+'/'+_srr+'.sra &&\n')
        fout.write('''/panfs/panfs.cluster/home/ks2073/p/2016/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri'  --split-files ./%s.sra  '''%_srr)
        fout.write(''' && rm ./%s.sra  \n'''%_srr)
        for ele in [_srr+'_1.fastq',_srr+'_2.fastq']:
            fout.write('java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic/trimmomatic.jar ')
            fout.write(' SE -threads 6 -phred33 ')
            fout.write(ele +'.gz ' + ele+'.trim ')
            fout.write('ILLUMINACLIP:/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50 \n ')
            fout.write('/panfs/panfs.cluster/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i ' + ele +'.trim -o ' + ele + '.fa -r -Q33 \n' )
            fout.write('rm ' + ele +'.gz\n')
            fout.write('rm ' + ele +'.trim\n')
        fout.write('cat ' + _srr+'_2.fastq.fa >> ' + _srr+'_1.fastq.fa && rm ' + _srr+'_2.fastq.fa\n')
#        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/util/insilico_read_normalization.pl --seqType fa --JM 30G --max_cov 250 --output ./Trinity --CPU 12 --single '  + _srr+'_1.fastq.fa\n')
#        fout.write('mv ./Trinity/'+ _srr+'_1.fastq.fa.normalized_K25_C250_pctSD200.fa ../../RNA/'+_srr+'.fa\n')
#        fout.write('rm -rf ../'+_srr+'\n')
        fout.close()

def reads2transcriptsAp20170613():
    '''
    generate running scripts
    '''
    filenames = open('list.txt').read().split()
    for filename in filenames:
        fout = open('2017DmReads2Transcripts_'+filename,'w')
        fout.write('''
    #!/bin/bash
    #PBS -q batch
    #PBS -N %s
    #PBS -l nodes=1:ppn=12
    #PBS -l walltime=120:00:00
    #PBS -m abe -M atps@outlook.com
    #PBS -j oe
    
    module reset
    module load python
    module load perl
    module load gcc-4.9.2
    module load jdk/1.8.0_45
    module load bowtie2
    
    cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/RNA
    /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Chrysalis/ReadsToTranscripts -i %s -f /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/Trinity/bundled_iworm_contigs.fasta -o ../reads/%s  -t 12 -max_mem_reads 50000000 
    /bin/sort -T . -S 30G -k 1,1n ../reads/%s > ../reads/%s.sort
    rm ../reads/%s
        '''%(filename,filename,filename,filename,filename,filename))
        fout.close()

def reads2transcriptsTc20170613():
    '''
    generate running scripts
    '''
    filenames = open('list.txt').read().split()
    for filename in filenames:
        fout = open('2017TcReads2Transcripts_'+filename,'w')
        fout.write('''
    #!/bin/bash
    #PBS -q batch
    #PBS -N %s
    #PBS -l nodes=1:ppn=12
    #PBS -l walltime=120:00:00
    #PBS -m abe -M atps@outlook.com
    #PBS -j oe
    
    module reset
    module load python
    module load perl
    module load gcc-4.9.2
    module load jdk/1.8.0_45
    module load bowtie2
    
    cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Tc/RNA
    /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Chrysalis/ReadsToTranscripts -i %s -f /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Tc/Trinity/bundled_iworm_contigs.fasta -o ../reads/%s  -t 12 -max_mem_reads 50000000 
    /bin/sort -T . -S 30G -k 1,1n ../reads/%s > ../reads/%s.sort
    rm ../reads/%s
        '''%(filename,filename,filename,filename,filename,filename))
        fout.close()

def generateQsubForTrinityFinal20170614():
    #Ap
    seqs = range(0,158780)
    fout = open('recursive_trinity.cmds','w')
    for seqn in seqs:
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/trinityout{0} --CPU 1 --max_memory 2G --seqType fa --trinity_complete --full_cleanup  \n'.format(seqn))
    fout.close()
    
    #Tc
    seqs = range(0,154807)
    fout = open('recursive_trinity.cmds','a')
    for seqn in seqs:
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Tc/seqs/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Tc/TrinityOne/trinityout{0} --CPU 1 --max_memory 2G --seqType fa --trinity_complete --full_cleanup  \n'.format(seqn))
    fout.close()

def collectFastaAp20170222():
    import os
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/'
    fout = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/ApTrinity.fa','w')
    from Bio import SeqIO
    files = os.listdir(folder)
    for file in files:
        if file[-6:] == '.fasta':
            head = file.split('trinityout')[1].split('.')[0]
            ls = list(SeqIO.parse(open(folder+file),'fasta'))
            for seq in ls:
                fout.write('>Ap'+head+seq.id+'\n'+str(seq.seq)+'\n')
    fout.close()
def collectFastaTc20170222():
    import os
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Tc/TrinityOne/'
    fout = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TcTrinity.fa','w')
    from Bio import SeqIO
    files = os.listdir(folder)
    for file in files:
        if file[-6:] == '.fasta':
            head = file.split('trinityout')[1].split('.')[0]
            ls = list(SeqIO.parse(open(folder+file),'fasta'))
            for seq in ls:
                fout.write('>Tc'+head+seq.id+'\n'+str(seq.seq)+'\n')
    fout.close()

