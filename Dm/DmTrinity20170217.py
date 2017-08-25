# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 10:49:08 2017

@author: k
"""
def getDcSRR():
    '''
    return dcSRR shown below
    '''
    dcSRR = {}
    dcSRR['egg'] = '''SRR023199 SRR023659 SRR023663 SRR023671 SRR023747 SRR023755 SRR024013 SRR027110 SRR034220 SRR035417 SRR035418 SRR1197337 SRR1197370 SRR488290 SRR488297 SRR767606
    SRR023604 SRR023644 SRR023694 SRR023727 SRR023831 SRR024010 SRR024221 SRR034225 SRR035396 SRR035397 SRR1197334 SRR1197367 SRR767616
    SRR023602 SRR023643 SRR023703 SRR023724 SRR023832 SRR024012 SRR034226 SRR035221 SRR035398 SRR035399 SRR070438 SRR1197332 SRR1197369 SRR488302 SRR767613
    SRR023540 SRR023641 SRR023725 SRR023734 SRR023736 SRR023822 SRR024011 SRR034227 SRR035400 SRR035401 SRR1197331 SRR767618
    SRR023600 SRR023707 SRR023715 SRR023720 SRR023751 SRR023826 SRR034228 SRR035402 SRR1197330 SRR1197365 SRR488303 SRR767605
    SRR023503 SRR023599 SRR023673 SRR023682 SRR023686 SRR023699 SRR023700 SRR024215 SRR027114 SRR034229 SRR035222 SRR035403 SRR1197327 SRR1197363 SRR767622
    SRR023502 SRR023660 SRR023705 SRR023722 SRR023745 SRR027112 SRR034221 SRR1197336 SRR1197368 SRR488291 SRR488298 SRR767626
    SRR023506 SRR023665 SRR023684 SRR023718 SRR023728 SRR027109 SRR034230 SRR1197329 SRR1197364 SRR488304 SRR767620
    SRR023538 SRR023685 SRR023711 SRR023733 SRR023834 SRR024015 SRR034231 SRR035404 SRR1197328 SRR1197366 SRR767625
    SRR023539 SRR023669 SRR023696 SRR023746 SRR023836 SRR024014 SRR034222 SRR035220 SRR035405 SRR035406 SRR1197338 SRR488292 SRR767609
    SRR023504 SRR023654 SRR023668 SRR023688 SRR023691 SRR023732 SRR024009 SRR024217 SRR027113 SRR034223 SRR035407 SRR035408 SRR1197333 SRR488300 SRR767610
    SRR023549 SRR023596 SRR023603 SRR023657 SRR023701 SRR023749 SRR023750 SRR023754 SRR023759 SRR024219 SRR034224 SRR035409 SRR1197335 SRR488301 SRR767615
    '''.split()
    dcSRR['L1'] = '''SRR023597 SRR023646 SRR023661 SRR023666 SRR023706 SRR023835 SRR035410 SRR1197325 SRR1197426 SRR1509512 SRR1509513 SRR488305 SRR767624
    '''.split()
    dcSRR['L2'] = '''SRR023542 SRR023670 SRR023719 SRR023761 SRR023824 SRR024016 SRR035223 SRR035411 SRR035412 SRR1197324 SRR1197425 SRR488293 SRR767623
    '''.split()
    dcSRR['L3'] = '''SRR023507 SRR023649 SRR023677 SRR023731 SRR023760 SRR027111 SRR1197326 SRR1197424 SRR488306 SRR767627
    SRR023546 SRR023608 SRR023638 SRR023640 SRR023697 SRR023726 SRR023758 SRR035413 SRR1197312 SRR1197392 SRR1197412
    SRR023505 SRR023676 SRR023683 SRR023690 SRR023692 SRR023742 SRR027108 SRR1197308 SRR1197388 SRR1197467 SRR488294
    SRR023729 SRR023739 SRR023782 SRR023865 SRR026433 SRR035219 SRR035392 SRR1197307 SRR1197387 SRR1197408
    '''.split()
    dcSRR['pupa'] = '''SRR023667 SRR023721 SRR023743 SRR023785 SRR023829 SRR026431 SRR1197287 SRR1197420 SRR488296 SRR767614
    SRR023738 SRR023744 SRR023748 SRR023783 SRR024216 SRR026430 SRR035218 SRR035391 SRR1197285 SRR1197419 SRR767611
    SRR023609 SRR023653 SRR023723 SRR023735 SRR023740 SRR023827 SRR035416 SRR1197286 SRR1197416 SRR767612
    SRR023687 SRR023695 SRR023704 SRR023784 SRR023837 SRR026432 SRR035217 SRR1197290 SRR1197423 SRR1509514 SRR1509515 SRR767617
    SRR023544 SRR023639 SRR023647 SRR023689 SRR023716 SRR023833 SRR035414 SRR1197289 SRR1197422 SRR488295 SRR767608
    SRR023541 SRR023656 SRR023708 SRR023756 SRR023830 SRR035415 SRR1197288 SRR1197421 SRR767604
    '''.split()
    dcSRR['adult'] = '''SRR023595 SRR023662 SRR023680 SRR023702 SRR023737 SRR023828 SRR024220 SRR035393 SRR1197317 SRR1197418 SRR767621
    SRR023548 SRR023598 SRR023674 SRR023681 SRR023712 SRR023757 SRR029232 SRR029236 SRR1197314 SRR1197394 SRR1197411
    SRR023547 SRR023607 SRR023645 SRR023651 SRR023717 SRR023730 SRR029230 SRR029234 SRR1197313 SRR1197393 SRR1197414 SRR1197427 SRR1197428 SRR1197433 SRR1197434 SRR1197464 SRR1197465 SRR1197470 SRR1197471
    SRR023543 SRR023648 SRR023664 SRR023678 SRR023741 SRR023825 SRR035394 SRR1197315 SRR1197415 SRR767619
    SRR023601 SRR023650 SRR023698 SRR023714 SRR023752 SRR023823 SRR024218 SRR035395 SRR1197311 SRR1197391 SRR1197413
    SRR023550 SRR023605 SRR023606 SRR023642 SRR023658 SRR023672 SRR023679 SRR023713 SRR029176 SRR029231 SRR029233 SRR029235 SRR1197316 SRR1197417 SRR1197429 SRR1197430 SRR1197431 SRR1197432 SRR1197435 SRR1197436 SRR1197468 SRR1197469 SRR767607
    '''.split()
    return dcSRR

def downloadSRRconvertNormalize20170212():
    '''
    generate qsub script to download SRR files
    '''
    dcSRR = getDcSRR()
    for key in dcSRR:
        for _srr in dcSRR[key]:
            fout = open('2017DmSRR_%s_%s.txt'%(key,_srr),'w')
            fout.write('#!/bin/bash\n')
            fout.write('#PBS -q batch\n')
            fout.write('#PBS -N %s%s\n'%(key,_srr))
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
            fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/util/insilico_read_normalization.pl --seqType fa --JM 30G --max_cov 250 --output ./Trinity --CPU 12 --single '  + _srr+'_1.fastq.fa\n')
            fout.write('mv ./Trinity/'+ _srr+'_1.fastq.fa.normalized_K25_C250_pctSD200.fa ../../RNA/'+_srr+'.fa\n')
            fout.write('rm -rf ../'+_srr+'\n')
            fout.close()
    
    srrs = []
    for key in dcSRR:
        srrs += dcSRR[key]
    print(len(srrs))
    
    #download from EBI
    for key in dcSRR:
        for _srr in dcSRR[key]:
            fout = open('2017DmSRR_%s_%s.txt'%(key,_srr),'w')
            fout.write('#!/bin/bash\n')
            fout.write('#PBS -q batch\n')
            fout.write('#PBS -N %s%s\n'%(key,_srr))
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
            fout.write('cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/SRR\n')
            fout.write('mkdir '+_srr+'\ncd ' + _srr +'\n')
            for ele in [_srr+'_1.fastq',_srr+'_2.fastq',_srr+'.fastq']:
                fout.write('wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+_srr[:6]+'/'+_srr+'/'+ele+'.gz &&\n')
                fout.write('java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic/trimmomatic.jar ')
                fout.write(' SE -threads 6 -phred33 ')
                fout.write(ele +'.gz ' + ele+'.trim ')
                fout.write('ILLUMINACLIP:/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50 \n ')
                fout.write('/panfs/panfs.cluster/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i ' + ele +'.trim -o ' + ele + '.fa -r -Q33 \n' )
                fout.write('rm ' + ele +'.gz\n')
                fout.write('rm ' + ele +'.trim\n')
            fout.write('cat ' + _srr+'_2.fastq.fa >> ' + _srr+'_1.fastq.fa && rm ' + _srr+'_2.fastq.fa\n')
            fout.write('cat ' + _srr+'.fastq.fa >> ' + _srr+'_1.fastq.fa && rm ' + _srr+'.fastq.fa\n')
            fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/util/insilico_read_normalization.pl --seqType fa --JM 30G --max_cov 250 --output ./Trinity --CPU 12 --single '  + _srr+'_1.fastq.fa\n')
            fout.write('mv ./Trinity/'+ _srr+'_1.fastq.fa.normalized_K25_C250_pctSD200.fa ../../RNA/'+_srr+'.fa\n')
            fout.write('rm -rf ../'+_srr+'\n')
            fout.close()



def checkInchwormResult20170220():
    f = r"F:\Insects\Drosophila_melanogaster\2017GeneModeling\inchwormK27L31.fa"
    from Bio import SeqIO
    fout = f+'2'
    fout = open(fout,'w')
    for seq in SeqIO.parse(open(f),'fasta'):
        fout.write('>'+seq.description+'\n'+str(seq.seq)+'\n')
    fout.close()
    
    fo = open(f+'2')
    fout = open(f+'min53','w',1000000)
    line1 = fo.readline()
    line2 = fo.readline()
    while line1:
        if len(line2) >53:
            fout.write(line1)
            fout.write(line2)
        line1 = fo.readline()
        line2 = fo.readline()
    fout.close()
    
    #lengths
    fo = open(f+'2')
    import numpy as np
    l = np.array([int(ele[:-1].split()[-1]) for ele in fo if ele[0] == '>'])
    
    from collections import Counter
    lcount = Counter(l)
    fout = open('list3.txt','w')
    for key,keycount in lcount.items():
        fout.write('%d\t%d\n'%(key,keycount))
    fout.close()

def reads2transcripts20170221():
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

def readsSplitBasedOnalignment20170221():
    '''
    files generated from reads2transcripts20170221, split to small files 
    '''
    import heapq
    
    import os
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/reads/'
    files = os.listdir(folder)
    f = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/'
    
    fos = [open(folder+file,'r',80000000) for file in files if file[-4:] == 'sort']
    fout = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/monitor.txt','a')
    lines = heapq.merge(*fos,key = lambda x:int(x.split()[0]))
    
    start = 43639
    #start = 2
    for line in lines:
        eles = line.split()
        m = int(eles[0])
        if m >= start:
            break
    
    f2 = open(f+str(m),'w',80000000)
    seq = '>'+'s\n' + eles[-1]+'\n'
    f2.write(seq)
    n = m
    count = 1
    
    for line in lines:
        eles = line.split()
        m = int(eles[0])
        seq = '>'+'s\n' + eles[-1]+'\n'
        if m == n:
            f2.write(seq)
            count += 1
        else:
            f2.close()
            fout.write(str(n)+'\t'+str(count)+'\n')
            count = 1
            f2 = open(f+str(m),'w',80000000)
            f2.write(seq)
            n = m
    
    #        if n >10:
    #            break
    f2.close()
    fout.write(str(n)+'\t'+str(count)+'\n')
    fout.close()


def collectSomeInformation20170222():
    folder = "D:\\P\\3Language\\Xiaolong\\python\\temp\\"
    import os
    files = os.listdir(folder)
    fout = open('list.txt','w')
    for file in files:
        c = open(folder+file).read()
        fout.write(file.split('.')[0]+'\t'+c)
    fout.close()

def generateQsubForTrinityFinal20170222():
    seqs = range(40000,50000)
    fout = open('recursive_trinity.cmds','w')
    for seqn in seqs:
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/trinityout{0} --CPU 1 --max_memory 2G --seqType fa --trinity_complete --full_cleanup  \n'.format(seqn))
    fout.close()
    
    seqs = range(40000,50000)
    fout = open('recursive_trinity.cmds','w')
    for seqn in seqs:
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/trinityout{0} --CPU 12 --max_memory 30G --seqType fa --trinity_complete --full_cleanup --KMER_SIZE 27 &&\n'.format(seqn))
        fout.write('cat {0} >>/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/monitorTrinity.txt && \n'.format(seqn))
        fout.write('mv  /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{0}  /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqsUsed/{0} \n'.format(seqn))
    fout.close()
    
    
    seqs = range(30000,40000)
    fout = open('toremove.txt','w')
    for seqn in seqs:
        fout.write('rm /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{0}\n'.format(seqn))
    fout.close()
    

def collectFasta20170222():
    import os
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/'
    fout = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/DmTrinity.fa','a')
    from Bio import SeqIO
    files = os.listdir(folder)
    for file in files:
        if file[-6:] == '.fasta':
            head = file.split('trinityout')[1].split('.')[0]
            ls = list(SeqIO.parse(open(folder+file),'fasta'))
            for seq in ls:
                fout.write('>Dm'+head+seq.id+'\n'+str(seq.seq)+'\n')
    fout.close()

def convertReads2component2ReadsFiles20170226():
    '''
    kmer length 27, transcripts cannot be assembled properly, many transcripts are too short. 
    decide to run with kmer 25. But I've deleted the .fasta files. convert the reads2component file back to reads files.
    '''
    ls = open('list.txt').read().split()
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/reads/'
    folderout = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/RNA/'
    for ele in ls:
        filename = ele.split('.')[0]
        fout = open(filename,'w')
        fout.write('#!/panfs/panfs.cluster/opt/python/3.5.0/gcc/bin/python3\n')
        fout.write('#PBS -q batch\n')
        fout.write('#PBS -N %s\n'%filename)
        fout.write('#PBS -l nodes=1:ppn=1\n')
        fout.write('#PBS -l walltime=120:00:00\n')
        fout.write('#PBS -m abe -M atps@outlook.com\n')
        fout.write('#PBS -j oe\n')
        fout.write('''
    fo = open('%s%s')
    fout = open('%s%s.fa','w')
    for line in fo:
        fout.write('>s\\n'+line.split()[-1]+'\\n')
    fout.close()
    fo.close()'''%(folder,ele,folderout,filename))
        fout.close()

def generateCursiveTrinityBasedOnReadsNumber20170301():
    ls = open('monitor.txt').readlines()
    fout1 = open('recursive_trinity.txt','w')

    for ele in ls:
        n,c = ele.split()
        n = int(n)
        c = int(c)
        f = n // 1000
        if c > 10000000:
            fout2 = open('20170301qsub'+str(n)+'.txt','w')
            fout2.write('''
    #!/bin/bash
    #PBS -q batch
    #PBS -N ones
    #PBS -l nodes=1:ppn=12
    #PBS -l walltime=120:00:00
    #PBS -m abe -M atps@outlook.com
    #PBS -j oe
    
    module reset
    module load python
    module load perl
    module load gcc-4.9.2
    module load jdk/1.8.0_45
    module load bowtie2\n''')
            fout2.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{1}/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/trinityout{0} --CPU 12 --max_memory 30G --seqType fa --trinity_complete --full_cleanup  \n'.format(n,f))
            fout2.close()
        else:
            fout1.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/seqs/{1}/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/Dme/TrinityOne/trinityout{0} --CPU 1 --max_memory 2G --seqType fa --trinity_complete --full_cleanup  \n'.format(n,f))
    fout1.close()
    fout2.close()