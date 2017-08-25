# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 11:17:55 2015

@author: k
"""

txt1="""#!/bin/bash

#PBS -q batch
#PBS -N RSEM
#name you want to give your job 
#PBS -l nodes=1:ppn=12
#request 1 nodes w/8 processors per node
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/bowtie2-2.2.3/ 
export PATH
PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/samtools-0.1.19/
export PATH
cd /scratch/ks2073/RSEM/RSEMTemp

module load gcc-4.7.2

"""
mylist = open("list.txt","r").read().split()
for ele in mylist:
    fo = open(ele+".pbs","w")
    fo.write(txt1)
    fo.write("mkdir /scratch/ks2073/RSEM/RSEMTemp/out"+ele+"\n")
    fo.write("/home/ks2073/p/rsem-1.2.15/rsem-calculate-expression  --no-qualities --bowtie2   --samtools-sort-mem 2G -p 12 --no-bam-output  /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA/ERR"+ele+"_1.fastqTrim.fa50,/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA/ERR"+ele+"_2.fastqTrim.fa50 /scratch/ks2073/RSEM/RSEMTemp/DNAseq/ProDB /scratch/ks2073/RSEM/RSEMTemp/out"+ele+"/"+ele)
    fo.close()


txt1="""#!/bin/bash

#PBS -q batch
#PBS -N RSEM3
#name you want to give your job 
#PBS -l nodes=1:ppn=12
#request 1 nodes w/8 processors per node
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/bowtie2-2.2.3/ 
export PATH
PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/samtools-0.1.19/
export PATH
cd /scratch/ks2073/RSEM/RSEMTemp3

module load gcc-4.7.2

"""


mylist = open("list.txt","r").read().split()
for ele in mylist:
    fo = open(ele+".pbs","w")
    fo.write(txt1)
    fo.write("mkdir /scratch/ks2073/RSEM/RSEMTemp3/out"+ele+"\n")
    fo.write("/home/ks2073/p/rsem-1.2.15/rsem-calculate-expression  --no-qualities --bowtie2   --samtools-sort-mem 2G -p 12 --no-bam-output  /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA/ERR"+ele+".fastqTrim.fa50 /scratch/ks2073/RSEM/RSEMTemp3/DNAseq/ProDB /scratch/ks2073/RSEM/RSEMTemp3/out"+ele+"/"+ele)
    fo.close()

    