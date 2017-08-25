# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:07:05 2015

@author: k
"""

txt1 = """#!/bin/bash

#PBS -q batch
#PBS -N Trim
#name you want to give your job 
#PBS -l nodes=1:ppn=1
#request 1 nodes w/8 processors per node  
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

cd /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq

"""

downloads = list(open("list.txt", "r"))
filenames =[]
for ele in downloads:
    filenames.append(ele.split("/")[-1][:-1].split(".")[0])
for num in range(len(filenames)):
    fo = open("PBS"+filenames[num]+".pbs","w")
    fo.write(txt1+"wget "+downloads[num] + "\n"+ \
    "gzip -d " + filenames[num]+".fastq.gz \n" +\
    "java -jar /panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.0/trinity-plugins/Trimmomatic/trimmomatic.jar SE -threads 12 " +\
    filenames[num]+".fastq " + filenames[num]+".fastqTrim" +\
    " ILLUMINACLIP:/panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.0/trinity-plugins/Trimmomatic/adapters/XCTruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50" +\
    "\nrm " +filenames[num]+".fastq \n" +\
    "#/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i " +filenames[num]+".fastqTrim -o " +\
    filenames[num]+".fastqTrim.fa -Q 33" +\
    "\n#rm " + filenames[num]+".fastqTrim\n")
    fo.close()


txt2 = """#!/bin/bash

#PBS -q batch
#PBS -N Trim
#name you want to give your job 
#PBS -l nodes=1:ppn=1
#request 1 nodes w/8 processors per node  
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

cd /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq

"""




downloads = list(open("list.txt", "r"))
filenames =[]
for ele in downloads:
    filenames.append(ele.split("/")[-1][:-1].split(".")[0])
writelist = []
for num in range(len(filenames)):
    writelist.append("wget "+downloads[num] + "\n"+ \
    "gzip -d " + filenames[num]+".fastq.gz \n" +\
    "java -jar /panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.0/trinity-plugins/Trimmomatic/trimmomatic.jar SE -threads 12 " +\
    filenames[num]+".fastq " + filenames[num]+".fastqTrim" +\
    " ILLUMINACLIP:/panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.0/trinity-plugins/Trimmomatic/adapters/XCTruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:50" +\
    "\nrm " +filenames[num]+".fastq \n" +\
    "mv " +filenames[num]+".fastqTrim ../RNA/\n\n")

dic_write = {}
for num in range(12):
    dic_write[num] = open("AgDownload"+str(num)+".txt","w")
    dic_write[num].write(txt1)


for num in range(116):
    dic_write[num % 12].write(writelist[num])

for num in range(12):
    dic_write[num].close()
    


"""
download and decompress only. each to a individual file. 
20160210
"""
txt3 = """#!/bin/bash

#PBS -q bigmem
#PBS -N Trim
#name you want to give your job 
#PBS -l nodes=1:ppn=1
#request 1 nodes w/8 processors per node  
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

cd /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq

"""
downloads = list(open("list.txt", "r")) #each line looks like: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR593/ERR593540/ERR593540.fastq.gz
filenames =[] # list of filename. here it should be like ERR593540
for ele in downloads:
    filenames.append(ele.split("/")[-1][:-1].split(".")[0])
writelist = [] #list of commands to work. here is: wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR593/ERR593540/ERR593540.fastq.gz && gzip -d ERR593540.fastq.gz
for num in range(len(filenames)):
    writelist.append("wget "+downloads[num][:-1] + " && "+ \
    "gzip -d " + filenames[num]+".fastq.gz \n")

for ele in range(len(filenames)):
    fo = open("download_"+filenames[ele]+".txt","w")
    fo.write(txt3+writelist[ele])
    fo.close()


def downloadInOneCommand():
    """
    download and decompress only. each to a individual file. 
    20160210
    """
    txt3 = """#!/bin/bash

#PBS -q bigmem
#PBS -N Trim
#name you want to give your job 
#PBS -l nodes=1:ppn=1
#request 1 nodes w/8 processors per node  
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

cd /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq

"""
    downloads = list(open("list.txt", "r")) #each line looks like: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR593/ERR593540/ERR593540.fastq.gz
    filenames =[] # list of filename. here it should be like ERR593540
    for ele in downloads:
        filenames.append(ele.split("/")[-1][:-1].split(".")[0])
    writelist = [] #list of commands to work. here is: wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR593/ERR593540/ERR593540.fastq.gz && gzip -d ERR593540.fastq.gz
    for num in range(len(filenames)):
        writelist.append("wget "+downloads[num][:-1] + " && "+ \
        "gzip -d " + filenames[num]+".fastq.gz ")
    
    fo = open("download.txt","w")
    fo.write(txt3)
    for i in range(0,len(writelist),8):
        for j in range(i,i+min(8,len(writelist)-i)):
            fo.write(writelist[j]+" & ")
        fo.write("\n")
    fo.close()
downloadInOneCommand()


"""
Trim based on quality score. Phred33
based on fastQC, the first ~10 is is not so good, trim the first 15 bases off for each read.
SLIDINGWINDOW:4:30 LEADING:20 TRAILING:20 MINLEN:50
slidingwindow size 4, mean quality score is 30.
leading min quality score 20, trailing min 20
fastq-to-fasta, do not keep sequences with unknown (N) nucleotides. rename sequence identifiers to numbers
"""

myfiles = open("list.txt", "r").read().split()

txt3 = """#!/bin/bash

#PBS -q phi
#PBS -N Trim
#name you want to give your job 
#PBS -l nodes=1:ppn=8
#request 1 nodes w/8 processors per node  
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

cd /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq

"""
for ele in myfiles:
    fo = open("trim_2fa_"+ele+".txt","w")
    fo.write(txt3)
    fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar SE -phred33 -threads 12 " + \
    ele+" " + ele+"Trim " +\
    "ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30 HEADCROP:15 LEADING:20 TRAILING:20 MINLEN:50" +\
    "\n")
    fo.write("/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -r -Q33 -i " +ele+"Trim -o ../RNA/" +ele+".fa\n")
    fo.close()



