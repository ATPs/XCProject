# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 17:11:01 2015

@author: k
"""

filelist = open("list.txt").readlines() #list.txt stores the full path to files.
filenames =[]
folder = filelist[0][:filelist[0].rfind("/")+1]
folder2 =folder[:-1]+"d/"

for ele in filelist:
    filenames.append(ele.split("/")[-1][:-1])

for ele in filenames:
    fo = open("CufflinksTopHat"+ele+".txt", "w")
    fo.write("#!/bin/bash\n")
    fo.write("#PBS -q phi\n")
    fo.write("#PBS -N CT" +ele+"\n")
    fo.write("#PBS -l nodes=1:ppn=24\n")
    fo.write("#PBS -l walltime=120:00:00\n")
    fo.write("#PBS -m abe -M atps@outlook.com\n")
    fo.write("#PBS -j oe\n")
    fo.write("cd /scratch/ks2073/Transcriptome/AnophelesGambiae\n")
    fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/samtools-0.1.19/\n")
    fo.write("export PATH\n")
    fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/bowtie2-2.2.3/\n")
    fo.write("export PATH\n")
    fo.write("module load tophat/2.0.12\n")
    fo.write("module load cufflinks/2.2.1\n")
    fo.write("tophat -p 48 --library-type fr-firststrand --read-realign-edit-dist 0 -o ./Cuff/"+ele+"/ ./OGS/AgOGS "+folder+ele+"\n")
    fo.write("mv "+folder+ele+" "+folder2+"\n")
    fo.write("cufflinks -q -p 48 -b ./OGS/AgOGS.fa -u -o ./Cuff/"+ele+"/ ./Cuff/"+ele+"/accepted_hits.bam\n")
    fo.close()







#for running with phi and STAR

filelist = open("list.txt").readlines()
filenames =[]
folder = filelist[0][:filelist[0].rfind("/")+1]
folder2 =folder[:-1]+"d/"

for ele in filelist:
    filenames.append(ele.split("/")[-1][:-1])

for ele in filenames:
    fo = open("CufflinksTopHat"+ele+".txt", "w")
    fo.write("#!/bin/bash\n")
    fo.write("#PBS -q phi\n")
    fo.write("#PBS -N CT" +ele+"\n")
    fo.write("#PBS -l nodes=1:ppn=24\n")
    fo.write("#PBS -l walltime=120:00:00\n")
    fo.write("#PBS -m abe -M atps@outlook.com\n")
    fo.write("#PBS -j oe\n")
    fo.write("cd /scratch/ks2073/Transcriptome/AnophelesGambiae\n")
    fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/samtools-0.1.19/\n")
    fo.write("export PATH\n")
    fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/bowtie2-2.2.3/\n")
    fo.write("export PATH\n")
    fo.write("module load tophat/2.0.12\n")
    fo.write("module load cufflinks/2.2.1\n")
    fo.write("mv "+folder+ele+" "+folder2+"\n")
    fo.write("mkdir ./Cuff/"+ele + "/\n")
    fo.write("/home/ks2073/p/2015/STAR/bin/Linux_x86_64/STAR --runThreadN 48 --runMode alignReads --readFilesIn "+folder2+ele\
    + " --genomeDir /scratch/ks2073/Transcriptome/AnophelesGambiae/OGS_STAR --outReadsUnmapped Fastx --outSAMtype  BAM   Unsorted   --outSAMstrandField intronMotif "\
    + " --outFileNamePrefix ./Cuff/"+ele + "/\n")
    fo.write("samtools sort -m 4G -o ./Cuff/"+ele + "/accepted_hits.bam -O bam -@ 24 -T ./Cuff/"+ele + "/_STARtmp/ " +\
    "./Cuff/"+ele + "/Aligned.out.bam\n")
    fo.write("rm ./Cuff/"+ele + "/Aligned.out.bam\n")
    fo.write("rm "+folder2+ele+"\n")
    fo.write("cufflinks -q -p 48 -u -o ./Cuff/"+ele+"/ ./Cuff/"+ele+"/accepted_hits.bam\n")
    fo.write("rm ./Cuff/"+ele+"/accepted_hits.bam\n")
    fo.close()



