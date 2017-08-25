# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 20:15:03 2016

@author: k
"""
def AgDownloadSTARCufflinksAlltogether():
    mylist = open("list.txt").readlines()
    dicNameLink ={}
    for ele in mylist:
        filename = ele.split("/")[-1]
        if "_" in filename:
            errname = filename.split("_")[0]
        else:
            errname = filename.split(".")[0]
        if errname not in dicNameLink:
            dicNameLink[errname] = []
        dicNameLink[errname].append(ele[:-1])
    
    #dicNameLink['ERR537787']
    
    for ele in dicNameLink:
        fo = open("Ag20160308_"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
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
        fo.write("module load jdk\n")
        fo.write("module load samtools\n")
        
        if len(dicNameLink[ele]) == 1:
            
            fo.write("cd RNA\n")
            fo.write("wget -q " + dicNameLink[ele][0] +"\n")
            fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar SE -threads 12 " +\
            ele+".fastq.gz " + ele+".fastqTrim.gz" +\
            " ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.0/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30 HEADCROP:15 LEADING:20 TRAILING:20 MINLEN:50 ")
            fo.write(" && mkdir ../Cuff/"+ele+" && rm "+ele+".fastq.gz "+"\n\n")
            
            fo.write("/home/ks2073/p/2016/STAR/bin/Linux_x86_64_static/STAR  --runThreadN 12  --runMode alignReads ")
            fo.write("--genomeDir /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/GenomeChrStar/STAR/ ")
            fo.write("--readFilesIn " + ele+".fastqTrim.gz --readFilesCommand gunzip -c ")
            fo.write("--outFileNamePrefix ../Cuff/"+ele+"/ --outSAMstrandField intronMotif --outSAMtype BAM Unsorted ")
            fo.write("--outReadsUnmapped fastx \n")
            
            
        else:
            fo.write("cd RNA\n")
            fo.write("wget -q " + dicNameLink[ele][0] +"\n")
            fo.write("wget -q " + dicNameLink[ele][1] +"\n")
            fo.write("java -jar /home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar PE -phred33 -threads 12 ")
            fo.write(ele+"_1.fastq.gz "+ele+"_2.fastq.gz "+ele+"_1fastqTrimPaired.fq.gz "+ele+"_1fastqTrimUnpaired.fq.gz ")
            fo.write(ele+"_2fastqTrimPaired.fq.gz "+ele+"_2fastqTrimUnpaired.fq.gz ")
            fo.write("ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30 HEADCROP:15 LEADING:20 TRAILING:20 MINLEN:50")
            fo.write(" && mkdir ../Cuff/"+ele+" && rm "+ele+"_1.fastq.gz "+" && rm "+ele+"_2.fastq.gz \n\n")
                    
            fo.write("/home/ks2073/p/2016/STAR/bin/Linux_x86_64_static/STAR  --runThreadN 12  --runMode alignReads ")
            fo.write("--genomeDir /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/GenomeChrStar/STAR/ ")
            fo.write("--readFilesIn " + ele+"_1fastqTrimPaired.fq.gz ")
            fo.write(ele+"_2fastqTrimPaired.fq.gz ")
            fo.write(" --readFilesCommand gunzip -c ")
            fo.write("--outFileNamePrefix ../Cuff/"+ele+"/ --outSAMstrandField intronMotif --outSAMtype BAM Unsorted ")
            fo.write("--outReadsUnmapped fastx \n\n")
            
        fo.write("samtools sort -m 2G -@ 12 -O bam -T ../Cuff/"+ ele +"/AgsortTemp -o ../Cuff/"+ele+"/Aligned.sort.bam ../Cuff/" + ele+"/Aligned.out.bam ")
        fo.write("&& rm ../Cuff/"+ele+"/Aligned.out.bam \n\n")
        fo.write("cufflinks -q -p 12 -b /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/GenomeChrStar/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa ")
        fo.write("-o ../Cuff/"+ele +"/ ../Cuff/"+ele+"/Aligned.sort.bam\n\n")
        
        fo.write("cp ../Cuff/"+ele +"/transcripts.gtf ../Cuff_gtf/"+ele+".gtf\n")
        fo.close()
    

def rsemAg():
    mylist = open("list.txt").readlines()
    dicNameLink ={}
    for ele in mylist:
        filename = ele.split("/")[-1]
        if "_" in filename:
            errname = filename.split("_")[0]
        else:
            errname = filename.split(".")[0]
        if errname not in dicNameLink:
            dicNameLink[errname] = []
        dicNameLink[errname].append(ele[:-1])

    for ele in dicNameLink:
        fo = open("Ag20160420RSEM_"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/RSEM/RSEMTemp\n")
        fo.write("module load gcc-4.9.2\n")
        fo.write("module load python\n")
        fo.write("module load jdk\n")
        fo.write("module load perl\n")
        fo.write("module load bowtie2\n")
        fo.write("mkdir "+ele+"\n")
        fo.write("cd " + ele +"\n")
        
        if len(dicNameLink[ele]) == 1:
        
            fo.write("/home/ks2073/p/pigz-2.3.3/pigz -d -p 12 " + dicNameLink[ele][0] +"\n")
            fileFQ = dicNameLink[ele][0].rsplit(".",1)[0]
            fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/RSEM-1.2.29/rsem-calculate-expression --no-qualities -p 12 --bowtie2 --no-bam-output %s /panfs/panfs.cluster/scratch/ks2073/RSEM/RSEMTemp/DNAseq/ProDB %s\n"%(fileFQ,ele))
            fo.write("cp *.isoforms.results ../result/\n")
            fo.write("/home/ks2073/p/pigz-2.3.3/pigz -p 12  " + fileFQ +"\n")
            
        else:
            fo.write("/home/ks2073/p/pigz-2.3.3/pigz -d -p 12 " + dicNameLink[ele][0] +"\n")
            fo.write("/home/ks2073/p/pigz-2.3.3/pigz -d -p 12 " + dicNameLink[ele][1] +"\n")
            fileFQ1 = dicNameLink[ele][0].rsplit(".",1)[0]
            fileFQ2 = dicNameLink[ele][1].rsplit(".",1)[0]
            fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/RSEM-1.2.29/rsem-calculate-expression --no-qualities -p 12 --bowtie2 --no-bam-output %s,%s /panfs/panfs.cluster/scratch/ks2073/RSEM/RSEMTemp/DNAseq/ProDB %s\n"%(fileFQ1,fileFQ2,ele))
            fo.write("cp *.isoforms.results ../result/\n")
            fo.write("/home/ks2073/p/pigz-2.3.3/pigz -p 12  " + fileFQ1 +"\n")
            fo.write("/home/ks2073/p/pigz-2.3.3/pigz -p 12  " + fileFQ2 +"\n")
            
        fo.close()

def decompressNormalizeAg():
    """
    use khmer to normalize reads. For New assembly by Trinity or Bridger
    """
    lsFq = open("list.txt").read().split()
    for key in lsFq:
        ele = key.split("/")[-1].split(".")[0]
        fo = open("20160324AgNormal"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load gcc-4.9.2\n")
        fo.write("module load python\n")
        fo.write("module load jdk\n")
        fo.write("module load perl\n")
        fo.write("gzip -d "+key+"\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA\n")
        fo.write("/panfs/panfs.cluster/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta -i "+key[:-3] +" -Q 33 -r -o "+key[:-3].split("fastq")[0]+".fa\n")
#        fo.write("/panfs/panfs.cluster/opt/khmer/gcc/1.4.1/scripts/normalize-by-median.py -k 27 ")
#        fo.write(" -o /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/"+ele+"Norm ")
#        fo.write(" -N 8 -C 100 -M 2e10 "+key[:-3]+"\n")
        fo.close()
        
def compressAg():
    """
    use khmer to normalize reads. For New assembly by Trinity or Bridger
    """
    lsFq = open("list.txt").read().split()
    for key in lsFq:
        ele = key.split("/")[-1].split(".")[0]
        fo = open("20160324AgNormal"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=10:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load gcc-4.9.2\n")
        fo.write("module load python\n")
        fo.write("module load jdk\n")
        fo.write("module load perl\n")
        fo.write("/home/ks2073/p/pigz-2.3.3/pigz -p 12  " +key+"\n")
        fo.close()


def AgNcbiDownload():
    """
    download based srr numbers
    """
    lsSRR = open("list.txt").read().split()
    for ele in lsSRR:
        fo = open("20160323NCBI_srr_download"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/\n")
        fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/sratoolkit.2.5.7-ubuntu64/bin/fastq-dump --split-files " +ele +"\n")
        fo.close()

def AgNcbiDownloadFastQC():
    """
    download based srr numbers calculate fastQC
    """
    lsSRR = open("list.txt").read().split()
    for ele in lsSRR:
        fo = open("20160323NCBI_srr_FastQC"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/\n")
        fo.write("/home/ks2073/p/FastQC0.11.3/fastqc " +ele +".fastq -o ../fastQC/\n")
        fo.close()

def AgNCBIFastQTrim():
    mylist = open("list.txt").read().split()
    for ele in mylist:
        fo = open("Ag20160323RSEM_"+ele.split("/")[-1]+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N "+ele.split("/")[-1]+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load gcc-4.9.2\n")
        fo.write("module load python\n")
        fo.write("module load jdk\n")
        fo.write("module load perl\n")
        fo.write("module load bowtie2\n")
        
        fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar SE -phred33 -threads 12 " +\
        ele+" " + ele+".fastqTrim" +\
        " ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30  LEADING:20 TRAILING:20 MINLEN:50 \n")
        fo.close()


def AgRSEMNCBi():
    """
    RSEM for reads downloaded from NCBI. Need to trim those reads before use
    """
    mylist = open("list.txt").readlines()
    dicNameLink ={}
    for ele in mylist:
        filename = ele.split("/")[-1]
        if "_" in filename:
            errname = filename.split("_")[0]
        else:
            errname = filename.split(".")[0]
        if errname not in dicNameLink:
            dicNameLink[errname] = []
        dicNameLink[errname].append(ele[:-1])
    
        for ele in dicNameLink:
            fo = open("Ag20160323RSEM_"+ele+".txt","w")
            fo.write("#!/bin/bash\n")
            fo.write("#PBS -q batch\n")
            fo.write("#PBS -N "+ele+"\n")
            fo.write("#PBS -l nodes=1:ppn=12\n")
            fo.write("#PBS -l walltime=120:00:00\n")
            fo.write("#PBS -m abe -M atps@outlook.com\n")
            fo.write("#PBS -j oe\n")
            fo.write("cd /panfs/panfs.cluster/scratch/ks2073/RSEM/RSEMTemp\n")
            fo.write("module load gcc-4.9.2\n")
            fo.write("module load python\n")
            fo.write("module load jdk\n")
            fo.write("module load perl\n")
            fo.write("module load bowtie2\n")
            fo.write("mkdir "+ele+"\n")
            fo.write("cd " + ele +"\n")
            
            if len(dicNameLink[ele]) == 1:
                fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar SE -phred33 -threads 12 /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/" +\
                ele+".fastq " + ele+".fastqTrim" +\
                " ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:30 \n")

                fileFQ = ele+".fastqTrim"
                fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/RSEM-1.2.29/rsem-calculate-expression -p 12 --bowtie2 --no-bam-output %s /panfs/panfs.cluster/scratch/ks2073/RSEM/RSEMTemp/DNAseq/ProDB %s\n"%(fileFQ,ele))
                
            else:
                fo.write("java -jar /home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar PE -phred33 -threads 12 /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/")
                fo.write(ele+"_1.fastq /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/"+ele+"_2.fastq "+ele+"_1fastqTrimPaired "+ele+"_1fastqTrimUnpaired ")
                fo.write(ele+"_2fastqTrimPaired "+ele+"_2fastqTrimUnpaired ")
                fo.write("ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:30\n")


                fileFQ1 = ele+"_1fastqTrimPaired"
                fileFQ2 = ele+"_2fastqTrimPaired"
                fileFQ3 = ele+"_1fastqTrimUnpaired"
                fileFQ4 = ele+"_2fastqTrimUnpaired"
                fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/RSEM-1.2.29/rsem-calculate-expression -p 12 --bowtie2 --no-bam-output %s,%s,%s,%s /panfs/panfs.cluster/scratch/ks2073/RSEM/RSEMTemp/DNAseq/ProDB %s\n"%(fileFQ1,fileFQ2,fileFQ3,fileFQ4,ele))
            
            
            fo.write("cp *.isoforms.results ../result/\n")
                
            fo.close()

def changeFormatAg():
    lsFq = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/list.txt").read().split()
    fout = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/allRNA.fa","w")
    
    def baseEncode(number, base=36):
        """
        given a int number, return a string based on base
        10:A, 36:Z, 37:Z1, under default setting
        """
        if not isinstance(number, int):
            raise TypeError('number must be an integer')
        alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
        if base > 62 or base <=1:
            print("base should be between 2 and 62")
            return None
        sign = ""
        if number < 0:
            sign = "-"
            number = -number
        alphabet = alphabet[:base+1]
        if 0 <= number and number <base:
            return sign+alphabet[number]
        numberbase=""
        while number != 0:
            number, i = divmod(number, base)
            numberbase = alphabet[i] + numberbase
        return sign+numberbase
    
    for file in lsFq:
        fo = open(file,"r")
        mylist = fo.readlines()
        fo.close()
        for num in range(0,len(mylist),4):
            head = "1"
            fout.write(">"+head+"\n"+mylist[num+1])
        del mylist
#            if num >1000:
#                break
#        break
    fout.close()
            
def changeFastq2fastaAndRenameToNumbers():
    mylist = open("list.txt").readlines()
    dicNameLink ={}
    for ele in mylist:
        filename = ele.split("/")[-1]
        errname = filename.split(".")[0]
        if errname not in dicNameLink:
            dicNameLink[errname] = []
        dicNameLink[errname].append(ele[:-1])
    textcode = '''
    def baseEncode(number, base=36):
        """
        given a int number, return a string based on base
        10:A, 36:Z, 37:Z1, under default setting
        """
        if not isinstance(number, int):
            raise TypeError('number must be an integer')
        alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
        if base > 62 or base <=1:
            print("base should be between 2 and 62")
            return None
        sign = ""
        if number < 0:
            sign = "-"
            number = -number
        alphabet = alphabet[:base+1]
        if 0 <= number and number <base:
            return sign+alphabet[number]
        numberbase=""
        while number != 0:
            number, i = divmod(number, base)
            numberbase = alphabet[i] + numberbase
        return sign+numberbase
    
    def changeFormat(filename,outname):
        fo = open(filename,"r",2000000000)
        fout = open(outname,"w",2000000000)
        line1 = fo.readline()
        line2 = fo.readline()
        line3 = fo.readline()
        line4 = fo.readline()
        count = 0
        while line2:
            count += 1
            fout.write(">"+baseEncode(count)+"\\n"+line2)
            line1 = fo.readline()
            line2 = fo.readline()
            line3 = fo.readline()
            line4 = fo.readline()
        fout.close()
        fo.close()
    '''


    for ele in dicNameLink:
        fo = open("Ag20160419_"+ele+".txt","w")
        fo.write("#!/usr/bin/env python\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=3\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("import os\n")
        fo.write(textcode)
        fo.write("os.chdir('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA/')\n")
        fo.write("changeFormat('%s','%s.fa')\n"%(dicNameLink[ele][0],ele))
        fo.close()
