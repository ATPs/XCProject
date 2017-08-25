"""
for paired end reads
"""

mylist = open("list.txt").readlines()
mynames = set()
for ele in mylist:
    mynames.add(ele.split("/")[-1].split("_")[0])

print(len(mynames))

txt1="""#!/bin/bash

#PBS -q bigmem
#PBS -N Trimatic
#name you want to give your job 
#PBS -l nodes=1:ppn=1
#request 1 nodes w/8 processors per node  
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

cd /scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq

"""
for ele in mynames:
    fo = open("Trim_"+ele+".txt","w")
    fo.write(txt1)
    fo.write("java -jar /home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/trimmomatic.jar PE -phred33 -threads 12 ")
    fo.write(ele+"_1.fastq "+ele+"_2.fastq "+ele+"_1fastqTrimPaired.fq "+ele+"_1fastqTrimUnpaired.fq ")
    fo.write(ele+"_2fastqTrimPaired.fq "+ele+"_2fastqTrimUnpaired.fq ")
    fo.write("ILLUMINACLIP:/home/ks2073/p/2015/trinityrnaseq-2.1.1/trinity-plugins/Trimmomatic/adapters/XCTruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30 HEADCROP:15 LEADING:20 TRAILING:20 MINLEN:50")
    fo.write("\n")
    fo.write("/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta  -Q33 -i ")
    fo.write(ele+"_1fastqTrimPaired.fq -o ../RNA/"+ele+"_1.fa\n")
    fo.write("/home/ks2073/p/fastx_toolkit_0.0.13/fastq_to_fasta  -Q33 -i ")
    fo.write(ele+"_2fastqTrimPaired.fq -o ../RNA/"+ele+"_2.fa\n")
    fo.close()




"""
create TopHat Cufflinks commard for paired end reads. Only use those paired reads. 
"""
for ele in mynames:
    fo = open("CufflinksTopHat"+ele+".txt","w")
    fo.write("#!/bin/bash\n")
    fo.write("#PBS -q phi\n")
    fo.write("#PBS -N CT"+ele+"\n")
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
    fo.write("tophat -p 48 --read-realign-edit-dist 0 -o ./Cuff/"+ele+" ./OGS/AgOGS ./RNA/"+ele+"_1.fa ./RNA/"+ele+"_2.fa\n")
    fo.write("cufflinks -q -p 48 -b ./OGS/AgOGS.fa -u -o ./Cuff/"+ele+" ./Cuff/"+ele+"/accepted_hits.bam\n")
    fo.write("cp ./Cuff/"+ele+"/transcripts.gtf ./Cuff_gtf/"+ele+".gtf\n")
    fo.write("mv ./RNA/"+ele+"_1.fa ./RNAd/\n")
    fo.write("mv ./RNA/"+ele+"_2.fa ./RNAd/\n")
    fo.close()
    

    