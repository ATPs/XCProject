# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 17:08:50 2016

@author: k
"""

def csvopener(filename, separator = None):
    """
    given a filename of csv, return a list of dictionaries by open the csv file to dictionary
    """
    import unicodecsv
    with open(filename,"rb") as f:
        if separator is None:
            reader = unicodecsv.DictReader(f)
        else:
            reader = unicodecsv.DictReader(f, delimiter = separator)
        return list(reader)

def msPrepareFilenameListDc():
    '''
    order of files were stored in file 20160707librariesAndFileNamesInfo.txt, in the following format
    Num	Type	ShortName	File1	File2
    1	paired	H-L2-D1	G9-1f.txt	G9-1r.txt
    2	paired	H-L3-D1	G9-2f.txt	G9-2r.txt
    3	paired	H-L4-12h	G9-3f.txt	G9-3r.txt
    ...
    67	single	An-A-M-S	ERR1011037.fastq	
    return a list of dictionary using unicodecsv.DictReader
    '''
    filename = 'D:\\mine\\OneDrive\\Lab\\Jiang Lab Pubilic\\ManducaTranscriptome\\20160425ManducaTranscriptomeManuscript\\20160707librariesAndFileNamesInfo.txt'
    ls_libraries = csvopener(filename, '\t')
    print(len(ls_libraries))
        


def MsNcbiDownload():
    """
    download based srr numbers
    """
    lsSRR = open("list.txt").read().split()
    for ele in lsSRR:
        fo = open("20160608NCBI_srr_download"+ele+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N "+ele+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/\n")
        fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/sratoolkit.2.5.7-ubuntu64/bin/fastq-dump --split-files " +ele +"\n")
        fo.close()

def unmappedReads():
    f_unmappedReads = r'E:\store_for_D\Transcriptome\unMappedReads\ManduaUnmappedReadsUnique.fq'
    templist = open(f_unmappedReads).readlines()
    return templist
    
def trimReads20160607(ls_libraries):
    """
    trim all reads with Trimmomatic
    trim paired as paired, trim single end as single end.
    min length 50
    """
    for ele in ls_libraries:
        fo = open("20160711Ms_"+ele['ShortName']+".txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N "+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/RNA\n")
        fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/samtools-0.1.19/\n")
        fo.write("export PATH\n")
        fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/bowtie2-2.2.3/\n")
        fo.write("export PATH\n")
        fo.write("module load jdk\n")
        fo.write("module load samtools\n")
        
        if ele['File2'] == '':
            
            fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic.jar SE -phred33 -threads 12 " +\
            ele['File1']+" ../RNA-trim/" + ele['File1'] +\
            " ILLUMINACLIP:/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/SE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:10 TRAILING:10 MINLEN:50 ")
            
            
            
            
        else:
            fo.write("java -jar /panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/trimmomatic.jar PE -phred64 -threads 12 ")
            fo.write(ele['File1']+" "+ele['File2']+" ../RNA-trim/"+ele['File1']+"paired.fq "+" ../RNA-trim/"+ele['File1']+"unpaired.fq ")
            fo.write(" ../RNA-trim/"+ele['File2']+"paired.fq "+" ../RNA-trim/"+ele['File2']+"unpaired.fq ")
            fo.write("ILLUMINACLIP:/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-2.2.0/trinity-plugins/Trimmomatic-0.32/adapters/PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:10 TRAILING:10 MINLEN:50")
            
                    
            
       
        fo.write("\n")
        fo.close()


def msTopHatCufflinks(ls_libraries):
    """
    """
    for ele in ls_libraries:
        fo = open("20160712Ms_"+ele['ShortName']+"TopHatCufflinks.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N TopCuff"+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/RNA-trim/\n")
        fo.write("module load jdk\n")
        fo.write("module load samtools\n")
        fo.write("module load cufflinks\n")
        fo.write("module load bowtie2\n")
        fo.write("/panfs/panfs.cluster/opt/tophat/2.0.12/prebuilt/tophat -p 12 --read-realign-edit-dist 0 -o ../TopHat/"+ele['ShortName']+"/")
        fo.write(" /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/GenomeMs/msBowtie2DB ")
        if ele['File2'] == '':
            file_reads = ele['File1']
            fo.write(file_reads)
            
            
            
            
            
        else:
            file_reads_p1 = ele['File1']+"paired.fq"
            file_reads_p2 = ele['File2']+"paired.fq"
            file_reads_s1 = ele['File1']+"unpaired.fq"
            file_reads_s2 = ele['File2']+"unpaired.fq"
            fo.write(file_reads_p1 +" " + file_reads_p2+","+file_reads_s1+","+file_reads_s2)
                    
            
        fo.write("\ncufflinks -p 12 -b /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/GenomeMs/ms_scaffolds.fa -u -o ../Cufflinks/")
        fo.write(ele['ShortName']+"/ " + "../TopHat/"+ele['ShortName']+"/accepted_hits.bam\n")
        fo.write("\n")
        fo.close()

def msCuffQuant(ls_libraries):
    """
    """
    for ele in ls_libraries:
        fo = open("20160713Ms_"+ele['ShortName']+"CuffQuant.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N CuffQuant"+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/Cuffquant/\n")
        fo.write("module load jdk\n")
        fo.write("module load samtools\n")
        fo.write("module load cufflinks\n")
        fo.write("module load bowtie2\n")
        fo.write("cuffquant -p 12 -b /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/GenomeMs/ms_scaffolds.fa -u -o ")
        fo.write(ele['ShortName']+"/ /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/gff/20160713ManducaSexta_merged.gtf " + "../TopHat/"+ele['ShortName']+"/accepted_hits.bam\n")
        fo.write("\n")
        fo.close()



def msBamIndexStat20160719(ls_libraries):
    """
    """
    for ele in ls_libraries:
        fo = open("20160712Ms_"+ele['ShortName']+"bamIndexStat.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N TopCuff"+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/bamStat/\n")
        fo.write("module load jdk\n")
        fo.write("module load samtools\n")
        fo.write("module load cufflinks\n")
        fo.write("module load bowtie2\n")
        
                    
            
        fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/samtools-1.3/samtools index ")
        fo.write("../TopHat/"+ele['ShortName']+"/accepted_hits.bam\n")
        fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/samtools-1.3/samtools idxstats ")
        fo.write("../TopHat/"+ele['ShortName']+"/accepted_hits.bam >>"+ele['ShortName']+".SamIndexStat.txt\n")
        fo.write("\n")
        fo.close()


def msBamIndexStatDataDealWithCoverage20160719(ls_libraries):
    """
    combine the result to a single file.
    """
    folder = "E:\\Lab\\works\\20121023_manduca_genome_new\\20120818Genome Ms\\scaffold\\20160719coverage\\"
    import csv
    dc_files = {}
    for ele in ls_libraries:
        dc_files[ele['ShortName']] = list(csv.reader(open(folder + ele['ShortName']+".SamIndexStat.txt",'r'),delimiter = '\t'))
    fout = open("20160719MsGenomeContigReadsSummary.txt", "w")
    fout.write('contig\tlength\t')
    for ele2 in ls_libraries:
        fout.write(ele2['ShortName']+'\t')
    fout.write('\n')
    for num in range(len(dc_files[ele['ShortName']])):
        fout.write(dc_files[ele['ShortName']][num][0] + '\t' + dc_files[ele['ShortName']][num][1] + '\t')
        for ele2 in ls_libraries:
            fout.write(dc_files[ele2['ShortName']][num][2] + '\t')
        fout.write("\n")
    fout.close()

def msRNATrimSurvivingRate20160719(ls_libraries):
    """
    get the overall surviving rate by Trimmomatic
    """
    folder = "E:\\store_for_D\\RNAlab\\20160711TrimRecord\\"
    fout = open("20160719TrimmomaticSurrivingRateMs.txt","w")
    import glob
    for ele in ls_libraries:
        filename = glob.glob(folder + ele['ShortName'] +".o*")
        filename = filename[0]
        ls_content = open(filename).readlines()
        if ele['Type'] == 'paired':
            line_survive = ls_content[13]
            line_survive_ele = line_survive.split()
            reads_total = int(line_survive_ele[3]) *2
            reads_pair = int(line_survive_ele[6]) *2
            reads_f = int(line_survive_ele[11])
            reads_r = int(line_survive_ele[16])
            overall_rate = (reads_f + reads_r + reads_pair)/float(reads_total)
        else:
            line_survive = ls_content[7]
            line_survive_ele = line_survive.split()
            overall_rate = int(line_survive_ele[4])/float(line_survive_ele[2])
        fout.write(ele["ShortName"]+"\t"+str(overall_rate)+"\n")
    
    fout.close()


def msRNATopHatSurvivingRate20160719(ls_libraries):
    """
    get the overall surviving rate by TopHat
    """
    folder = "E:\\store_for_D\\RNAlab\\20160712TopHatRecord\\"
    fout = open("20160719TopHatSurrivingRateMs.txt","w")
    import glob
    for ele in ls_libraries:
        filename = glob.glob(folder + ele['ShortName'] +".txt")
        filename = filename[0]
        ls_content = open(filename).readlines()
        if ele['Type'] == 'paired':
            line_survive = ls_content[12]
            line_survive_ele = line_survive.split()
            overall_rate = line_survive_ele[0]
        else:
            line_survive = ls_content[4]
            line_survive_ele = line_survive.split()
            overall_rate = line_survive_ele[0]
        fout.write(ele["ShortName"]+"\t"+str(overall_rate)+"\n")
    
    fout.close()


def msRNASequencingDepthEachposition20160720(ls_libraries):
    """
    generate a scripts for get the sequencing depth of each base in the genome. 
    """
    fo = open("msRNASequencingDepthEachposition20160720.txt",'w')
    fo.write("#!/bin/bash\n")
    fo.write("#PBS -q batch\n")
    fo.write("#PBS -N depthSamtools\n")
    fo.write("#PBS -l nodes=1:ppn=12\n")
    fo.write("#PBS -l walltime=120:00:00\n")
    fo.write("#PBS -m abe -M atps@outlook.com\n")
    fo.write("#PBS -j oe\n")
    fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/TopHat\n")
    fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/samtools-1.3/samtools depth -a -a -m 100000000 ")
    for ele in ls_libraries:
        fo.write(' ./'+ele['ShortName']+"/accepted_hits.bam  ")
    fo.write(">>../20160720sequencingDepth.txt\n")
    fo.close()


def msGenomeDescription20160721(ls_libraries):
    """
    get some information about the newest genome.
    """
    filename = r"E:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\scaffold\ms_scaffolds.fa"
    import largeFastaFile
    lsscf = list(largeFastaFile.open_fasta_to_list(filename))
    #find bases used in this file
    basesUsed=set()
    for ele in lsscf:
        basesUsed.update(list(ele.seq))
    print(basesUsed)
    #answer is ATCGN
    genomeseq = ''.join(str(ele.seq) for ele in lsscf)
    print('length of genome is %d'%len(genomeseq))
    #genome length is 419424057
    genomeseq.count('A')
    genomeseq.count('T')
    genomeseq.count('C')
    genomeseq.count('G')
    genomeseq.count('N')


def msGenomeTopHatCoverage20160722(ls_libraries):
    """
    count bases mapped by TopHat in each of the 67 libraries
    """
    filename = r'E:\store_for_D\Transcriptome\20160721GenomeCoverageByRNAseq\20160720sequencingDepth.txt'
    fo = open(filename,'r',1000000000)
    dcFiles = {}#split the file to 67 small files, so that I can read the whole file into memory
    for _num in range(67):
        dcFiles[_num] = open(filename+str(_num),"w",10000000)
    for ele in fo:
        elements = ele.split()
        for _num in range(67):
            dcFiles[_num].write(elements[_num+2] + "\n")
    for _num in range(67):
        dcFiles[_num].close()
    fo.close()#split the file to 67 small files, each line with bases number of each position of genome
    
    #20160723 count total aligned bases
    def countBases(filebases):
        """
        count bases, given a filebases which is a filename of the 67 small files.
        """
        ls_bases = open(filebases)
        count = 0
        for ele in ls_bases:
            count += int(ele[:-1])
        return count
    countBases(filename + '0')
    for num in range(1,17):
        print(num,countBases(filename+str(num)))
    for num in range(17,33):
        print(num,countBases(filename+str(num)))
    for num in range(33,49):
        print(num,countBases(filename+str(num)))
    for num in range(49,67):
        print(num,countBases(filename+str(num)))
    
    
        
    
    #get the first two column to a new file
    filename = r'E:\store_for_D\Transcriptome\20160721GenomeCoverageByRNAseq\20160720sequencingDepth.txt'
    fo = open(filename,'r',1000000000)
    fout = open(filename + 'baseID', 'w',100000000)
    for ele in fo:
        elements = ele.split()
        fout.write(elements[0] + '\t'+elements[1] + '\n')
    fout.close()
    fo.close()
    
    
    #generate a file store contig name and length in the order in the sequencingDepth file
    filename = r'E:\store_for_D\Transcriptome\20160721GenomeCoverageByRNAseq\20160720sequencingDepth.txtbaseID'
    fo = open(filename,'r',1000000000)
    fout = open(filename + 'contigs', 'w',100000000)
    lscontigs = []
    stcontig = set()
    for ele in fo:
        elements = ele.split()
        contig = elements[0]
        if contig not in stcontig:
            lscontigs.append([contig,1])
            stcontig.add(contig)
        else:
            lscontigs[-1][-1] += 1
    for ele in lscontigs:
        fout.write(ele[0] + '\t' + str(ele[1])+'\n')
    fout.close()
    fo.close()




def msGenomeTopHatExpressedRegion20160724(ls_libraries):
    """
    get the length of genome with coverage of 1
    """
    filename = r'E:\store_for_D\Transcriptome\20160721GenomeCoverageByRNAseq\20160720sequencingDepth.txt'
    fo = open(filename + '66','r')
    ls = [2**14,]
    count = 0
    for ele in fo:
        if ele[0] != '0':
            if int(ele[:-1]) > ls[0]:
                count += 1
    print(count)
    
    #below is a for a .py file we will used later
    import argparse
    parser = argparse.ArgumentParser(description = 'output a file with base numbers of given range')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename, content lower boundary\tbaseCount')
    parser.add_argument('lsLowerValues', type = int, nargs = '+', help = 'list of values as lower boundary, like [1], or [1,2,4,8]')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    lsLower = arg.lsLowerValues
    lsCount = [0 for _i in range(len(lsLower))]
    fo = open(inputfile,'r')
    fout = open(outfile,'w')
    for ele in fo:
        if ele[0] != '0':
            coverage = int(ele[:-1])
            for n in range(len(lsLower)):
                if coverage >= lsLower[n]:
                    lsCount[n] += 1
    for n in range(len(lsLower)):
        fout.write('%d\t%d\n'%(lsLower[n], lsCount[n]))
    fout.close()
    fo.close()
    
    #below is for generating scripts for HPC
    for ele in range(67):
        fo = open("20160724Ms_"+str(ele)+"_SeqDepth.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N SeqDepth"+str(ele)+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/\n")
        fo.write("module load python\n")
        fo.write("python sequenceDepth.py ./data/20160720sequencingDepth.txt" + str(ele)+' ./result/' + str(ele)+'.txt')
        fo.write(" 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288\n")
        fo.close()
    
    
    #below is a for a .py file we will used later. use BPKM
    import argparse
    parser = argparse.ArgumentParser(description = 'output a file with base numbers of given range')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename, content lower boundary\tbaseCount')
    parser.add_argument('totalbase', type = int, help = 'total number of bases')
    parser.add_argument('lsLowerValues', type = int, nargs = '+', help = 'list of values as lower boundary, like [1], or [1,2,4,8]')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    total = arg.totalbase
    lsLower = arg.lsLowerValues
    lsLower2 = [i * total / 1e9 for i in lsLower]
    lsCount = [0 for _i in range(len(lsLower))]
    fo = open(inputfile,'r')
    fout = open(outfile,'w')
    for ele in fo:
        if ele[0] != '0':
            coverage = int(ele[:-1])
            for n in range(len(lsLower)):
                if coverage >= lsLower2[n]:
                    lsCount[n] += 1
    for n in range(len(lsLower)):
        fout.write('%d\t%d\n'%(lsLower[n], lsCount[n]))
    fout.close()
    fo.close()
    
    #below is for generating scripts for HPC
    lsLen = '2491948975 2992678379 1626115160 3018987962 510159596 314106550 315667388 1829214301 2326374942 972029428 910644130 2932068874 2968756239 5008954446 3530755689 498122359 3058127635 2256130035 1879143968 345778722 315220651 320783326 316384573 380244745 2735729409 252385830 1185284026 338536327 212598355 348403523 217038294 185678701 424433863 3785387913 2937345468 3834197107 3666785430 2437241236 1235871474 1303091800 2123525841 2386039119 417411040 2255302344 374806306 2658085256 347263609 1692494168 3001278049 1433808690 2803397562 2211962973 919673347 1015475452 832273726 1028228849 937534144 738506392 726197360 680908873 618155276 790736561 643037274 686585213 631448788 401034318 401080868'.split()
    for ele in range(67):
        fo = open("20160724Ms_"+str(ele)+"_SeqDepth.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N SeqDepth"+str(ele)+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/\n")
        fo.write("module load python\n")
        fo.write("python sequenceDepthBPKM.py ./data/20160720sequencingDepth.txt" + str(ele)+' ./result/' + str(ele)+'.txt')
        fo.write(' '+ lsLen[ele])
        fo.write(" 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288\n")
        fo.close()
    
    
    #combine output of TopHat coverage
    folder = 'E:\\store_for_D\\Transcriptome\\20160721GenomeCoverageByRNAseq\\basesBPKM\\'
    lslsResults = []
    for num in range(67):
        lsResults = open(folder+str(num)+'.txt').readlines()
        lslsResults.append(lsResults)
    fout = open(folder+'20160725TopHatBaseCoverage.txt','w')
    for num in range(67):
        for num2 in range(len(lsResults)):
            fout.write(lslsResults[num][num2].split()[1] +'\t')
        fout.write('\n')
    fout.close()
    
    #count translated length based on TopHat all
    folder = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/data'
    dc_Depthfile = {}
    for num in range(67):
        dc_Depthfile[num] = open(folder + '20160720sequencingDepth.txt' + str(num),'r')
    count = 0
    for num in range(419424057):
        depth = 0
        for num2 in range(67):
            depth += int(dc_Depthfile[num2].readline()[:-1])
        if depth > 0:
            count += 1
    fout = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/20160705depth.txt','w')
    fout.write(str(count))
    fout.close()
    
    #20160726 convert TopHat base depth to BPKM
    ## python script
    import argparse
    parser = argparse.ArgumentParser(description = 'output a file with base numbers of given range')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename')
    parser.add_argument('totalbase', type = int, help = 'total number of bases')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    total = arg.totalbase
    fo = open(inputfile,'r')
    fout = open(outfile,'w')
    for ele in fo:
        if ele[0] == '0':
            fout.write(ele)
        else:
            coverage = int(ele[:-1])
            bpkm = int(round(coverage *1e9 / total))
            fout.write(str(bpkm) + '\n')
    fout.close()
    fo.close()
    
    
    
    lsLen = '2491948975 2992678379 1626115160 3018987962 510159596 314106550 315667388 1829214301 2326374942 972029428 910644130 2932068874 2968756239 5008954446 3530755689 498122359 3058127635 2256130035 1879143968 345778722 315220651 320783326 316384573 380244745 2735729409 252385830 1185284026 338536327 212598355 348403523 217038294 185678701 424433863 3785387913 2937345468 3834197107 3666785430 2437241236 1235871474 1303091800 2123525841 2386039119 417411040 2255302344 374806306 2658085256 347263609 1692494168 3001278049 1433808690 2803397562 2211962973 919673347 1015475452 832273726 1028228849 937534144 738506392 726197360 680908873 618155276 790736561 643037274 686585213 631448788 401034318 401080868'.split()
    file1 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/data/20160720sequencingDepth.txt'
    file2 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/dataBPKM/20160720sequencingDepth.txt'
    for ele in range(67):
        fo = open("20160726Ms_"+str(ele)+"_SeqDepth.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N SeqDepth"+str(ele)+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load python\n")
        fo.write("python /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/sequenceDepthCovert2BPKM.py ")
        fo.write(file1 + str(ele)+ " " +file2 + str(ele)+ " " + lsLen[ele] +'\n')
        fo.close()
    
    #20161020 convert TopHat base depth to BPKM. Remove mt and rRNA scaffolds and scaffolds short than 200kbp first,
    # then calculate bpkm
    ## python script
    import argparse
    parser = argparse.ArgumentParser(description = 'output a file with base numbers of given range')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    filecontig = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/20160720sequencingDepth.IDcontigs'
    contig_toRemove = ['AIXA01038380.1','AIXA01038379.1','AIXA01038378.1','AIXA01032915.1','AIXA01021581.1','AIXA01037114.1','AIXA01021582.1']
    contignames = []
    contiglens = []
    for ele in open(filecontig):
        ele1,ele2 = ele.split()
        contigname = ele1.split('|')[3]
        contiglen = int(ele2)
        contignames.append(contigname)
        contiglens.append(contiglen)
    
    fo = open(inputfile,'r',2000000000)
    fout = open(outfile,'w',1000000)
    fouttemp = open(outfile+'temp','w',1000000)
    total = 0
    for num in range(len(contignames)):
        if contignames[num] not in contig_toRemove and contiglens[num] > 200000:
            for n in range(contiglens[num]):
                depth = fo.readline()
                fouttemp.write(depth)
                depth = int(depth[:-1])
                total += depth
        else:
            for n in range(contiglens[num]):
                depth = fo.readline()
    fouttemp.close()
    fo.close()
    fo = open(outfile+'temp','r',2000000000)
    for ele in fo:
        if ele[0] == '0':
            fout.write(ele)
        else:
            coverage = int(ele[:-1])
            bpkm = int(round(coverage *1e9 / total))
            fout.write(str(bpkm) + '\n')
    fout.close()
    fo.close()
    
    
    #generateScripts
    file1 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/data/20160720sequencingDepth.txt'
    file2 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/dataBPKM/20160720sequencingDepth.txt'
    for ele in range(67):
        fo = open("20161021Ms_"+str(ele)+"_SeqDepth.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N SeqDepth"+str(ele)+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load python\n")
        fo.write("python /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/20161021convert2BPKMnoMtNo_rRNAlongerthan200kbp.py ")
        fo.write(file1 + str(ele)+ " " +file2 + str(ele) +'\n')
        fo.close()

    
    #20160726Counter of BPKM TopHat aligned bases
    import argparse
    parser = argparse.ArgumentParser(description = 'output Counter')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    fo = open(inputfile,'r')
    fout = open(outfile,'w')
    count = {}
    for ele in fo:
        bpkm = int(ele[:-1])
        if bpkm not in count:
            count[bpkm] = 0
        count[bpkm] += 1
    count = list(count.items())
    count.sort()
    for m, n in count:
        fout.write('%d\t%d\n'%(m,n))
    fout.close()
    fo.close()
    
    
    file1 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/dataBPKM/20160720sequencingDepth.txt'
    file2 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/result/20160720sequencingDepth.txt'
    for ele in range(67):
        fo = open("20160726Ms_"+str(ele)+"_SeqDepth.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q express\n")
        fo.write("#PBS -N SeqDepth"+str(ele)+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load python\n")
        fo.write("python /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/BPKMcounter.py ")
        fo.write(file1 + str(ele)+ " " +file2 + str(ele) +'\n')
        fo.close()
    
    # I got files with bpkm and count of bases with that bpkm(rounded to int). file looks like: 0\t381372537\n1\t10738781\n
    # I want to devide the bases to N groups, each group with same sum of bpkm, and count number of bases
    folder = 'E:\\store_for_D\\Transcriptome\\20160721GenomeCoverageByRNAseq\\basesBPKM\\BPKMconterWithoutMt\\'
    GROUPS = 10
    lsCount = []
    for eleFile in range(67):
        filename = folder + '20160720sequencingDepth.txt' + str(eleFile)
        fo = open(filename)
        bpkmV = []
        bpkmC = []
        for ele in fo:
            bpkmValue, bpkmCount = ele.split()
            bpkmV.append(int(bpkmValue))
            bpkmC.append(int(bpkmCount))
        fo.close()
        totalBPKM = 0
        for num in range(len(bpkmC)):
            totalBPKM += bpkmV[num] * bpkmC[num]
        targetBPKMs = [totalBPKM/GROUPS * (num+1) for num in range(GROUPS)]
        lsBases = [0,]
        for targetBPKM in targetBPKMs:
            bases = 0
            totalbpkm = 0
            for num in range(len(bpkmC)):
                totalbpkm += bpkmV[num] * bpkmC[num]
                bases += bpkmC[num]
                if totalbpkm >= targetBPKM:
                    bases = bases - (totalbpkm - targetBPKM) / bpkmV[num]
                    lsBases.append(bases)
                    break
        lsCount.append([lsBases[num+1] - lsBases[num] for num in range(GROUPS)])
    fout = open('20161017MsTopHatAlignmentDistributionDepthNoMt.txt', 'w')
    for num in range(GROUPS):
        for num2 in range(67):
            fout.write(str(lsCount[num2][num])+'\t')
        fout.write('\n')
    fout.close()
    
    ###top 400 to Top 1600 to 401, 6400 to 1601, ... average bpkm. 10 groups
    ls_boundary = [400* 2 ** num for num in range(0,21,1)]
    GENOMESIZE = 419424057
    ls_boundary[-1] = GENOMESIZE
    ls_boundarylen = [ls_boundary[0]]+[ ls_boundary[n+1] - ls_boundary[n] for n in range(len(ls_boundary)-1)]
    lsAveBPKM = []
    folder = 'E:\\store_for_D\\Transcriptome\\20160721GenomeCoverageByRNAseq\\basesBPKM\\BPKMconterWithoutMt_rRNA_shorterThan200k\\'
    for eleFile in range(67):
        filename = folder + '20160720sequencingDepth.txt' + str(eleFile)
        fo = open(filename)
        bpkmV = []
        bpkmC = []
        bpkmBoundary =[0 for i in range(len(ls_boundary))]
        for ele in fo:
            bpkmValue, bpkmCount = ele.split()
            bpkmV.append(int(bpkmValue))
            bpkmC.append(int(bpkmCount))
        fo.close()
        bpkmV.reverse()
        bpkmC.reverse()
        baseSum = 0
        bpkmSum = 0
        ls_boundaryNum = 0
        for num in range(len(bpkmV)):
            baseSum += bpkmC[num]
            bpkmSum += bpkmV[num] * bpkmC[num]
            if baseSum >= ls_boundary[ls_boundaryNum]:
                bpkmBoundary[ls_boundaryNum] += bpkmSum - (baseSum - ls_boundary[ls_boundaryNum]) * bpkmV[num]
                bpkmSum = (baseSum - ls_boundary[ls_boundaryNum]) * bpkmV[num]                
                ls_boundaryNum += 1
        lsAveBPKM.append([bpkmBoundary[_n]/ls_boundarylen[_n] for _n in range(len(bpkmBoundary))])
    fout = open('20161021MsTopHatAlignmentDistributionDepthTopBasesNoMt_rRNAover200.txt', 'w')
    for num in range(len(lsAveBPKM[0])):
        for num2 in range(67):
            fout.write(str(lsAveBPKM[num2][num])+'\t')
        fout.write('\n')
    fout.close()
    
            
            
        
    
    
    ##20160726 hot regions in the genome
    import argparse
    parser = argparse.ArgumentParser(description = 'output Counter')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    fo = open(inputfile,'r')
    fout = open(outfile,'w')
    sep = 100000
    count = 0
    num = 0
    for ele in fo:
        num += 1
        count += int(ele[:-1])
        if num == sep:
            fout.write("%f\n"%(float(count)/num))
            count = 0
            num = 0            
    fout.write("%.1f\n"%(float(count)/num))
    fout.close()
    fo.close()
    
    
    file1 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/dataBPKM/20160720sequencingDepth.txt'
    file2 = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/result/20160720sequencingDepth.txt'
    for ele in range(67):
        fo = open("20160726Ms_"+str(ele)+"_SeqDepth.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N SeqDepth"+str(ele)+"\n")
        fo.write("#PBS -l nodes=1:ppn=1\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("module load python\n")
        fo.write("python /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/BPKMcounter.py ")
        fo.write(file1 + str(ele)+ " " +file2 + str(ele) +'\n')
        fo.close()
    
    ##20160729 hot regions in the genome data dealwith
    folder = 'E:\\store_for_D\\Transcriptome\\20160721GenomeCoverageByRNAseq\\basesBPKM\\BPKMgenomeHotRegionsNoMt\\'
    filelist =[folder + '20160720sequencingDepth.txt'+str(_n) for _n in range(67)]
    bpkmRegionAll =[open(_file).readlines() for _file in filelist]
    fout = open('20161017GenomeHotRegionsNoMt.txt','w')
    for ele in bpkmRegionAll:
        for ele1 in ele:
            fout.write(ele1[:-1]+'\t')
        fout.write('\n')
    fout.close()
    
    #Genome HotRegions 20161022 only use scaffolds longer than 1m. sep = 100kd
    import argparse
    parser = argparse.ArgumentParser(description = 'output Counter')
    parser.add_argument('inputfilename', type = str, help = 'input filename')
    parser.add_argument('outputfilename', type = str, help = 'output filename')
    arg = parser.parse_args()
    inputfile = arg.inputfilename
    outfile = arg.outputfilename
    filecontig = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/20160720sequencingDepth.IDcontigs'
    contig_toRemove = ['AIXA01038380.1','AIXA01038379.1','AIXA01038378.1','AIXA01032915.1','AIXA01021581.1','AIXA01037114.1','AIXA01021582.1']
    contignames = []
    contiglens = []
    for ele in open(filecontig):
        ele1,ele2 = ele.split()
        contigname = ele1.split('|')[3]
        contiglen = int(ele2)
        contignames.append(contigname)
        contiglens.append(contiglen)
    
    
    
    fo = open(inputfile,'r',1000000000)
    fout = open(outfile,'w',100000000)
    sep = 100000
    count = 0
    num = 0
    skip = 0
    skip_short1m = 0
    keep1m = 0
    for nn in range(len(contignames)):
        if contignames[nn] in contig_toRemove or contiglens[nn] <= 200000:
            for n in range(contiglens[nn]):
                skip += 1
        elif contiglens[nn] <= 1000000:
            for n in range(contiglens[nn]):
                fo.readline()
                skip_short1m += 1
        else:
            for n in range(contiglens[nn]):
                ele = fo.readline()
                count += int(ele[:-1])
                num += 1
                keep1m += 1
                if num == sep:
                    fout.write("%f\n"%(float(count)/num))
                    count = 0
                    num = 0
    fout.write("%.1f\n"%(float(count)/num))
    fout.close()
    fo.close()
    
    print('short than 200k, 200k to 1000k, over 1000k:', skip, skip_short1m, keep1m)    
    
    ##20160729 hot regions in the genome data dealwith
    folder = 'E:\\store_for_D\\Transcriptome\\20160721GenomeCoverageByRNAseq\\basesBPKM\\GenomeHotRegionsLongerThan1M\\'
    filelist =[folder + '20160720sequencingDepth.txt'+str(_n) for _n in range(67)]
    bpkmRegionAll =[open(_file).readlines() for _file in filelist]
    fout = open('20161022GenomeHotRegionsLongerThan1M.txt','w')
    for ele in bpkmRegionAll:
        for ele1 in ele:
            fout.write(ele1[:-1]+'\t')
        fout.write('\n')
    fout.close()


def msCufflinksV4GTFtranslatedLength20160724():
    """
    from the gff file, get the length of transcriptable part of the genome
    """
    filename = r'E:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\GTF_GFF3\20160713ManducaSexta_merged.gtf'
    fo = open(filename,'r')
    lslsgff = []
    for ele in fo:
        eles = ele.split()
        lslsgff.append([eles[0], int(eles[3]), int(eles[4])])
    fo.close()
    dcContigGff = {}
    for ele in lslsgff:
        if ele[0] not in dcContigGff:
            dcContigGff[ele[0]] = []
        dcContigGff[ele[0]].append(ele[1:])
    dcContigTranslength = {}
    for ele in dcContigGff:
        dcContigTranslength[ele] = set()
        for eleStart, eleEnd in dcContigGff[ele]:
            dcContigTranslength[ele].update(list(range(eleStart,eleEnd+1)))
        dcContigTranslength[ele] = len(dcContigTranslength[ele])
    print(sum(dcContigTranslength.values()))

def msCufflinksV4GTFCDSLength20161023():
    """
    from the gff file, get the length of transcriptable part of the genome
    """
    filename = r"D:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\transcriptsAndTranslation\20160713ManducaCufflinksTranscripts.fa.transdecoder.genome.gff3"
    fo = open(filename,'r')
    lslsgff = []
    for ele in fo:
        if ele == '\n':
            continue
        eles = ele.split()
        if eles[2] == "CDS":
            lslsgff.append([eles[0], int(eles[3]), int(eles[4])])
    fo.close()
    dcContigGff = {}
    for ele in lslsgff:
        if ele[0] not in dcContigGff:
            dcContigGff[ele[0]] = []
        dcContigGff[ele[0]].append(ele[1:])
    dcContigTranslength = {}
    for ele in dcContigGff:
        dcContigTranslength[ele] = set()
        for eleStart, eleEnd in dcContigGff[ele]:
            dcContigTranslength[ele].update(list(range(eleStart,eleEnd+1)))
        dcContigTranslength[ele] = len(dcContigTranslength[ele])
    print(sum(dcContigTranslength.values()))
    
    #20170709 OGS2
    filename = r"D:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2_2016\ms_ogs.gff"
    fo = open(filename,'r',encoding='utf-8')
    lslsgff = []
    for ele in fo:
        if ele == '\n' or ele[0] == '#':
            continue
        eles = ele.split()
        if eles[2] == "CDS":
            lslsgff.append([eles[0], int(eles[3]), int(eles[4])])
    fo.close()
    dcContigGff = {}
    for ele in lslsgff:
        if ele[0] not in dcContigGff:
            dcContigGff[ele[0]] = []
        dcContigGff[ele[0]].append(ele[1:])
    dcContigTranslength = {}
    for ele in dcContigGff:
        dcContigTranslength[ele] = set()
        for eleStart, eleEnd in dcContigGff[ele]:
            dcContigTranslength[ele].update(list(range(eleStart,eleEnd+1)))
        dcContigTranslength[ele] = len(dcContigTranslength[ele])
    print(sum(dcContigTranslength.values()))

def DmCufflinksV4GTFCDSLength20161023():
    """
    from the gff file, get the length of transcriptable part of the genome
    """
    filename = r"D:\Insects\Drosophila_melanogaster\dmel-all-r6.16.gtf"
    fo = open(filename,'r')
    lslsgff = []
    targetwords = 'exon'.split()
    for ele in fo:
        if ele == '\n':
            continue
        eles = ele.split()
        if eles[2] in targetwords:
            lslsgff.append([eles[0], int(eles[3]), int(eles[4])])
    fo.close()
    dcContigGff = {}
    for ele in lslsgff:
        if ele[0] not in dcContigGff:
            dcContigGff[ele[0]] = []
        dcContigGff[ele[0]].append(ele[1:])
    dcContigTranslength = {}
    for ele in dcContigGff:
        dcContigTranslength[ele] = set()
        for eleStart, eleEnd in dcContigGff[ele]:
            dcContigTranslength[ele].update(list(range(eleStart,eleEnd+1)))
        dcContigTranslength[ele] = len(dcContigTranslength[ele])
    print(sum(dcContigTranslength.values()))

def msSTARAlignment20160801(ls_libraries):
    '''
    align the trimmed reads to genome with STAR, and with the gtf file of cufflinks V4. Two alignment
    '''
    for ele in ls_libraries:
        fo = open("20160801Ms_"+ele['ShortName']+"_STAR.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N STAR"+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=1:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/RNA-trim/\n")
        fo.write("module load jdk\n")
        fo.write("module load samtools\n")
        fo.write("module load cufflinks\n")
        fo.write("module load bowtie2\n")
        fo.write("mkdir ../STAR/"+ele['ShortName']+'\n')
        fo.write("/panfs/panfs.cluster/home/ks2073/p/2016/STAR/bin/Linux_x86_64_static/STAR  --runThreadN 12  --runMode alignReads --genomeDir ")
        fo.write("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/GenomeMs/STAR/ ")
        fo.write("--outFileNamePrefix ../STAR/"+ele['ShortName']+'/ --outFilterMultimapNmax 10000 ')
        fo.write('--outReadsUnmapped Fastx --outSAMtype None --outSAMmode None --twopassMode Basic --readFilesIn ')
        if ele['File2'] == '':
            file_reads = ele['File1']
            fo.write(file_reads)
        else:
            file_reads_p1 = ele['File1']+"paired.fq"
            file_reads_p2 = ele['File2']+"paired.fq"
            file_reads_s1 = ele['File1']+"unpaired.fq"
            file_reads_s2 = ele['File2']+"unpaired.fq"
            fo.write(file_reads_p1 +"," + file_reads_p2+","+file_reads_s1+","+file_reads_s2)
        
        
        fo.write("\n")
        fo.close()
        
    ##deal with the output
    folder = 'E:\\store_for_D\\Transcriptome\\unMappedReads20160801\\StarFinalResult\\'
    for ele in ls_libraries:
        fo = open(folder +ele['ShortName']+'.txt')
        templs = fo.readlines()
        fo.close()
        uniqueReads = int(templs[8].split()[-1])
        multipleReads = int(templs[23].split()[-1])
        print(ele['ShortName'], uniqueReads+multipleReads)


def msSTARunmappedReadsBlastn20160807(ls_libraries):
    """
    blast the unmapped reads of each libraries to ncbi nt database
    """
    for ele in ls_libraries:
        fo = open("20160807Ms_"+ele['ShortName']+"_blastnUnmap.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N STAR"+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/\n")
        fo.write("module load blast+\n")
        fo.write("module load python\n")
        fastqname = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta//STAR/'+ ele['ShortName'] + '/Unmapped.out.mate1'
        fastaname = fastqname+'.fa'
        fo.write('''python -c "from Bio import SeqIO;SeqIO.convert('%s','fastq','%s','fasta')"  \n'''%(fastqname, fastaname))
        fo.write("blastn -db /opt/ncbi/databases/nt -outfmt 6 -num_threads 12 -query ./STAR/")
        fo.write(ele['ShortName'] + '/Unmapped.out.mate1.fa -out ./blastUnmappedReads/'+ele['ShortName'] + '.txt ')
        fo.write(" -max_target_seqs 1")
        fo.write("\n")
        fo.close()
    
def msSTARunmappedReadsBlastnResultKeepFirstLine(ls_libraries):
    """
    for result in previous step, each query have maximum 1 target sequence, and more than one match against that target
    only keep the first match for each query and subject pair
    """
#    folder = ''
    for ele in ls_libraries:
        None

def tRNA_dataExtract20160810():
    """
    "E:\store_for_D\Transcriptome\20160722tRNAgenesManducaSexta\tRNAGenesNonePseudo20180810.txt"
    return the location of tRNA in the genome (file like 20160720sequencingDepth.txt0), begin from 1
    """
    f_contigs = r"E:\store_for_D\Transcriptome\20160721GenomeCoverageByRNAseq\20160720sequencingDepth.IDcontigs"
    ls_contigs = []#store contig id, length of contig
    for ele in open(f_contigs):
        ls_contigs.append((ele.split()[0], int(ele.split()[1])))
    
    f_tRNA = r"E:\store_for_D\Transcriptome\20160722tRNAgenesManducaSexta\tRNAGenesNonePseudo20180810.txt"
    ls_tRNA = []#store contig id, begin, end. if in the reverse strand, modify so that begin is smaller than end
    fopen_tRNA = open(f_tRNA)
    fopen_tRNA.readline()
    for ele in fopen_tRNA:
        eles = ele.split()
        ele_start = int(eles[2])
        ele_end = int(eles[3])
        if ele_start < ele_end:
            ls_tRNA.append([eles[0], ele_start, ele_end])
        else:
            ls_tRNA.append([eles[0], ele_end, ele_start])
    
    #ls_tRNA is a list of elements like '['gi|388482364|gb|JH668279.1|', 264831, 264902]'
    # add three field to each element, length, the begin and ends
    for ele in ls_tRNA:
        ele_length = ele[2] - ele[1] + 1
        ele.append(ele_length)
        num = 0
        for contig, contiglen in ls_contigs:
            if ele[0] != contig:
                num += contiglen
            else:
                break
        ele.append(num + ele[1])
        ele.append(num + ele[2])
    #Store ls_tRNA
    fout = open('tRNAgenomeLocation.txt','w')
    for ele in ls_tRNA:
        fout.write('%s\t%d\t%d\t%d\t%d\t%d\n'%(ele[0],ele[1],ele[2],ele[3],ele[4],ele[5],))
    fout.close()
    
    files = r'E:\store_for_D\Transcriptome\20160721GenomeCoverageByRNAseq\20160720sequencingDepth.txt'
    for num in range(67):
        filename = files + str(num)
        ls_genomeCov =map(int, open(filename).read().split())
        for ele in ls_tRNA:
            ele_bpkm = sum(ls_genomeCov[ele[4]-1:ele[5]]) /ele[3]
            ele.append(ele_bpkm)
    #save ls_tRNA
    
    fout = open('tRNAgenomeBPKMlibraries.txt','w')
    for ele in ls_tRNA:
        for ele2 in ele:
            fout.write(str(ele2) + '\t')
        fout.write('\n')
    fout.close()
    
    #conclusion, most tRNA genes were not sequenced in different libraries.


def unmappedBlastResultKeepOne20160811():
    """
    for all unmapped reads in each library, the blastn program only keeps one target. 
    Still, the read might map to multiple sites of one target. Only keep the first hit
    """
    folder1 = 'E:\\store_for_D\\Transcriptome\\unMappedReads20160801\\blastnResult20160811\\'
    folder2 = 'E:\\store_for_D\\Transcriptome\\unMappedReads20160801\\blastnResultKeepOne20160811\\'
    import os
    filenames = os.listdir(folder1)
#    for filename in filenames:
    def keepOne(filenameIn, filenameOut):
        lsLines = open(filenameIn).readlines()
        fout = open(filenameOut,'w')
        fout.write(lsLines[0])
        for num in range(len(lsLines)-1):
            line1 = lsLines[num]
            line2 = lsLines[num+1]
            if line1.split()[0] != line2.split()[0]:
                fout.write(line2)
        fout.close()
    
    for ele in filenames:
        keepOne(folder1+ele, folder2+ele)
    
    folder3 = 'E:\\store_for_D\\Transcriptome\\unMappedReads20160801\\blastnResultKeepOneEvalueFilter20160811\\'
    def filterEvalue(filenameIn, filenameOut, eMin):
        lsLines = open(filenameIn).readlines()
        fout = open(filenameOut,'w')
        num = 0
        for ele in lsLines:
            if float(ele.split()[-2]) <= eMin:
                fout.write(ele)
                num += 1
        fout.close()
        print(filenameIn.split('\\')[-1], num)
    
    for ele in filenames:
        filterEvalue(folder2+ele, folder3+ele,1e-6)

def unmappedReadsBlastResultDealWith20160811(ls_libraries):
    """
    deal with the result generated in previous step
    """
    folder3 = 'E:\\store_for_D\\Transcriptome\\unMappedReads20160801\\blastnResultKeepOneEvalueFilter20160811\\'
    from collections import Counter
    dc_giCounter = {}
    dc_giCounter ['total'] = Counter()
    def countgis(filename):
        gis = [line.split()[1].split('|')[-2] for line in open(filename)]
        return Counter(gis)
    for ele in ls_libraries:
        filename = folder3 + ele["ShortName"] + '.txt'
        key = ele["ShortName"]
        dc_giCounter[key] = countgis(filename)
        dc_giCounter['total'].update(dc_giCounter[key])
    print(len(dc_giCounter['total']))
    
    fout = open('20160811MsUnmappedReadsAccList.txt','w')
    giSort = dc_giCounter['total'].most_common()
    for giNo in giSort:
        giNo = giNo[0]
        fout.write(giNo)
        fout.write('\t%d'%(dc_giCounter['total'][giNo]))
        for ele in ls_libraries:
            fout.write('\t%d'%(dc_giCounter[ele['ShortName']][giNo]))
        fout.write('\n')
    fout.close()
    
    import seqExtract
    ls_acc = [ele[0] for ele in giSort]
    ls_acc[:3]
    seqExtract.downloadFastaFromNCBI(ls_acc[94800:])


def unmappedReadsBlastResultDealWith20160816(ls_libraries):
    """
    further deal with the result of previous step
    """
    folder = 'E:\\store_for_D\\Transcriptome\\unMappedReads20160801\\'
    f_ncbiAcc2name = folder + '20160816NCBI_Acc2Names.txt'
    f_ncbiAcc2nameClean = folder + '20160816NCBI_Acc2NamesClean.txt'
    
    #clean the file f_ncbiAcc2name
    fout = open(f_ncbiAcc2nameClean,'w')
    for ele in open(f_ncbiAcc2name):
        ntID, ntName = ele.split(' ',1)
        ntID = ntID.split('|')[3]
        fout.write(ntID +'\t'+ntName)
    fout.close()
    
    fo = open(f_ncbiAcc2nameClean)
    ls = []
    import re
    for ele in fo:
        ls.append(ele.split('\t'))
    
    for ele in ls:
        ele.append(0)
    for ele in ls:
        if ele[2] == 0:
            if re.search('ribosomal RNA|rRNA', ele[1], re.IGNORECASE):
                ele[2] = 'rRNA'
            elif re.search("phage", ele[1], re.IGNORECASE):
                ele[2] = 'phage'
            elif re.search("m.sexta|manduca| sexta", ele[1], re.IGNORECASE):
                ele[2] = 'manduca'
            elif re.search("Escherichia coli|e.coli|e. coli", ele[1], re.IGNORECASE):
                ele[2] = 'E.coli'
            elif re.search("Oryza ", ele[1], re.IGNORECASE):
                ele[2] = 'Oryza'
                
                
def boxplotOfReadsInfo20161010():
    P_readsnum = [32871436,40185672,22070754,44383880,29293474,32747108,15172318,17941016,45154416,45642828,73341684,56128290,42871586,32823050,28926682,37647944,17687930,58038170,45101558,64746514,55197306,37180536,19956562,21923010,33916950,38503332,37585098,40415404,24431546,44396480,20673204,40572690,34558776]
    S_readsnum = [11296498,6945200,6977117,11533022,7708329,6995402,7356017,7086348,8435005,6092462,7683464,4788970,8027650,4901125,4217238,10124358,9520822,8262108,7618332]
    H_readsnum = [21935947,23890715,20008392,24333223,22284518,17548013,17403333,16588623]
    An_readsnum = [8808794,11144629,9839501,9442389,8837391,9010179,9265614]
    
    P_Trim = [0.93125171,0.916323584,0.912327191,0.891705277,0.914489077,0.929461557,0.937802121,0.923821204,0.861757729,0.869672734,0.859531368,0.84645064,0.905439374,0.928598683,0.930866976,0.912000108,0.900118951,0.869616461,0.867810154,0.833773228,0.864633176,0.891224591,0.91289657,0.845746866,0.856192258,0.864064778,0.863320591,0.870116528,0.862261643,0.856988324,0.862764475,0.863119428,0.833191314]
    S_Trim = [0.939336509,0.94174221,0.94117685,0.932808764,0.941873654,0.942117837,0.932368427,0.948340527,0.944579997,0.960163067,0.994791542,0.991286644,0.986276058,0.997124334,0.994081434,0.957752778,0.933975239,0.938359799,0.94353554]
    H_Trim = [0.906167488,0.905654603,0.90262091,0.905925204,0.904307107,0.905180832,0.903169812,0.900980389]
    An_Trim = [0.890928202,0.894194594,0.862283667,0.941384961,0.939072742,0.61521164,0.617502197]
    
    P_star = [0.936200257,0.923104692,0.934278581,0.896410547,0.825215889,0.897920399,0.820799638,0.688874957,0.890580471,0.884550317,0.933447276,0.911330113,0.942062783,0.90011805,0.823580463,0.971830144,0.93621408,0.913659399,0.901809268,0.853023373,0.912549878,0.921938322,0.84948917,0.877035056,0.886647637,0.866852495,0.829978708,0.900797589,0.960678968,0.944081684,0.955898614,0.94121921,0.917228553]
    P_tophat = [0.86,0.859,0.865,0.829,0.74,0.818,0.728,0.597,0.807,0.799,0.853,0.804,0.86,0.814,0.76,0.879,0.835,0.8,0.797,0.763,0.817,0.833,0.75,0.784,0.811,0.79,0.759,0.817,0.856,0.843,0.856,0.854,0.834]
    S_star = [0.976955698,0.978654213,0.982505218,0.970800896,0.977193692,0.978910075,0.976198065,0.964853798,0.976682003,0.903448126,0.930643054,0.94365065,0.940532081,0.935003277,0.950878019,0.922558059,0.976329624,0.9836056,0.985079367]
    S_tophat = [0.943,0.942,0.943,0.908,0.934,0.938,0.917,0.923,0.936,0.846,0.869,0.878,0.863,0.871,0.868,0.858,0.921,0.948,0.947]
    H_star = [0.952862719,0.961814712,0.945891895,0.958394862,0.954828264,0.954209376,0.950869838,0.939167367]
    H_tophat = [0.907,0.921,0.904,0.915,0.912,0.912,0.906,0.894]
    An_star = [0.968819711,0.969399126,0.923438091,0.962671856,0.964285736,0.919511536,0.915904765]
    An_tophat = [0.866,0.867,0.825,0.833,0.821,0.781,0.757]
    
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams["figure.figsize"] = [4.0,8]
    pp = PdfPages('20161010boxPlotReadsManduca.pdf')
    bp = plt.boxplot([P_readsnum, S_readsnum, H_readsnum, An_readsnum],showmeans=1)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle = '-')
    plt.setp(bp['fliers'], color='red', marker='x')
    pp.savefig()
    plt.close()
    bp = plt.boxplot([P_Trim, S_Trim, H_Trim, An_Trim],showmeans=1)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle = '-')
    plt.setp(bp['fliers'], color='red', marker='x')
    pp.savefig()
    plt.close()
    plt.rcParams["figure.figsize"] = [8.0,8]
    bp = plt.boxplot([P_star,P_tophat, S_star,S_tophat, H_star, H_tophat, An_star,An_tophat],showmeans=1)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle = '-')
    plt.setp(bp['fliers'], color='red', marker='x')
    pp.savefig()
    plt.close()
    pp.close()
    
    P_mt = [0.02596841,0.016523229,0.017038096,0.007071722,0.045562044,0.033347557,0.030308836,0.010167476,0.085138735,0.10440908,0.038490179,0.010162706,0.010901262,0.007712684,0.005099659,0.0364107,0.023046312,0.028443672,0.138247718,0.083221463,0.041284423,0.011906248,0.052384216,0.027736562,0.042205324,0.044978242,0.092629041,0.141031644,0.034193828,0.044870226,0.03455992,0.111559523,0.052070047]
    S_mt = [0.160821777,0.150624702,0.065592876,0.010050767,0.04482501,0.090503396,0.084543448,0.203299243,0.188413389,0.184504034,0.111381217,0.176414051,0.206069168,0.190414056,0.153245804,0.061359455,0.059683097,0.154488786,0.169932484]
    H_mt = [0.029469164,0.021694182,0.019173312,0.029548486,0.019484688,0.020509627,0.017782133,0.031692828]
    An_mt = [0.0169601,0.035565305,0.189332768,0.021506215,0.022822548,0.023035482,0.022653059]
    P_cov = [0.138083076,0.153231037,0.093275096,0.10380462,0.074045259,0.132970506,0.081931309,0.148585573,0.159732118,0.131861952,0.078357475,0.088804968,0.158674074,0.11483249,0.124286095,0.142115782,0.093964026,0.172373439,0.222376858,0.203026661,0.189051261,0.145377109,0.156246216,0.127432674,0.146318908,0.126700215,0.078705483,0.108974722,0.167200576,0.232157484,0.174006273,0.196630877,0.14743465]
    S_cov = [0.025636961,0.020030825,0.01741485,0.076806593,0.034182729,0.035644288,0.054386363,0.019857716,0.021795776,0.048466264,0.041473601,0.035273589,0.043463756,0.0364829,0.033602848,0.053624368,0.046522532,0.015378507,0.024444113]
    H_cov = [0.098756169,0.107365344,0.087468681,0.101487634,0.091010376,0.087647247,0.084007122,0.077337297]
    An_cov = [0.088824402,0.092301766,0.079500867,0.101909581,0.100996772,0.085868031,0.082483583]
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams["figure.figsize"] = [4.0,8]
    pp = PdfPages('20161013boxPlotReadsManducaMTAndCoverageRatio.pdf')
    bp = plt.boxplot([P_mt, S_mt, H_mt, An_mt],showmeans=1)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle = '-')
    plt.setp(bp['fliers'], color='red', marker='x')
    pp.savefig()
    plt.close()
    bp = plt.boxplot([P_cov, S_cov, H_cov, An_cov],showmeans=1)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle = '-')
    plt.setp(bp['fliers'], color='red', marker='x')
    pp.savefig()
    plt.close()
    pp.close()
    
    P_rRNA = [0.893383046,0.846908302,0.925963935,0.950808401,0.927961331,0.891872813,0.9243705,0.929871101,0.91987379,0.69906462,0.903724864,0.906824503,0.889936045,0.903469797,0.299953422,0.549386286,0.840176354,0.884881995,0.881607329,0.939670895,0.918077799,0.917783031,0.926812455,0.924226838,0.932709873,0.938860402,0.661124069,0.894738085,0.79680663,0.873867424,0.820065212,0.801439673,0.899296529]
    S_rRNA = [0.541166665,0.401113372,0.24599823,0.146989978,0.548509987,0.458518726,0.22568331,0.609528755,0.452556784,0.027471965,0.03621188,0.032706334,0.025640167,0.013914126,0.01905529,0.031314414,0.33891093,0.110099565,0.200949928]
    H_rRNA = [0.286453188,0.320968978,0.298878744,0.182353499,0.322252289,0.331527347,0.272159379,0.345958605]
    An_rRNA = [0.481794018,0.358983442,0.262196348,0.457045279,0.326767103,0.098115074,0.104662807]
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams["figure.figsize"] = [4.0,8]
    pp = PdfPages('20161019boxPlotReadsRatioTorRNA.pdf')
    bp = plt.boxplot([P_rRNA, S_rRNA, H_rRNA, An_rRNA],showmeans=1)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle = '-')
    plt.setp(bp['fliers'], color='red', marker='x')
    pp.savefig()
    plt.close()
    pp.close()
    
    def plotFigureBoxPlot(filename='list.txt', outname = 'box-plot.pdf',ordered = ['P','S','H','An']):
        '''
        filename, file contains data, looks like
        P 0.1
        S 0.2
        first part group name, second part value.
        plot a box plot
        '''
        templs = open(filename).readlines()
        dc={}
        for ele in templs:
            groupname, value = ele.split()
            value = float(value)
            if groupname not in dc:
                dc[groupname] = []
            dc[groupname].append(value)
        
        lsplot = [dc[ele] for ele in ordered]
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        plt.rcParams["figure.figsize"] = [4.0,8]
        pp = PdfPages(outname)
        bp = plt.boxplot(lsplot,showmeans=1)
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['whiskers'], color='black', linestyle = '-')
        plt.setp(bp['fliers'], color='red', marker='x')
        pp.savefig()
        plt.close()
        pp.close()
    plotFigureBoxPlot()
    
    
def mappingscore():
    """
    mapping score calculation. Defined in the paper.
    Comparison of D. melanogaster and C. elegans developmental stages, tissues, and cells by modENCODE RNA-seq data
    """
    file_tsv = r"E:\Lab\works\2013ManducaUtilities\Transcriptome\RSEM\OGS2_20140430\ForTranscriptomePaper20161011\20161012TissueSpecific.txt"
    templs = open(file_tsv).readlines()
    ls_z =[[] for i in range(67)]
    for templine in templs[1:]:
        tempzs = templine.split()
        tempzs = [float(ele) for ele in tempzs[1:]]
        for i in range(67):
            ls_z[i].append(tempzs[i])
    
    n = 15543
    comparisonNum = 67*67
    dc_pairs ={}
    for _i in range(67):
        for _j in range(_i,67):
            a = sum(1 if ele >0 else 0 for ele in ls_z[_i])
            b = sum(1 if ele >0 else 0 for ele in ls_z[_j])
            ab = sum(1 if ls_z[_i][_n] > 0 and ls_z[_j][_n] > 0 else 0 for _n in range(len(ls_z[_i])))
            dc_pairs[(_i,_j)] = [a,b,ab]
            
    def mappingScore_geneComparison(n,a,b, ab,comparisonNum):
        """
        for hypergenometic testing in within species stage/tissue/cell comparison.
        check paper Comparison of D. melanogaster and C. elegans developmental stages, tissues, and cells by modENCODE RNA-seq data
        method part
        """
        from decimal import Decimal
        import math
        n = Decimal(n)
        a = Decimal(a)
        b = Decimal(b)
        ab = Decimal(ab)
        comparisonNum = Decimal(comparisonNum)
    
        def fl(value):
            value = int(value)
            return Decimal(sum(math.log(i) for i in range(1,value+1)))
        def pvalue3(i):
            return Decimal(math.e) ** (fl(a) - fl(a-i) + fl(b) - fl(b-i) -fl(i) + fl(n-a) - fl(n) + fl(n-b) - fl(n+i-a-b))
        
        p = Decimal()
        for i in range(int(ab), int(min(a,b))+1):
            p += pvalue3(i)
        return - float(p.log10() + comparisonNum.log10())
    
    templs = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/paired_info.txt').readlines()
    ls_paras = []
    for templine in templs:
        eles = templine.split()
        eles = [int(ele) for ele in eles]
        ls_paras.append(eles)
    
    def m(paras):
        l1, l2, n,a,b,ab,comparisonNum = paras
        score = mappingScore_geneComparison(n,a,b, ab,comparisonNum)
        return l1, l2, score
    
    from multiprocessing import Pool
    print(ls_paras[:3])
    if __name__ == '__main__':
        p = Pool(12)
        c = p.map(m, ls_paras)
        fout = open('/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/testResult.txt','w')
        for ele in c:
            fout.write(str(ele) +'\n')
        fout.close()
    
    ls_paras =[]
    for ele in dc_pairs:
        a,b,ab = dc_pairs[ele]
        ls_paras.append((n,a,b,ab,comparisonNum))

    import time
    time1 = time.time()
    for ele in dc_pairs:
        a,b,ab = dc_pairs[ele]
        mappingScore_geneComparison(n,a,b,ab,comparisonNum)
        time2 = time.time()
        print(time2-time1)
        break
    
    fout = open('paired_info.txt','w')#store dc_pairs
    for ele in dc_pairs:
        _lib1 = ele[0]
        _lib2 = ele[1]
        a,b,ab = dc_pairs[ele]
        fout.write('%d\t%d\t%d\t%d\t%d\t%d\t%d\n'%(_lib1,_lib2,n,a,b,ab,comparisonNum))
    fout.close()
    
    fout = open('paired_info_organized.txt','w') #store dc_pairs in table
    _lsinfo = [[0 for _i in range(68)] for _j in range(68)]
    for _i in range(67):
        for _j in range(67):
            if (_i,_j) in dc_pairs:
                _values = dc_pairs[(_i,_j)]
            else:
                _values = dc_pairs[(_j,_i)]
                _values = [_values[1],_values[0],_values[2]]
            _lsinfo[_i+1][0] = _values[0]
            _lsinfo[0][_j+1] = _values[1]
            _lsinfo[_i+1][_j+1] = _values[2]
    for _ele in _lsinfo:
        fout.write('\t'.join(str(_e2) for _e2 in _ele))
        fout.write('\n')
    fout.close()
    
    fout = open(r"E:\Lab\works\2013ManducaUtilities\Transcriptome\RSEM\OGS2_20140430\ForTranscriptomePaper20161011\20161012ZscoreResultForFigure.txt",'w')
    _lsScore = open(r"E:\Lab\works\2013ManducaUtilities\Transcriptome\RSEM\OGS2_20140430\ForTranscriptomePaper20161011\20161012TissueSpecificResult.txt").readlines()
    dc_score = {}
    for ele in _lsScore:
        _i, _j, _S = ele.split()
        _i = int(_i)
        _j = int(_j)
        _S = float(_S)
        dc_score[(_i,_j)] = _S
    _lsinfo = [[0 for _i in range(68)] for _j in range(68)]
    for _i in range(67):
        for _j in range(67):
            if (_i,_j) in dc_pairs:
                _values = dc_score[(_i,_j)]
            else:
                _values = dc_score[(_j,_i)]
            _lsinfo[_i+1][0] = _i+1
            _lsinfo[0][_j+1] = _j+1
            _lsinfo[_i+1][_j+1] = _values
    for _ele in _lsinfo:
        fout.write('\t'.join(str(_e2) for _e2 in _ele))
        fout.write('\n')
    fout.close()
        
    

def noncodingGenesCufflinks20161019():
    fileNt = r"E:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\transcriptsAndTranslation\20160713ManducaCufflinksTranscripts.fa"
    filepep = r"E:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\transcriptsAndTranslation\20160713ManducaCufflinksTranscripts.fa.transdecoder.pep"
    import largeFastaFile
    lsNt = largeFastaFile.open_fasta_to_list(fileNt)
    lspep = largeFastaFile.open_fasta_to_list(filepep)
    lsCodingT = []
    for ele in lspep:
        codingName = ele.id.split(':')[2]
        lsCodingT.append(codingName)
    lsCodingT = list(set(lsCodingT))
    print(len(lsCodingT))
    lsNoncoding = []
    for ele in lsNt:
        if ele.id not in lsCodingT:
            lsNoncoding.append(ele)
    lsNoncoding = set(lsNoncoding)
    print(len(lsNoncoding))
    
    largeFastaFile.saveFastaListToFile(lsNoncoding,fileNt+'.noncodingT')
    lsNoncodingGenes = []
    for ele in lsNoncoding:
        geneID = ele.description.split()[1]
        if geneID not in lsNoncodingGenes:
            lsNoncodingGenes.append(geneID)
    print(len(lsNoncodingGenes))
    
    lsCodingGenes = set()
    lsGenes = set()
#    lsNoncodingT = [ele.id for ele in lsNoncoding]
    
    for ele in lsNt:
        geneID = ele.description.split()[1]
        if ele.id in lsCodingT:
            lsCodingGenes.add(geneID)
        lsGenes.add(geneID)
    print(len(lsCodingGenes), len(lsGenes))
    
    lsNoncodingGenesStict = lsGenes - lsCodingGenes
    fout = open('list.txt','w')
    fout.write('\n'.join(lsNoncodingGenesStict))
    fout.close()

    fout = open('list.txt','w')
    for ele in lsNt:
        fout.write(ele.description.split()[1]+' '+ele.id+'\n')
    fout.close()
    #No good, the file generated by cuffnorm are very bad. The average FPKM in each libraries are very different.

def noncodingGenesCufflinks20170609():
    fileNt = r"D:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\transcriptsAndTranslation\min60\20160713ManducaCufflinksTranscripts.fa"
    filepep = r"D:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\transcriptsAndTranslation\min60\20160713ManducaCufflinksTranscripts.fa.transdecoder.pep"
    import largeFastaFile
    lsNt = largeFastaFile.open_fasta_to_list(fileNt)
    lspep = largeFastaFile.open_fasta_to_list(filepep)
    lsCodingT = []
    for ele in lspep:
        codingName = ele.id.split(':')[2]
        lsCodingT.append(codingName)
    lsCodingT = list(set(lsCodingT))
    print(len(lsCodingT))
    lsNoncoding = []
    for ele in lsNt:
        if ele.id not in lsCodingT:
            lsNoncoding.append(ele)
    lsNoncoding = set(lsNoncoding)
    print(len(lsNoncoding))
    
    largeFastaFile.saveFastaListToFile(lsNoncoding,fileNt+'.noncodingT')
    lsNoncodingGenes = []
    for ele in lsNoncoding:
        geneID = ele.description.split()[1]
        if geneID not in lsNoncodingGenes:
            lsNoncodingGenes.append(geneID)
    print(len(lsNoncodingGenes))
    
    lsCodingGenes = set()
    lsGenes = set()
#    lsNoncodingT = [ele.id for ele in lsNoncoding]
    
    for ele in lsNt:
        geneID = ele.description.split()[1]
        if ele.id in lsCodingT:
            lsCodingGenes.add(geneID)
        lsGenes.add(geneID)
    print(len(lsCodingGenes), len(lsGenes))

def msRSEM20161019(ls_libraries):
    """
    Have to run RSEM for Cufflinks Gene models. Cuffnorm works bad.
    """
    for ele in ls_libraries:
        fo = open("20161019Ms_"+ele['ShortName']+"RSEM.txt","w")
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N TopCuff"+ele['ShortName']+"\n")
        fo.write("#PBS -l nodes=1:ppn=12\n")
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        fo.write("cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/RNA-trim/\n")
        fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/bowtie2-2.2.3/\n")
        fo.write("export PATH\n")
        fo.write("PATH=$PATH:/home/ks2073/p/20140314cufflinksNew/samtools-0.1.19/\n")
        fo.write("export PATH\n")
        fo.write("module load gcc-4.7.2\n")
        fo.write("mkdir /scratch/ks2073/RSEM/RSEMTemp/" + ele['ShortName'] + '\n')
        fo.write("/home/ks2073/p/rsem-1.2.15/rsem-calculate-expression  --bowtie2   --ci-memory 32000 -p 12 --no-bam-output ")
        if ele['File2'] == '':
            file_reads = ele['File1']
            fo.write(file_reads)
            
            
            
            
            
        else:
            file_reads_p1 = ele['File1']+"paired.fq"
            file_reads_p2 = ele['File2']+"paired.fq"
            file_reads_s1 = ele['File1']+"unpaired.fq"
            file_reads_s2 = ele['File2']+"unpaired.fq"
            fo.write(file_reads_p1 +"," + file_reads_p2+","+file_reads_s1+","+file_reads_s2)
                    
        fo.write(" /scratch/ks2073/RSEM/RSEMTemp/DNAseq/ProDB ")
        fo.write(" /scratch/ks2073/RSEM/RSEMTemp/" + ele['ShortName'] + '/'+ele['ShortName'] +'\n')
        fo.close()
    
    
def msRSEMdataCombine20161020(ls_libraries):
    '''
    combine the result to one table
    '''
    folderG = r"E:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\RSEM\GenesOri\\"
    filenames = []
    for ele in ls_libraries:
        filename = folderG + ele["ShortName"]+'.genes.results'
        filenames.append(filename)
    def getcsvColumn(filename,columnNum,separater = '\t'):
        '''
        return a list of columnNum given a filename
        '''
        templs = open(filename).readlines()
        columContent = []
        for ele in templs:
            columContent.append(ele.replace('\n','').split(separater)[columnNum])
        return columContent
    lsFPKMs = [getcsvColumn(filename,6) for filename in filenames]
    lsOtherInfo = [getcsvColumn(filename,num) for num in range(4)]
    fout = open('20161020CufflinksV4GeneFPKM.txt','w')
    lsAll = lsOtherInfo + lsFPKMs
    for num in range(len(lsFPKMs[0])):
        lineInfo = [lsAll[n][num] for n in range(len(lsAll))]
        fout.write('\t'.join(lineInfo)+'\n')
    fout.close()
    
    folderG = r"E:\Lab\works\20121023_manduca_genome_new\TranscriptAssembleNew\201607CufflinksV4NewGenome\RSEM\GenesOri\\"
    filenames = []
    for ele in ls_libraries:
        filename = folderG + ele["ShortName"]+'.genes.results'
        filenames.append(filename)
    def getcsvColumn(filename,columnNum,separater = '\t'):
        '''
        return a list of columnNum given a filename
        '''
        templs = open(filename).readlines()
        columContent = []
        for ele in templs:
            columContent.append(ele.replace('\n','').split(separater)[columnNum])
        return columContent
    lsFPKMs = [getcsvColumn(filename,4) for filename in filenames]
    lsOtherInfo = [getcsvColumn(filename,num) for num in range(4)]
    fout = open('20161020CufflinksV4GeneExpectedCount.txt','w')
    lsAll = lsOtherInfo + lsFPKMs
    for num in range(len(lsFPKMs[0])):
        lineInfo = [lsAll[n][num] for n in range(len(lsAll))]
        fout.write('\t'.join(lineInfo)+'\n')
    fout.close()
    
    folderG = "E:\\Lab\\works\\20121023_manduca_genome_new\\TranscriptAssembleNew\\201607CufflinksV4NewGenome\\RSEM\\IsoformsOri\\"
    filenames = []
    for ele in ls_libraries:
        filename = folderG + ele["ShortName"]+'.isoforms.results'
        filenames.append(filename)
    def getcsvColumn(filename,columnNum,separater = '\t'):
        '''
        return a list of columnNum given a filename
        '''
        templs = open(filename).readlines()
        columContent = []
        for ele in templs:
            columContent.append(ele.replace('\n','').split(separater)[columnNum])
        return columContent
    lsFPKMs = [getcsvColumn(filename,6) for filename in filenames]
    lsOtherInfo = [getcsvColumn(filename,num) for num in range(4)]
    fout = open('20161020CufflinksV4isoformsFPKM.txt','w')
    lsAll = lsOtherInfo + lsFPKMs
    for num in range(len(lsFPKMs[0])):
        lineInfo = [lsAll[n][num] for n in range(len(lsAll))]
        fout.write('\t'.join(lineInfo)+'\n')
    fout.close()
    
def topExpressedGenes20161026():
    filename = r"E:\Lab\works\2013ManducaUtilities\Transcriptome\RSEM\OGS2_20140430\20161026OGS2_FPKM4program.txt"
    fo = open(filename)
    fo.readline()
    lsls = [[] for i in range(67)]
    lsLines = []
    lsnames = []
    for ele in fo:
        lsLines.append(ele)
        templs = ele.split()
        lsnames.append(templs[0])
        for i in range(67):
            lsls[i].append(float(templs[i+1]))
    import numpy as np
    for i in range(67):
        lsls[i] = np.array(lsls[i])
    
    topNum = 3
    topIDs = set()
    for templs in lsls:
        topIDs = topIDs.union(set(templs.argsort()[::-1][:topNum]))
#            print([templs[i] for i in topIDs])
#            break
    print(len(topIDs))
    fout = open("list.txt",'w')
    topIDs = list(topIDs)
    topIDs.sort()
    for ele in topIDs:
        fout.write(lsLines[ele])
    fout.close()
    topnames = []
    for ele in topIDs:
        topnames.append(lsnames[ele])
    import largeFastaFile
    filenamepep = r"E:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2_2016\ms_ogs_proteins.fa"
    lspep = largeFastaFile.open_fasta_to_list(filenamepep)
    fout = open('seq.txt','w')
    for ele in lspep:
        if ele.id.split('-')[0] in topnames:
            fout.write('>'+ele.id+'\n'+str(ele.seq)+'\n')
    fout.close()

def goEnrichmentAnalysis20161027():
    import largeFastaFile
    filenamepep = r"E:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2_2016\ms_ogs_proteins.fa"
    lspep = largeFastaFile.open_fasta_to_list(filenamepep)
    #change list of gene names to list of transcript names for GO enrichment analysis
    set_genes = set(open("list.txt").read().split())
    fout = open("trancript.txt",'w')
    for ele in lspep:
        if ele.id.split('-')[0] in set_genes:
            fout.write(ele.id+'\n')
    fout.close()
    
    
def OGS2GeneNaming20170106():
    '''
    name OGS2 genes
    blast OGS2 proteins to NCBI nr. 
    download gene names from NCBI gene
    download gene id to accession id relation file
    '''
    f_ogs2blast = r"D:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2201404\GeneNaming\20170106OGS2blastNR_out.txt"
    import pandas as pd
    colnames = 'Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score'.split(', ')
    df_ogs2blast = pd.read_csv(f_ogs2blast,sep = '\t',header = None)
    df_ogs2blast.columns = colnames
    #Subject id looks like 'gi|768422635|ref|XP_011552185.1|', change to 'XP_011552185.1'
    df_ogs2blast['Subject id'] = df_ogs2blast['Subject id'].apply(lambda x: x.split('|')[3])  
    
    #filter, '% identity'>25, alignment length >50, e-value < 10^-6, bit score > 100
    df_ogs2blastFilter = df_ogs2blast[(df_ogs2blast['% identity'] > 25) & 
                                      (df_ogs2blast['alignment length'] > 50) & 
                                      (df_ogs2blast['e-value'] < 10**6) & 
                                      (df_ogs2blast['bit score'] > 100)]
    del df_ogs2blast
    
    #open gene accession file. store accessions that exist in df_ogs2blastFilter
    f_GeneAccession = r"D:\Insects\ncbi_geneNames\gene2accession"
    f_GeneAccessionMsexOGS2 = f_GeneAccession+'.MsexOGS2'
    fo = open(f_GeneAccession,'r')
    fout = open(f_GeneAccessionMsexOGS2,'w')
    line = fo.readline()
    fout.write(line)
    ogs2Accs = set(df_ogs2blastFilter['Subject id'])
    for line in fo:
        accession = line.split('\t')[5]
        if accession in ogs2Accs:
            fout.write(line)
    fout.close()
    fo.close()
    
    #open gene accession file of Manduca OGS2, store in dataframe
    f_GeneAccession = r"D:\Insects\ncbi_geneNames\gene2accession"
    f_GeneAccessionMsexOGS2 = f_GeneAccession+'.MsexOGS2'
    df_geneAcc = pd.read_csv(f_GeneAccessionMsexOGS2,sep = '\t')
    
    #open gene name file, store genes that exist in df_geneAcc
    ogs2Gene = set(df_geneAcc['GeneID'].astype(str))
    f_GeneName = r"F:\Insects\ncbi_geneNames\All_Data.gene_info"
    fo = open(f_GeneName,'r')
    fout = open(f_GeneName+'.MsexOGS2','w')
    line = fo.readline()
    fout.write(line)
    for line in fo:
        gene = line.split('\t')[1]
        if gene in ogs2Gene:
            fout.write(line)
    fout.close()
    fo.close()
    
    #open gene name file of Manduca OGS2, store in dataframe
    f_GeneName = r"D:\Insects\ncbi_geneNames\All_Data.gene_info"
    df_geneName = pd.read_csv(f_GeneName+'.MsexOGS2',sep = '\t')
    
    #combine df_geneAcc with df_geneName, get a dataframe with accession, geneID and geneName
    dftemp = pd.merge(df_geneAcc,df_geneName,on = 'GeneID',how = 'inner')
    df_accName = dftemp.loc[:,['GeneID','protein_accession.version','description']]
    df_accName = df_accName.drop_duplicates(subset = 'protein_accession.version')
    df_accName.columns = ['GeneID', 'Subject id', 'description']
    
    #give each transcript in df_ogs2blastFilter a name
    dftemp = pd.merge(df_ogs2blastFilter,df_accName,on = 'Subject id',how = 'left')
    df_transcript = dftemp.iloc[:,[0,1,2,3,10,11,12,13]]
    
    #add a column with OGS2 gene id
    df_transcript['OGS2ID'] = [ele.split('-')[0] for ele in df_transcript['Query id']]
    df_transcriptSort = df_transcript.sort(columns = ['OGS2ID','bit score'], ascending =[True, False])
    df_transcriptSortName = df_transcriptSort.dropna(axis = 0,subset = ['description'])
    df_transcriptBest = df_transcriptSortName.drop_duplicates(subset = 'OGS2ID')
    
    #save df_transcriptBest to file
    df_transcriptBest.to_csv('20170106OGS2NamingByNCBIgenes.csv')
        
    
    
    
    
    
    
    
def geneOntologyEnrichmentAnalysisTest20170111():
    '''
        #obodag is like a dictionary, with go like 'GO:0021833' as key and class GOTerm as element
    for element of GOTerm, element.id, name, namespace is the same as in the go-basic.obo file shown below.
    element.children, parent is the relationship of GO from the is_a statement
    
    id: GO:0021833
    name: cell-matrix adhesion involved in tangential migration using cell-cell interactions
    namespace: biological_process
    def: "The interaction of a cell and the extracellular matrix involved in the directed tangential movement of cells mediated by cell-cell interactions in the developing cerebral cortex." [GO_REF:0000021, GOC:ascb_2009, GOC:cls, GOC:dgh, GOC:dph, GOC:jid, GOC:mtg_15jun06, GOC:tb, PMID:12626695]
    comment: This term was added by GO_REF:0000021.
    is_a: GO:0031589 ! cell-substrate adhesion
    relationship: part_of GO:0021823 ! cerebral cortex tangential migration using cell-cell interactions
    '''
    from goatools.obo_parser import GODag
    f_gobasic = r"F:\Insects\20170111go-basic.obo"
    obodag = GODag(f_gobasic) #read in the go-basic file
    
    from goatools.associations import read_ncbi_gene2go
    geneid2gos_mouse = read_ncbi_gene2go(r"F:\Insects\20170111gene2go",taxids = [10090])
    print("{N:,} annotated mouse genes".format(N=len(geneid2gos_mouse)))
    #geneid2gos is a dictionary, looks like {328035:{'GO:0006629', 'GO:0006631', 'GO:0006633', 'GO:0016020', 'GO:0016021', 'GO:0016491', 'GO:0055114'},...}
    
    from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
    #GeneID2nt_mus is a dictionary, with geneID as key, and a class NtData as value.
    '''
    3. Initialize a GOEA object
    The GOEA object holds the Ontologies, Associations, and background.
    Numerous studies can then be run withough needing to re-load the above items.
    In this case, we only run one GOEA.
    '''
    from goatools.go_enrichment import GOEnrichmentStudy
    goeaobj = GOEnrichmentStudy(
        GeneID2nt_mus.keys(), # List of mouse protein-coding genes
        geneid2gos_mouse, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    geneid2symbol = {}
    din_xlsx = r"D:\P\3Language\Anaconda\Lib\site-packages\goatools-master\tests\data\nbt_3102\nbt.3102-S4_GeneIDs.xlsx"
    import xlrd
    book = xlrd.open_workbook(din_xlsx)
    pg = book.sheet_by_index(0)
    for r in range(pg.nrows):
        symbol, geneid, pval = [pg.cell_value(r,c) for c in range(pg.ncols)]
        if geneid:
            geneid2symbol[int(geneid)] = symbol
    
    geneids_study = geneid2symbol.keys()
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh <0.05]
    
    #save file
    goeaobj.wr_xlsx("nbt3102.xlsx", goea_results_sig)
    goeaobj.wr_txt("nbt3102.txt", goea_results_sig)

    
def geneOntologyEnrichmentAnalysisMsexta20170111():
    '''
    gene ontology enrichment analysis for Transcriptome paper Fig. 7
    '''
    from goatools.obo_parser import GODag
    f_gobasic = r"F:\Insects\20170111go-basic.obo"
    obodag = GODag(f_gobasic) #read in the go-basic file
    
    #read in GeneID to go file
    f_gene2go = r"E:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2201404\20170108OGS2Gene2GeneOntology.txt"
    fo = open(f_gene2go)
    from collections import defaultdict
    gene2go = defaultdict(set)
    fo.readline()
    for ele in fo:
        _gene, _go = ele[:-1].split('=')
        gene2go[_gene].add('GO:'+_go)
    fo.close()
    
    geneIDsAll = set()
    f_geneIDsall = r"E:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2201404\20160720ManducaOGS2Gene2Transcript.txt"
    fo = open(f_geneIDsall,'r')
    for ele in fo:
        _gene = ele.split()[0]
        geneIDsAll.add(_gene)
    fo.close()
    
    geneIDsHigh = set()
    geneIDsHighGroups = {}
    f_geneIDgroups = r"E:\Lab\works\2013ManducaUtilities\Transcriptome\RSEM\OGS2_20140430\ForTranscriptomePaper20161011\2017OGS2TopExpressedGeneIDtoGroups.txt"
    fo = open(f_geneIDgroups,'r')
    for ele in fo:
        _gene,_group = ele.split()
        if _group not in geneIDsHighGroups:
            geneIDsHighGroups[_group] = set()
        geneIDsHigh.add(_gene)
        geneIDsHighGroups[_group].add(_gene)
    fo.close()
    
    #Initialize a GOEA object
    from goatools.go_enrichment import GOEnrichmentStudy
    goeaobj = GOEnrichmentStudy(
        geneIDsAll, # List of all genes
        gene2go, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    #test geneIDsHigh against geneIDsAll
    goea_results_all = goeaobj.run_study(geneIDsHigh)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh <0.01]
    #save file
    goeaobj.wr_xlsx("20170116HighExpressedGenesEnrichedAgainstAll.xlsx", goea_results_sig)
    
    #Initialize a GOEA object, use geneIDsHigh as background
    goeaobj2 = GOEnrichmentStudy(
        geneIDsHigh, # List of all genes
        gene2go, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    #test geneIDsHighGroups against goeaobj2
    dc_goea_results_all = {}
    dc_goea_results_sig = {}
    for ele in geneIDsHighGroups:
        dc_goea_results_all[ele] = goeaobj2.run_study(geneIDsHighGroups[ele])
        dc_goea_results_sig[ele] = [r for r in dc_goea_results_all[ele] if (r.p_fdr_bh <0.05 and r.enrichment =='e')]
        goeaobj.wr_xlsx("20170116TissueGenesEnrichedAgainstAllHigh_group"+ ele +".xlsx", dc_goea_results_sig[ele])
    
    fout = open("20170116TissueGenesEnrichedAgainstAllHigh_groupAll.txt",'w')
    for ele in dc_goea_results_sig:
        for ele2 in dc_goea_results_sig[ele]:
            fout.write(ele +'\t'+ele2.__str__()+'\n')
    fout.close()

def getHiearchicalClusterTree20170629():
    '''
    extract tree for Msex2
    '''
    ls_cdt = open(r"C:\Users\k\Downloads\20161023OGS2MCOTonlyCufflinksV4FPKMmin100zCombine.cdt").readlines()
    ls_gtr = open(r"C:\Users\k\Downloads\20161023OGS2MCOTonlyCufflinksV4FPKMmin100zCombine.gtr").readlines()
    ls_cdtM = [e for e in ls_cdt if 'Msex2.' in e]
    print('len of total',len(ls_cdt),'len of ls_cdtM',len(ls_cdtM))
    
    open(r"C:\Users\k\Downloads\20161023OGS2MCOTonlyCufflinksV4FPKMmin100zCombine2.cdt",'w').write(''.join(ls_cdtM))
    
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\Downloads\20161023OGS2MCOTonlyCufflinksV4FPKMmin100zCombine.gtr",sep='\t',header=None)
    nodesall = [e.split()[0] for e in ls_cdt if e[:4] == 'GENE']
    nodeskeep = [e.split()[0] for e in ls_cdtM if e[:4] == 'GENE']
    nodesRemove = [e for e in nodesall if e not in nodeskeep]
    nodesRemoved =[]
    for e in nodesRemove:
        line = (df[df.loc[:,1]==e].index.tolist() +df[df.loc[:,2]==e].index.tolist())[0]
        node = df.loc[line,0]
        anothernode = [_i for _i in df.loc[line,[1,2]] if _i != e][0]
#        linelen = df.loc[line,3]
        df2 = df.drop(line,axis=0)
        line2 = (df2[df2.loc[:,1]==node].index.tolist() +df2[df2.loc[:,2]==node].index.tolist())[0]
#        df2.loc[line2,3] = df2.loc[line2,3]+linelen
        for n in [1,2]:
            if df2.loc[line2,n] == node:
                df2.loc[line2,n] = anothernode
        df = df2
        nodesRemoved.append(e)
#        if len(nodesRemoved)>2000:
#            break
    fout = open(r"C:\Users\k\Downloads\20161023OGS2MCOTonlyCufflinksV4FPKMmin100zCombine3.cdt",'w')
    for e in ls_cdt:
        if e.split()[0] not in nodesRemoved:
            fout.write(e)
    fout.close()
    df.to_csv(r"C:\Users\k\Downloads\20161023OGS2MCOTonlyCufflinksV4FPKMmin100zCombine3.gtr",sep='\t',header=False,index=False)
    
    
    
    
def geneOntologyEnrichmentAnalysisMsexta2017630():
    '''
    gene ontology enrichment analysis for Transcriptome paper Fig. 7
    '''
    from goatools.obo_parser import GODag
    f_gobasic = r"D:\Insects\20170111go-basic.obo"
    obodag = GODag(f_gobasic) #read in the go-basic file
    
    #read in GeneID to go file
    f_gene2go = r"D:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2201404\20170108OGS2Gene2GeneOntology.txt"
    fo = open(f_gene2go)
    from collections import defaultdict
    gene2go = defaultdict(set)
    fo.readline()
    for ele in fo:
        _gene, _go = ele[:-1].split('=')
        gene2go[_gene].add('GO:'+_go)
    fo.close()
    
    geneIDsAll = set()
    f_geneIDsall = r"D:\Lab\works\20121023_manduca_genome_new\20120818Genome Ms\OGS2201404\20160720ManducaOGS2Gene2Transcript.txt"
    fo = open(f_geneIDsall,'r')
    for ele in fo:
        _gene = ele.split()[0]
        geneIDsAll.add(_gene)
    fo.close()
    
    geneIDsHigh = set()
    geneIDsHighGroups = {}
    f_geneIDgroups = r"D:\Lab\works\2013ManducaUtilities\Transcriptome\RSEM\OGS2_20140430\ForTranscriptomePaper20161011\2017OGS2TopExpressedGeneIDtoGroups2.txt"
    fo = open(f_geneIDgroups,'r')
    for ele in fo:
        _gene,_group = ele.split()
        if _group not in geneIDsHighGroups:
            geneIDsHighGroups[_group] = set()
        geneIDsHigh.add(_gene)
        geneIDsHighGroups[_group].add(_gene)
    fo.close()
    
    #Initialize a GOEA object
    from goatools.go_enrichment import GOEnrichmentStudy
    goeaobj = GOEnrichmentStudy(
        geneIDsAll, # List of all genes
        gene2go, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    #test geneIDsHigh against geneIDsAll
    goea_results_all = goeaobj.run_study(geneIDsHigh)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh <0.01]
    #save file
    goeaobj.wr_xlsx("20170630HighExpressedGenesEnrichedAgainstAll.xlsx", goea_results_sig)
    
    #Initialize a GOEA object, use geneIDsHigh as background
    goeaobj2 = GOEnrichmentStudy(
        geneIDsHigh, # List of all genes
        gene2go, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    #test geneIDsHighGroups against goeaobj2
    dc_goea_results_all = {}
    dc_goea_results_sig = {}
    for ele in geneIDsHighGroups:
        dc_goea_results_all[ele] = goeaobj2.run_study(geneIDsHighGroups[ele])
        dc_goea_results_sig[ele] = [r for r in dc_goea_results_all[ele] if (r.p_fdr_bh <0.05 and r.enrichment =='e')]
        goeaobj.wr_xlsx("20170630TissueGenesEnrichedAgainstAllHigh_group"+ ele +".xlsx", dc_goea_results_sig[ele])
    
    fout = open("20170630TissueGenesEnrichedAgainstAllHigh_groupAll.txt",'w')
    for ele in dc_goea_results_sig:
        for ele2 in dc_goea_results_sig[ele]:
            fout.write(ele +'\t'+ele2.__str__()+'\n')
    fout.close()

def getInforForhumanGenome20170709():
    '''
    get some information from the human genome
    '''
    #get genome length
    from Bio import SeqIO
    f_hsgenome = r"D:\Insects\human\2017RefSeq\GRCh38_latest_genomic.fna"
    l_hsgenome = 0
    for e in SeqIO.parse(open(f_hsgenome),'fasta'):
        l_hsgenome += len(e.seq)
    print(l_hsgenome)
    
    #get gene length
    f_hsGTF = r"D:\Insects\human\2017RefSeq\GRCh38_latest_genomic.gff"
    filename = f_hsGTF
    fo = open(filename,'r')
    lslsgff = []
    targetwords = 'CDS'.split()
    for ele in fo:
        if ele == '\n' or ele[0] =='#':
            continue
        eles = ele.split()
        if eles[2] in targetwords:
            lslsgff.append([eles[0], int(eles[3]), int(eles[4])])
    fo.close()
    dcContigGff = {}
    for ele in lslsgff:
        if ele[0] not in dcContigGff:
            dcContigGff[ele[0]] = []
        dcContigGff[ele[0]].append(ele[1:])
    dcContigTranslength = {}
    for ele in dcContigGff:
        dcContigTranslength[ele] = set()
        for eleStart, eleEnd in dcContigGff[ele]:
            dcContigTranslength[ele].update(list(range(eleStart,eleEnd+1)))
        dcContigTranslength[ele] = len(dcContigTranslength[ele])
    print(sum(dcContigTranslength.values()))