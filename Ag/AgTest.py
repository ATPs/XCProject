# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 13:10:11 2016

@author: k
"""

def AgHaobo100SPSPHcheck():
    from Bio import SeqIO
    lsSP100 = list(SeqIO.parse(open(r'D:\mine\OneDrive\Lab\Jiang Lab Pubilic\ProteomeAg\20160208Ag100SPSPH\20160324AgSPSPH100Haobo.txt'),"fasta"))
    dcAgP42 = SeqIO.to_dict(SeqIO.parse(open(r"D:\mine\OneDrive\Lab\Jiang Lab Pubilic\ProteomeAg\PEST_PEPTIDES.fasta"),"fasta"))
    dcAgP43 = SeqIO.to_dict(SeqIO.parse(open(r"F:\Insects\Anopheles_gambiae\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"),"fasta"))
    ls42only=[]
    ls42diff =[]
    for key in dcAgP42:
        if key not in dcAgP43:
            ls42only.append(key)
        else:
            if dcAgP42[key].seq != dcAgP43[key].seq:
                ls42diff.append(key)
    print(len(ls42only),len(ls42diff))
    
    fout = open("AgSPsCompare.txt","w")
#    fout2 = open("AgSPsCompareDiff.txt","w")
    for ele in lsSP100:
        key = ele.id
        seq1 = str(ele.seq)
        fout.write(">%s\n"%key)
        fout.write(seq1+"\n")
        if key in dcAgP42:
            seq2 = str(dcAgP42[key].seq)
            if seq1 == seq2:
                fout.write("AgamP4.2, same\n")
            else:
                fout.write("AgamP4.2, different\n"+seq2+"\n")
#                fout2.write(">"+key+"\n"+seq1+"\n"+"AgamP4.2, different\n"+seq2+"\n")
        else:
            fout.write("not found in AgamP4.2\n")
        if key in dcAgP43:
            seq3 = str(dcAgP43[key].seq)
            if seq1 == seq3:
                fout.write("AgamP4.3, same\n")
            else:
                fout.write("AgamP4.3, different\n"+seq3+"\n")
        else:
            fout.write("not found in AgamP4.3\n")
    fout.close()
#    fout2.close()

def AgTrinityRename(filename,keywords="",outfile=None):
    """
    name looks like TRINITY_DN59_c0_g1_i1, change to keywords+geneNumber+isoformNumber
    this one, it may change to 58_i1
    if outfile is None, outfile = filename + "rename"
    """
    if outfile is None:
        outfile = filename+"rename"
    from Bio import SeqIO
    fo = open(filename,"r",1000000000)
    fout = open(outfile,"w",1000000000)
    genes = {}
    for seq in SeqIO.parse(fo,"fasta"):
        genename = seq.id.rsplit("_",1)[0]
        isoformname = seq.id.rsplit("_",1)[1]
        if genename not in genes:
            genes[genename] = 0
            genes[genename] = keywords + str(len(genes))
        fout.write(">"+genes[genename]+isoformname+"\n"+str(seq.seq)+"\n")
    fo.close()
    fout.close()
#AgTrinityRename("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/XTrinity.Trinity.fasta","O")
#AgTrinityRename("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/2LTrinity.Trinity.fasta","2L")
#AgTrinityRename("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/2RTrinity.Trinity.fasta","2R")
#AgTrinityRename("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/3LTrinity.Trinity.fasta","3L")
#AgTrinityRename("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/3RTrinity.Trinity.fasta","3R")

def AglongpepSplit():
    filename ='/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/20160406TrinityCombined98.fa.transdecoder_dir/longest_orfs.pep'
    from Bio import SeqIO
    num = 100
    filedic ={}
    scriptsdic = {}
    towrite="""#!/bin/bash

#PBS -q batch 
#PBS -N Transdecoder
#name you want to give your job 
#PBS -l nodes=1:ppn=12
#request 1 nodes w/8 processors per node
#PBS -l walltime=12:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe


cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity

cp /panfs/panfs.cluster/opt/interproscan/5.17-56.0-64-bit/prebuilt/data/pfam/28.0/pfam_a.hmm /dev/shm/pfam_a.hmm
cp /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/20160406TrinityCombined98.fa.transdecoder_dir/longest_orfs.pep%d /dev/shm/longest_orfs.pep%d
module load hmmer
module load blast+

hmmpress /dev/shm/pfam_a.hmm

hmmscan --cpu 1 --domtblout pfam.domtblout%d /dev/shm/pfam_a.hmm /dev/shm/longest_orfs.pep%d >>hmmerlog.txt%d 
blastp -db /panfs/panfs.cluster/scratch/ks2073/Blast/blastDB/2016Athropoda/uniprotkb_arthropoda -query /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/20160406TrinityCombined98.fa.transdecoder_dir/longest_orfs.pep%d -out /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Trinity/20160406TrinityCombined98.fa2Arthropodaformat6.txt%d  -evalue 0.000001 -outfmt 6 -max_target_seqs 5 -num_threads 12 -task blastp-fast -threshold 23 


    """
    for i in range(num):
        filedic[i] = open(filename + str(i),"w")
        scriptsdic[i] = open(filename + str(i)+".txt","w")
        scriptsdic[i].write(towrite%(i,i,i,i,i,i,i))
    count =0
    for seq in SeqIO.parse(open(filename),"fasta"):
        count +=1
        filedic[count % num].write(">"+seq.id+"\n"+str(seq.seq)+"\n")
    for i in range(num):
        filedic[i].close()
        scriptsdic[i].close()
AglongpepSplit()



def getAgProtein20160415():
    from Bio import SeqIO
    lsfa = list(SeqIO.parse(open("seqs.txt"),"fasta"))
    fout = open("seqout.txt","w")
    for ele in lsfa:
        fout.write(">"+ele.id+"\n"+str(ele.seq)+"\n")
    fout.close()
    ls_name = open("list.txt").read().split("\n")
    dc_fa = SeqIO.to_dict(SeqIO.parse(open("seqout.txt"),"fasta"))
    fout = open("seqout2.txt","w")
    for ele in ls_name:
        fout.write(">"+ele+"\n"+str(dc_fa[ele].seq)+"\n")
    fout.close()


def AgProteome():
    fo = open(r'F:\Insects\AgAdultMassSpect\proteinGroups.txt')
    mylist = fo.readlines()
    fout = open(r'F:\Insects\AgAdultMassSpect\proteinGroupsModify.txt','w')
    for ele in mylist:
        lse = ele.split('\t')
        for num in range(len(lse)):
            if num <55+26:
                fout.write(lse[num]+'\t')
        fout.write("\n")
                    
    fout.close()
    fo.close()


#20160708
def agGeneNames():
    '''
    get a file with protein ID and its names from the OGS
    '''
    fo = open('F:\Insects\Anopheles_gambiae\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa')
    fout = open('F:\Insects\Anopheles_gambiae\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3_ID2Name.txt','w')
    mylist = fo.readlines()
    for line in mylist:
        if line[0] == '>':
            pep_id = line[1:].split(' ')[0]
            pep_name = line[1:].split(' ',1)[1].split('|')[0]
            fout.write(pep_id +'\t'+pep_name+'\n')
    fout.close()
    fo.close()
        