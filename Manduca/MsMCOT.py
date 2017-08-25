# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 15:42:21 2015

@author: k
"""

def removeSeqsInFile(filename,listname):
    """
    given a file of fasta format, and a list file is names to be removed
    output a filename_out file, to remove those un-wanted sequences
    """
    from Bio import SeqIO
    mylist = open(listname).read().split()
    myfas = list(SeqIO.parse(filename,"fasta"))
    fo = open(filename+"_out","w")
    for ele in myfas:
        if ele.id not in mylist:
            fo.write(">"+ele.description+"\n"+str(ele.seq)+"\n")
    fo.close()
    return None


filename = r"E:\store_for_D\Transcriptome\20150105MCOT31166proteinSequenceNewIDfasta.txt"
listname = "list.txt"
removeSeqsInFile(filename, listname)


folder = "E:\store_for_D\Transcriptome\\"
filename = "20150105MCOT30303proteinSequenceNewIDfastaRemoveOtherSpecies.txt"
from Bio import SeqIO
myfasta = list(SeqIO.parse(open(folder+filename),"fasta"))
for i in range(0,len(myfasta),5000):
    fo = open(folder+"part"+str(i)+".txt","w")
    for j in range(i,min(i+5000,len(myfasta))):
        SeqIO.write(myfasta[j],fo,"fasta")
    fo.close()

def MCOT20160404():
    folder = "E:\\store_for_D\\Transcriptome\\"
    file_mcot = "20150105MCOT30303proteinSequenceNewIDfastaRemoveOtherSpecies.txt"
    from Bio import SeqIO
    