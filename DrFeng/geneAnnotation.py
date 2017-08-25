# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:03:03 2016

@author: k
"""
folder = "E:\\store_for_D\\DrFeng\Assemblies\\TrinityTranslation\\"
fileBlast = "BlastTrinityPepNoneWithin2Nrformat6.txt"

def getGIs(fileBlast6,uniquekeep = 1, fileout = "gis.txt"):
    """
    fileBlast6 is the file in Blast output format 6.
    for each query, only keep the top uniquekeep number of hits.
    output gis.txt to store unique gi numbers
    """
    mylist = open(fileBlast6,"r").readlines()
    querys = {}
    giset = set()
    fout = open(fileout,"w")
    for ele in mylist:
        gi = ele.split()[1].split("|")[1]
        query = ele.split()[0]
        if query not in querys:
            querys[query] = 1
        else:
            querys[query] += 1
        if querys[query] <= uniquekeep:
            giset.add(gi)
    for gi in giset:
        fout.write(gi+"\n")
    fout.close()

getGIs(folder+fileBlast, 1, folder+"gis.txt")
        
def getFasta(fileGI,fileout = "gis.fasta", outfmt = "fasta"):
    """
    given a file of GIs, retrieve the information from NCBI, store the result
    in format of outfmt in fileout
    """
    myGIs = open(fileGI).read().split()
    gilist = [",".join(myGIs[i:i+500]) for i in range(0,len(myGIs),500)]
    from Bio import Entrez
    import time
    fout = open(fileout,"w")
    Entrez.email = "ks2074@gmail.com"
    for ele in gilist:
        handle = Entrez.efetch(db = "protein", id = ele, rettype = outfmt, retmode = "text")
        fout.write(handle.read())
        time.sleep(3)
    fout.close()

getFasta(folder+"gis.txt", folder+"gis.fasta")

def TrinityAnnotation(Trinityfa,fileBlast6,GIfa,output):
    """
    output a file, with TrinityID, and most similar gene
    """
    from Bio import SeqIO
    myTrinity = list(SeqIO.parse(Trinityfa,"fasta"))
    mygislist = list(SeqIO.parse(GIfa,"fasta"))
    mygis = {}
    for giseq in mygislist:
        mygis[giseq.id.split("|")[1]] = giseq
    mylist = open(fileBlast6,"r").readlines()
    fout = open(output,"w")
    sqdic = {}
    queryset = set()
    for ele in mylist:
        subject = ele.split()[1] +"|" + ele.split()[3]+" Identity: " + ele.split()[2]
        query = ele.split()[0]
        if query not in queryset:
            sqdic[query]=subject
            queryset.add(query)
    for ele in myTrinity:
        if ele.id not in queryset:
            fout.write(ele.id +"\t"+"NA\n")
        else:
            fout.write(ele.id +"\t"+mygis[sqdic[ele.id].split("|")[1]].description+" QL:"+\
            str(len(ele.seq)) +" SL: "+str(len(mygis[sqdic[ele.id].split("|")[1]].seq)) +\
            " ML: "+ sqdic[ele.id].split("|")[-1] +"\n")
    fout.close()

TrinityAnnotation(Trinityfa=folder+"TrinityPepNoneWithin2.fasta", \
fileBlast6=folder+fileBlast, GIfa=folder+"gis2.fasta", output = folder+"Trinity2GI.txt")
        
def fa2tab(filename):
    """
    given a file in fasta format, output filename.tab, sequence stored in tab format
    """
    from Bio import SeqIO
    fo = open(filename,"r")
    fout = open(filename+".tab","w")
    for seq in SeqIO.parse(fo, "fasta"):
        SeqIO.write(seq,fout,"tab")
    fo.close()
    fout.close()

fa2tab(folder+"TrinityPepNoneWithin2.fasta")