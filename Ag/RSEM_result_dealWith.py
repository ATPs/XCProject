# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:28:58 2015

@author: k
"""

folder = "D:\\Insects\\Anopheles_gambiae\\SPSPH\\RSEM\\"

import os
filelist = os.listdir(folder)
filelist = [i for i in filelist if "isoforms" in i]

def fileTable2ListOfList(filename, separator = None):
    """
    given a file, return a list of list
    each line is split by the default separators of python, if separator is not provided
    """
    mylist = open(filename).readlines()
    mylist2 = []
    if separator == None:
        for ele in mylist:
            mylist2.append(ele.split())
    else:
        for ele in mylist:
            mylist2.append(ele.split(separator))
    return mylist2

def extractColumnFromListOfList(mylist, column):
    """
    given a list of list, return a values of a given column, in the list format
    column number begins from 0
    """
    myreturn =[]
    for ele in mylist:
        myreturn.append(ele[column])
    return myreturn

def writeListOflist(myLoL, filename):
    """
    given a list of list, like below
    transcript_id	gene_id	length	effective_length	expected_count	TPM	FPKM	IsoPct
    Trinity100010|m.28666	Trinity100010|m.28666	894	845.52	215707.26	1094.61	1244.81	100.00
    Trinity100011|m.28669	Trinity100011|m.28669	336	287.52	23.22	0.35	0.39	100.00
    Trinity1000179|m.274915	Trinity1000179|m.274915	354	305.52	0.00	0.00	0.00	0.00
    Trinity1000344|m.274931	Trinity1000344|m.274931	327	278.52	14.96	0.23	0.26	100.00
    each column is stored in a list in myLoL, the out put is the same as above
    output to filename
    """
    fo = open(filename, "w")
    for num in range(len(myLoL[0])):
        towrite = ""
        for num2 in range(len(myLoL)):
            towrite += myLoL[num2][num]+"\t"
        towrite += "\n"
        fo.write(towrite)
    fo.close()

myTowrite=[]
myTowrite.append([])
myTowrite.append([])
for myfile in filelist:
    mylist = fileTable2ListOfList(folder+myfile)
    myTowrite.append(extractColumnFromListOfList(mylist, 6))
    myTowrite[-1][0] += myfile
myTowrite[0] = extractColumnFromListOfList(mylist, 0)
myTowrite[1] = extractColumnFromListOfList(mylist, 2)

myTowrite2 = []
myTowrite2.append(myTowrite[0])
myTowrite2.append(myTowrite[1])

fileorder1 = "ERR505238	ERR505234	ERR593540	ERR505235	ERR505236	ERR489305	ERR489300	ERR489298	ERR593537	ERR537786	ERR537790	ERR489297	ERR593539	ERR537778	ERR476715	ERR476718	ERR593542	ERR476720	ERR489301	ERR537785	ERR537784	ERR505239	ERR537779	ERR593541	ERR505237	ERR505240	ERR489303	ERR489299	ERR489302	ERR593538	ERR537782	ERR537787	ERR537781	ERR489304	ERR593536	ERR476717	ERR476719	ERR593543	ERR476716	ERR489306	ERR537788".split()
fileorder2 = "SRR513191	SRR513197	SRR513203	SRR513192	SRR513198	SRR513204	SRR513193	SRR513199	SRR513205	SRR513189	SRR513195	SRR513201	SRR513190	SRR513196	SRR513202	SRR513194	SRR513200	SRR513206	SRR1727555	SRR1727718	SRR1727719	SRR1727720	SRR1727721	SRR1727722	SRR1727723	SRR1727724	SRR1727725	SRR1727726	SRR520427	SRR520428	SRR1179005	SRR1179006	SRR1179007	SRR1179008	SRR1179009	SRR1179010	SRR1179011	SRR1179012".split()
fileorder3 = "ERR588643	ERR588645	ERR588651	ERR588639	ERR588650	ERR588659	ERR588664	ERR588647	ERR588654	ERR588667	ERR588642	ERR588657	ERR588669	ERR588666	ERR588663	ERR588649	ERR588655	ERR588641	ERR588648	ERR588661	ERR588653	ERR588660	ERR588656	ERR588662	ERR588640	ERR588646	ERR588652	ERR588644	ERR588658	ERR588665	ERR588638	ERR588668 ERR840657	ERR840656".split()
fileorder = fileorder1 + fileorder2+fileorder3
for ele in fileorder:
    for templist in myTowrite:
        if ele in templist[0]:
            myTowrite2.append(templist)
writeListOflist(myTowrite2,"AgamP4.3isoformsRSEM_FPKMall.txt")

def AgSPSPHs():
    folder = "F:\\Insects\\Anopheles_gambiae\\SPSPH\\"
    from Bio import SeqIO
    lsO = list(SeqIO.parse(open(folder+"SPSPH_aa.txt"),"fasta"))
    lsMCOT = list(SeqIO.parse(open(folder+"MCOTSPSPH_aa.txt"),"fasta"))
    print(len(lsO),len(lsMCOT))
    
    lsAll = []
    for ele1 in lsMCOT:
        seq1 = str(ele1.seq)
        if seq1[-1] == "*":
            seq1 = seq1[:-1]
        seq1in = False
        for ele2 in lsO:
            seq2 = str(ele2.seq)
            if seq1 == seq2:
                seq1in = True
                break
        if not seq1in:
            lsAll.append(ele1)
    print(len(lsAll))
    
    lsAll =  lsO + lsAll 
    print(len(lsAll))
    
    import re    #find all AGAP ids in lsAll
    lsAGAP =[]
    for ele in lsAll: 
        description = ele.description.split("uni:")[0]
        lsAGAP += re.findall("AGAP\d*-\w{2}",description)
    print(len(lsAGAP))
    lsAGAP = list(set(lsAGAP))
    lsAGAP.sort()
    print(len(lsAGAP))
    
    fout = open(folder+"SPSPH_combinedAA.txt","w")
    lsAllfinished =[]
    for ele in lsAll:
        description = ele.description.split("uni:")[0]
        agapids = re.findall("AGAP\d*-\w{2}",description)
        if len(agapids) == 0:
            fout.write(">"+ele.description+"\n"+str(ele.seq)+"\n")
            lsAllfinished.append(ele)
        else:
            agapid = agapids[0]
            if agapid not in lsAGAP:
                fout.write(">"+ele.description+"\n"+str(ele.seq)+"\n")
                lsAllfinished.append(ele)
#    fout.close()
    for ele2 in lsAGAP:
        for ele in lsAll:
            if ele not in lsAllfinished:
                description = ele.description.split("uni:")[0]
                agapid = re.findall("AGAP\d*-\w{2}",description)[0]
                if agapid == ele2:
                    fout.write(">"+ele.description+"\n"+str(ele.seq)+"\n")
                    lsAllfinished.append(ele)
    fout.close()

    
#    lsAll = list(SeqIO.parse(open(folder+"SPSPH_combinedAA.txt"),"fasta"))
#    lsSeq =[]
#    for ele in lsAll:
#        seq = str(ele.seq)
#        if seq[-1] == "*":
#            seq = seq[:-1]
#        lsSeq.append(seq)
    
#    import XCBlastpLike
#    dcKmer = XCBlastpLike.seqs2kmerdic(lsSeq,30)
#    tempout = XCBlastpLike.seqFindtargets(lsSeq[0],dcKmer)[:10]
#    fout = open(folder+"checke.txt","w")
#    for ele in tempout:
#        ele = lsAll[ele[0]]
#        fout.write(">"+ele.id+"\n"+str(ele.seq)+"\n")
#    fout.close()
#    
#    lsSort =[]
#    for num in range(len(lsSeq)):
        
    


def AgSPSPHgetExpression(f_list, fasta = False,outfile = "fpkm_selected.txt"):
    """
    given a list of ids, get expression data
    """
    folder = "D:\\Insects\\Anopheles_gambiae\\RSEM\\"
    f_mcot = "MCOTisoformsRSEM_FPKMall.txt"
    f_ogs = "AgamP4.3isoformsRSEM_FPKM_all.txt"
    if not fasta:
        ls_ids = open(f_list).read().split()
    else:
        from Bio import SeqIO
        ls_fas = list(SeqIO.parse(open(f_list),"fasta"))
        ls_ids =[i.id for i in ls_fas]
    ls_fpkm = open(folder+f_mcot).readlines() + open(folder+f_ogs).readlines()
    dc_fpkm ={}
    for ele in ls_fpkm:
        if ele[:4] == "MCOT":
            key = ele.split()[0]
            dc_fpkm[key] = ele
        elif ele[:4] == "AGAP":
            key = ele.split()[0].replace("R","P")
            dc_fpkm[key] = ele
    fout = open(outfile,"w")
    for ele in ls_ids:
        fout.write(dc_fpkm[ele])
    fout.close()

f_list = r"F:\Insects\Anopheles_gambiae\SPSPH\SPSPH_combinedAA.txt"
fasta = True
AgSPSPHgetExpression(r"F:\Insects\Anopheles_gambiae\SPSPH\SPSPH_combinedAA.txt",True)  
    
ls = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/test.fq").readlines()
fout = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq/test.fa","w")
for num in range(0,len(ls),4):
    fout.write(">"+ls[num][1:]+ls[num+1])
fout.close()

def getBasecoverage():
    filename = "/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Aligned.sort.bam.pileup"
    fout = open("/panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/STAR/Aligned.sort.bam.pileup.des","w")
    fout.write("#ch\tbase\tA\tC\tG\tT\n")
    with open(filename) as f:
        for line in f:
            linsplit = line.split()
            fout.write(linsplit[0]+"\t"+linsplit[1]+"\t%d\t%d\t%d\t%d\n"%(\
            linsplit[4].count('A')+linsplit[4].count('a'),\
            linsplit[4].count('C')+linsplit[4].count('c'),\
            linsplit[4].count('G')+linsplit[4].count('g'),\
            linsplit[4].count('T')+linsplit[4].count('t')))
    fout.close()
getBasecoverage()
            