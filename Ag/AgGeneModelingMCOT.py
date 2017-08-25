# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 12:50:50 2016

@author: k
"""
from Bio import SeqIO

fastafiledic ={}
fastafiledic["O"] = "E:\\store_for_D\\Insects Info\\Anopheles gambiae\\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"
folder = "E:\\store_for_D\\Insects Info\\Anopheles gambiae\\Assemblies\\"
fastafiledic["C"] = folder + "AgCufflinks20160301.fa.transdecoder.pep"
fastafiledic["TS"] = folder + "TrinitySS.Trinity.fasta.transdecoder.pep"
fastafiledic["BS"] = folder + "BridgerSS.Bridger.fasta.transdecoder.pep"
fastafiledic["TU"] = folder + "TrinityUnmapped.Trinity.fasta.transdecoder.pep"
fastafiledic["BU"] = folder + "BridgerUnmap.Bridger.fasta.transdecoder.pep"
fastafiledic["NC"] = folder + "AgCufflinks20160310.fa.transdecoder.pep"

def lenCount(genelengths):
    """
    return numbers of genes in length range <128	128-511	512-1023	1024-2047	2048-4095	4096-8191	8092-16383	>=16384
    """
    base =32
    return (len([i for i in genelengths if i <base*2]),\
    len([i for i in genelengths if i >=base*2 and i<base*4]),\
    len([i for i in genelengths if i >=base*4 and i<base*8]),\
    len([i for i in genelengths if i >=base*8 and i<base*16]),\
    len([i for i in genelengths if i >=base*16 and i<base*32]),\
    len([i for i in genelengths if i >=base*32 and i<base*64]),\
    len([i for i in genelengths if i >=base*64 and i<base*128]),\
    len([i for i in genelengths if i >=base*128 and i<base*256]),\
    len([i for i in genelengths if i >=base*256]) )

#Count AgOGS gene transcripts number
filename = "E:\\store_for_D\\Insects Info\\Anopheles gambiae\\Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.3.fa"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.id.split("-")[0]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #13793
print(len(transcriptnames)) #15656


#Count Cufflinks gene transcripts number
filename = folder + "AgCufflinks20160301.fa"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.description.split()[1]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #37033
print(len(transcriptnames)) #66110

#Count NewCufflinks gene transcripts number
filename = folder + "AgCufflinks20160310.fa"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.description.split()[1]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #36019
print(len(transcriptnames)) #74559


#Count TrinitySS gene transcripts number
filename = folder + "TrinitySS.Trinity.fasta"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.id.rsplit("_",1)[0]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #338372
print(len(transcriptnames)) #480264

#Count BridgerSS gene transcripts number
filename = folder + "BridgerSS.Bridger.fasta"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.id.rsplit("_",1)[0]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #197561
print(len(transcriptnames)) #339635

#Count TrinityUnmap gene transcripts number
filename = folder + "TrinityUnmapped.Trinity.fasta"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.id.rsplit("_",1)[0]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #82861
print(len(transcriptnames)) #95386

#Count BridgerUnmap gene transcripts number
filename = folder + "BridgerUnmap.Bridger.fasta"
mylist = list(SeqIO.parse(open(filename),"fasta"))
genenames = set()
transcriptnames = set()
genelen ={}
for ele in mylist:
    transcriptnames.add(ele.id)
    genename = ele.id.rsplit("_",1)[0]
    genenames.add(genename)
    if genename not in genelen:
        genelen[genename] = len(ele.seq)
    else:
        if len(ele.seq) > genelen[genename]:
            genelen[genename] = len(ele.seq)
genelengths = list(genelen.values())
print("average gene length is ", sum(genelengths)/len(genelengths))
print(lenCount(genelengths))
print(len(genenames)) #86648
print(len(transcriptnames)) #132764

#Count number of proteins
for key in fastafiledic:
    mylist = list(SeqIO.parse(open(fastafiledic[key]),"fasta"))
    print(key, len(mylist))

#Count non-redundant proteins. Save non-redundant protein to new file. 
#non-redundant: protein is not the same with each other
#aligned together, the 2% mismatch allowed, no protein is part of other protein.
import largeFastaFile
mylisdic ={}
for key in ["O"]:
    mylist = list(SeqIO.parse(open(fastafiledic[key]),"fasta"))
    mylistn = largeFastaFile.fasta_within_seq_big_withError(mylist)
    largeFastaFile.saveFastaListToFile(mylistn,fastafiledic[key]+"nonredun")
    print(key, len(mylist), len(mylistn))

#20160311 new cufflink file
mylist = list(SeqIO.parse(open(fastafiledic["NC"]),"fasta"))
mylistn = largeFastaFile.fasta_within_seq_big_withError(mylist)
largeFastaFile.saveFastaListToFile(mylistn,fastafiledic["NC"]+"nonredun")
print("NC", len(mylist), len(mylistn))

#the computer shut down before I saw the results. Check the number of proteins in each file
for key in fastafiledic:
    mylist = list(SeqIO.parse(open(fastafiledic[key]),"fasta"))
    mylistu = largeFastaFile.fasta_uni_keepone(mylist)
    mylistn = list(SeqIO.parse(open(fastafiledic[key]+"nonredun"),"fasta"))
    print(key,len(mylist),len(mylistu),len(mylistn))

#combine TSn (TrinitySS non-redundant) and TUs (TrinityUnmap non-redundant) together to generate Tn (Trinity non-redundant)
mylist = list(SeqIO.parse(open(fastafiledic["TU"]+"nonredun"),"fasta"))
for ele in mylist: #rename element in TU, add letter "U" in front of id and description
    ele.id = "U"+ele.id
    ele.description = "U"+ele.description
print(len(mylist))
mylist = mylist + list(SeqIO.parse(open(fastafiledic["TS"]+"nonredun"),"fasta"))
print(len(mylist))
mylistn = largeFastaFile.fasta_within_seq_big_withError(mylist)
print(len(mylistn))
largeFastaFile.SaveFastaListToFile(mylistn,folder+"20160304TrinityNonredun.fa")
    
#combine BSn (BridgerSS non-redundant) and BUs (BridgerUnmap non-redundant) together to generate Bn (Bridger non-redundant)
mylist = list(SeqIO.parse(open(fastafiledic["BU"]+"nonredun"),"fasta"))
for ele in mylist: #rename element in TU, add letter "U" in front of id and description
    ele.id = "U"+ele.id
    ele.description = "U"+ele.description
print(len(mylist))
mylist = mylist + list(SeqIO.parse(open(fastafiledic["BS"]+"nonredun"),"fasta"))
print(len(mylist))
mylistn = largeFastaFile.fasta_within_seq_big_withError(mylist)
print(len(mylistn))
largeFastaFile.SaveFastaListToFile(mylistn,folder+"20160304BridgerNonredun.fa")



#20160307 file fastafileNRdic, dic store file location of non-redundant protein seqs.
fastafileNRdic ={}
fastafileNRdic["O"] = "E:\\store_for_D\\Insects Info\\Anopheles gambiae\\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.2.pepnonredun"
folder = "E:\\store_for_D\\Insects Info\\Anopheles gambiae\\Assemblies\\"
fastafileNRdic["C"] = folder + "AgCufflinks20160301.fa.transdecoder.pepnonredun"
fastafileNRdic["T"] = folder + "20160304TrinityNonredun.fa"
fastafileNRdic["B"] = folder + "20160304BridgerNonredun.fa"
fastafileNRdic["N"] = folder + "AgCufflinks20160310.fa.transdecoder.pepnonredun"

#open 4 fasta files to list, store in fastadic
fastadic ={}
for key in fastafileNRdic:
    fastadic[key] = largeFastaFile.open_fasta_to_list(fastafileNRdic[key])
    print(key,len(fastadic[key]))

import largeFastaFile as lf
import XCBlastpLike as bp
#open 4 fasta files to namelist and seqlist, store in namedic and seqdic
namedic ={}
seqdic ={}
for key in fastafileNRdic:
    namedic[key], seqdic[key] =lf.openFasta2NameSeqlist(fastafileNRdic[key])
    print(len(namedic[key]))

#generate dickmers for Cuff
dickmerC = bp.seqs2kmerdic(seqdic["C"])
#Count how many proteins in O the same as in C
kmerlen =5
same = 0
for seq in seqdic["O"]:
    targets = bp.seqFindtargets(seq,dickmerC)
    for seq2id, seq2_counts in targets:
        seqkmer = bp.seq2kmer(seq)
        if seq2_counts == len(seqkmer):
            if seq in seqdic["C"][seq2id]:
                same += 1
                break
print(len(seqdic["O"]), same)

Otuple = tuple(seqdic["O"])
Ctuple = tuple(seqdic["C"])
same =[]
samein = []
import time
time0 = time.time()
for num in range(len(seqdic["O"])):
#    if num%1000 ==0:
#            print (num, same, samein,time.time()-time0)
    for num2 in range(len(seqdic["C"])):
        if Otuple[num] == Ctuple[num2]:
            same.append((num,num2))
            samein.append((num,num2))
            break
        elif Otuple[num] in Ctuple[num2]:
            samein.append((num,num2))
            break
print(len(same), len(samein))

sameO =[i[0] for i in same]
for num in range(len(seqdic["O"])):
    if num not in sameO:
        print(namedic["O"][num])
        break





#build database use protein sequence of Ag Cufflinks, TrinitySS, BridgerSS, TrinityUnmap, BridgerUnmap, AgOGS
folderblastdb = "E:\\Lab\\fastaDB\\Ag\\"
import localBlast
import largeFastaFile

dcBuild = {}
dcBuild["O"] = folderblastdb+"OGS\\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"
dcBuild["B"] = folderblastdb+"BridgerNR\\20160304BridgerNonredun.fa"
dcBuild["T"] = folderblastdb+"TrinityNR\\20160304TrinityNonredun.fa"
dcBuild["C"] = folderblastdb +"CuffN\\AgCufflinks20160310.fa.transdecoder.pepnonredun"

localBlast.blastDBbuild(dcBuild["O"],dcBuild["O"])
localBlast.blastDBbuild(dcBuild["B"],dcBuild["B"])
localBlast.blastDBbuild(dcBuild["T"],dcBuild["T"])
localBlast.blastDBbuild(dcBuild["C"],dcBuild["C"])

from itertools import product
for keys in product(dcBuild,repeat = 2):#4files, blast to each other, including self blast
    print(keys)
    if keys != ("B","B") and keys != ("C","C")and keys != ("O","O"):
        localBlast.localblastp(dcBuild[keys[0]],dcBuild[keys[1]],folderblastdb+keys[0]+"2"+keys[1]+"blastp.txt",6,0.00001,4,"BLOSUM62",10)

#localBlast.localblastp(dcBuild["O"],dcBuild["C"],folderblastdb+"O2C.txt",6,0.00001,4,"BLOSUM62",10)
#localBlast.localblastp(dcBuild["O"],dcBuild["T"],folderblastdb+"O2T.txt",6,0.00001,4,"BLOSUM62",10)
#localBlast.localblastp(dcBuild["O"],dcBuild["B"],folderblastdb+"O2B.txt",6,0.00001,4,"BLOSUM62",10)
#localBlast.localblastp(dcBuild["C"],dcBuild["T"],folderblastdb+"O2T.txt",6,0.00001,4,"BLOSUM62",10)



#20160314 store files in dcBuild in dictionary
dcFa ={}
for key in dcBuild:
    dcFa[key] = largeFastaFile.open_fasta_to_list(dcBuild[key])
    print(key, len(dcFa[key]))

dcFalen ={}
for key in dcFa:
    dcFalen[key] =[]
    for ele in dcFa[key]:
        seq = str(ele.seq)
        if seq[-1] == "*":
            seq = seq[:-1]
        dcFalen[key].append(len(seq))
    print(key, len(dcFalen[key]),lenCount(dcFalen[key]),sum(dcFalen[key])/len(dcFalen[key]))
    

dcdcFa = {} # open each file to dictionary
for key in dcBuild:
    dcdcFa[key] = largeFastaFile.open_fasta_to_dic(dcBuild[key])
    print(key, len(dcdcFa[key]))


#compare based on blast6
key = "O2C"
key = "O2T"
key = "O2B"
key = "C2T"
key = "C2O"


folderBlast6 = "E:\\store_for_D\\Insects Info\\Anopheles gambiae\\Assemblies\\blast2eachother\\"
import XCBlast
mldic = {}
dcdcFa = {}

dcdcFa[key[0]] = largeFastaFile.open_fasta_to_dic(dcBuild[key[0]])
dcdcFa[key[2]] = largeFastaFile.open_fasta_to_dic(dcBuild[key[2]])
mldic[key] = XCBlast.mcotMLwithblast6(dcdcFa[key[0]],dcdcFa[key[2]],folderBlast6+key+"blastp.txt")
fout = open(folderBlast6+key+"ml.txt")
for ele in mldic[key]:
    fout.write(ele +"\t" + str(mldic[key][ele])+"\n")
fout.close()


dcFa1 = dcdcFa["O"]
dcFa2 = dcdcFa["C"]
file_blast6 = folderBlast6+"O2C"+"blastp.txt"


#after many test, doing local alignment is too slow, which may takes days to finish. I dicided to go back to the original MCOT method
myTrinity = largeFastaFile.open_fasta_to_list(fastafiledic["TS"].split(".transdecoder")[0]) +\
largeFastaFile.open_fasta_to_list(fastafiledic["TU"].split(".transdecoder")[0])
fout = open(folder+"20160314TrinityAll.fa","w")
for num in range(len(myTrinity)):
    fout.write(">T"+str(num+1)+"\n"+str(myTrinity[num].seq)+"\n")
fout.close()

myBridger = largeFastaFile.open_fasta_to_list(fastafiledic["BS"].split(".transdecoder")[0]) +\
largeFastaFile.open_fasta_to_list(fastafiledic["BU"].split(".transdecoder")[0])
fout = open(folder+"20160314BridgerAll.fa","w")
for num in range(len(myBridger)):
    fout.write(">B"+str(num+1)+"\n"+str(myBridger[num].seq)+"\n")
fout.close()

##20160406 New Trinity Assembly
def TrinityNewGetUniqueSeq():
    import largeFastaFile
    filename = r"F:\Insects\Anopheles_gambiae\TrinityNew\TrinityConbined\20160406TrinityCombined98.fa.transdecoder.pep"
    from Bio import SeqIO
    ls_t = list(SeqIO.parse(open(filename),"fasta"))
    ls_tNonwithin = largeFastaFile.fasta_within_seq_big_withError(ls_t)
    fout = open(filename+"NoWithin","w")
    for seq in ls_tNonwithin:
        fout.write(">"+seq.description+"\n"+str(seq.seq)+"\n")
    fout.close()



def Ag20160412MCOT2GetDNA():
    folder = "F:\\Insects\\Anopheles_gambiae\\MCOT2\\"
    f_pep = folder+"MCDpeptides.txt"
    from Bio import SeqIO
    ls_fa = list(SeqIO.parse(open(f_pep),"fasta"))
    ls_oriName =[]
    for ele in ls_fa:
        _oriID = ele.description.split("\t")[2]
        ls_oriName.append((ele.description, _oriID))
    dc_allCDS ={}
    f_cdsT = "F:\\Insects\\Anopheles_gambiae\\TrinityNew\\TrinityConbined\\20160406TrinityCombined98.fa.transdecoder.cds"
    folder2 = "F:\\Insects\\Anopheles_gambiae\\Assemblies\\"
    f_cdsC = folder2+"AgCufflinks20160310.fa.transdecoder.cds"
    f_cdsB = folder2 + "20160314BridgerAll.fa.transdecoder.cds"
    f_cdsO = "F:\\Insects\\Anopheles_gambiae\\Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.3Rename2pep.fa"
    dc_allCDS.update(SeqIO.to_dict(SeqIO.parse(open(f_cdsT),"fasta")))
    dc_allCDS.update(SeqIO.to_dict(SeqIO.parse(open(f_cdsC),"fasta")))
    dc_allCDS.update(SeqIO.to_dict(SeqIO.parse(open(f_cdsB),"fasta")))
    dc_allCDS.update(SeqIO.to_dict(SeqIO.parse(open(f_cdsO),"fasta")))
    fo_cds = open(folder +"MCDcds.txt","w")
    for ele in ls_oriName:
        _cds = str(dc_allCDS[ele[1]].seq)
        fo_cds.write(">"+ele[0]+"\n"+_cds+"\n")
    fo_cds.close()

def AgGenomeDescription20170228():
    '''
    get some information about the genome of Ag, including length, soft-masked region and NNNs
    '''
    f_ch = r"D:\Insects\Anopheles_gambiae\Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
    from Bio import SeqIO
    ls_ch = list(SeqIO.parse(open(f_ch),'fasta'))
    print('chromosome number is', len(ls_ch))
    
    for _ch in ls_ch:
        print(_ch.id,len(_ch.seq))
    
    bases = set()
    for _ch in ls_ch:
        bases.update(set(str(_ch.seq)))
    print('bases:\n', bases)
    
    from collections import Counter
    baseCount = Counter(''.join(str(_ch.seq) for _ch in ls_ch))
    for _base in baseCount:
        print(_base,'\t',baseCount[_base])
    
        
    #check NNNs
    genomeseq = ''.join(str(_ch.seq) for _ch in ls_ch)
    genomeseq = genomeseq.upper()
    import re
    NNNs = re.split('[ATCG]+',genomeseq)
    NNNs = NNNs[1:-1]
    NNNs = [len(_nnn) for _nnn in NNNs]
    NNNcount = Counter(NNNs)
    fout = open('list.txt','w')
    for _e in NNNcount:
        fout.write(str(_e)+'\t'+str(NNNcount[_e])+'\n')
    fout.close()
    del genomeseq
    
    #save NNNcount for figure. add distubtion of 0.3
    fout = open('list.txt','w')
    lsNNNcount = list(NNNcount.items())
    lsNNNcount = sorted(lsNNNcount,key = lambda x:(x[1],x[0]))
    counts = [_e/100 for _e in range(0,50,2)]
    
    for n,ele in enumerate(lsNNNcount):
        _klen = ele[0]
        _kcount = ele[1]
        _kcountAjust = _kcount + counts[n%len(counts)]
        fout.write('%d\t%f\n'%(_klen,_kcountAjust))
    fout.close()
    
    xmin = 19
    xmax = 45000
    fractions = 300
    import math
    xlogfold = math.log10(xmax/xmin)
    xlims = [20*10**(n/fractions * xlogfold) for n in range(fractions)]
    
    dcNNNcount = {}
    for ele in lsNNNcount:
        if ele[1] not in dcNNNcount:
            dcNNNcount[ele[1]] = []
        dcNNNcount[ele[1]].append(ele[0])
    def lsValueSplitter(ls,lims):
        '''
        given a list, return a list of list of values in the range provided by lims
        ls = [0,1,2,3,4,5,6,7]
        lims = [2,5,8]
        return [[0,1],[2,3,4],[5,6,7],[]]
        '''
        ls = sorted(ls)
        lims = [float('-inf')]+  lims+[ float('+inf')]
        lsls = [[ele for ele in ls if ele>=lims[n] and ele <lims[n+1]] for n in range(len(lims)-1)]
        return lsls
    for ele in dcNNNcount:
        dcNNNcount[ele] = lsValueSplitter(dcNNNcount[ele],xlims)
    
    pertb = 0.01
    dcX = {}
    dcY = {}
    for ele in dcNNNcount:
        dcX[ele] = []
        dcY[ele] = []
        for _ls in dcNNNcount[ele]:
            for _n, _v in enumerate(_ls):
                dcX[ele].append(_v)
                if _n %2 == 0:
                    dcY[ele].append(ele + pertb * _n)
                else:
                    dcY[ele].append(ele - pertb * _n)

    

    import numpy as np
    import matplotlib.pyplot as plt
    xs = dcX[1]
    ys = dcY[1]
#    plt.plot(xs,ys, marker = 'o',linestyle = '',markeredgewidth = 0, markersize = 1)
#    plt.plot(dcX[2],dcY[2], marker = 'o',linestyle = '',markeredgewidth = 0, markersize = 1)
    for n in range(1,11):
        plt.plot(dcX[n],dcY[n], marker = 'o',linestyle = '',markeredgewidth = 0, markersize = 1)
    plt.ylim(-2,10.5)
    plt.xscale('log')
    plt.savefig('test.pdf')
    plt.show()
    
    #plot NNNs directly
    NNNs = sorted(NNNs)
    NNNsFig = [_e for _e in NNNs if _e != 20 and _e != 1000]
    NNNsSplit = lsValueSplitter(NNNsFig,xlims)
    NNNsY = []
    pertb = 0.02
    for _ls in NNNsSplit:
        for _n, _v in enumerate(_ls):
            NNNsY.append(1 + pertb * _n)

    plt.figure(figsize = (8,2))
    plt.plot(NNNsFig,NNNsY, marker = 's',linestyle = '',markeredgewidth = 0, markersize = 1)
    plt.xscale('log')
    
    plt.savefig('NNNs.pdf')
    plt.show()
    
    
    #plot NNNs in chromosomes
    f_ch = r"D:\Insects\Anopheles_gambiae\Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
    from Bio import SeqIO
    ls_ch = list(SeqIO.parse(open(f_ch),'fasta'))
    for _seq in ls_ch:
        _seq.seq = str(_seq.seq).upper()
    print('chromosome number is', len(ls_ch))
    
    dcNNNs = {}#chromosome id as key, a list of (NNN starting site, length)
    import re
    for _seq in ls_ch:
        dcNNNs[_seq.id] = []
        for _m in re.finditer('N+',_seq.seq):
            _s = _m.start()
            _l = _m.end() - _m.start()
            dcNNNs[_seq.id].append((_s,_l))
    
    dcChlen = {}
    for _seq in ls_ch:
        dcChlen[_seq.id] = len(_seq.seq)
    
    keyOrder = ['2L','2R','3L','3R','X','UNKN','Y_unplaced']
    NNNinfo = dcNNNs['2L']#combine all chromosome info to one list
    move = dcChlen['2L']
    for key in keyOrder[1:]:
        _NNN = dcNNNs[key].copy()
        _NNN = [(_e[0] + move, _e[1]) for _e in _NNN]
        move += dcChlen[key]
        NNNinfo = NNNinfo + _NNN
    print(len(NNNinfo))
    
    import matplotlib.pyplot as plt
    plt.rcParams["figure.figsize"] = [11,8.5]
    
    #plot chromosome.
    ch_xloc = [0]
    for key in keyOrder:
        ch_xloc.append(ch_xloc[-1] + dcChlen[key])
    print(ch_xloc)
    colors = ['b','g','r','c','m','y','k']
    for _i in range(len(ch_xloc) -1):
        plt.plot(ch_xloc[_i:_i+2],[0,0],color = colors[_i])
    plt.savefig('ChNNNs.pdf')
    #the figure is too compact, decide to plot each chromosomes individually
    
    #plot chromosome individually
    plt.rcParams["figure.figsize"] = [11,5]
    span = 200000
    pertb = 0.02
    basey = pertb
    import math
    
    def getys(xs,xslims,basey=0,pertb = 0.01):
        '''
        given xs, xslims, basey value, return a list of ys for plot
        '''
        xxs = lsValueSplitter(xs,xslims)
        ys = []
        for _ls in xxs:
            for _n, _v in enumerate(_ls):
                ys.append(basey + pertb * _n)
        return ys
    for _ch in keyOrder:
        basey = basey + 1
        NNNs = dcNNNs[_ch]
        plt.plot([0,dcChlen[_ch]],[basey-pertb,basey-pertb],color = 'k')
        xslims = list(range(0,dcChlen[_ch],span))
    #    xs = [_e[0] for _e in NNNs]
    #    ys = getys(xs,xslims,basey)
    #    sizes = [math.sqrt(_e[1])/math.sqrt(20) for _e in NNNs]
    #    colors = []
    #    for ele in NNNs:
    #        if ele[1] == 20:
    #            colors.append((0,1.0,1.0,0.5))
    #        elif ele[1] == 1000:
    #            colors.append((1.0,1.0,0,0.5))
    #        else:
    #            colors.append((1.0,0,1.0,0.2))
    #    plt.scatter(xs,ys,color = colors, marker = 'o', s = sizes,linewidth = 0)
        dcxs = {20:[],1000:[],'other':[]}
        dcsizes = {20:[],1000:[],'other':[]}
        for ele in NNNs:
            if ele[1] == 20:
                dcxs[20].append(ele[0])
                dcsizes[20].append(1)
            elif ele[1] == 1000:
                dcxs[1000].append(ele[0])
                dcsizes[1000].append(1)
            else:
                dcxs['other'].append(ele[0])
                dcsizes['other'].append(math.sqrt(ele[1])/math.sqrt(20))
        dcys = {}
        for _key in dcxs:
            if _key == 'other':
                xslims = list(range(0,dcChlen[_ch],span*2))
                dcys[_key] = getys(dcxs[_key],xslims,basey-pertb,-pertb)
            else:
                dcys[_key] = getys(dcxs[_key],xslims,basey,pertb)
        for _key in dcxs:
            if _key == 20:
                colors = (0,0,1.0,0.5)
            elif _key == 1000:
                colors = (0.25,0.7,0.25,0.2)
            else:
                colors = (1.0,0,1.0,0.2)
            plt.scatter(dcxs[_key], dcys[_key],color = colors, marker = 'o', s = dcsizes[_key],linewidth = 0)
    plt.savefig('2Ltest.pdf')
    

def AgMCOT2CompareWithAgP44_20170329():
    '''
    compare proteins of AgMCOT2 with AgP44
    '''
    f_mcot = r"D:\Insects\Anopheles_gambiae\MCOT2\MCDpeptides.txt"
    f_45 = r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.5.fa"
    folder = "D:\\Insects\\Anopheles_gambiae\\MCOT2\\compareWIthOGS4.4\\"
    f_blast = r"D:\Insects\Anopheles_gambiae\MCOT2\compareWIthOGS4.4\2017AgMCOT2OGS4.4Tab.txt"
    
    runfile('D:/P/3Language/Xiaolong/python/XCProject/MCuNovoGeneSelector/MCuNovoGeneSelectorMain.py', wdir=folder)
    import os
    os.getcwd()
    
    dc = file_EditedBlast6_to_best_match(f_blast6=f_blast,f_query=f_mcot,f_subject=f_45,outputfile='20170329AgMCOT2OGS4.4.txt',min_identity=95,returndic = True)
    
    import largeFastaFile
    dc_mcot = largeFastaFile.open_fasta_to_dic(f_mcot)
    
    lsMcotonly=[]
    for key in dc_mcot:
        if key not in dc:
            lsMcotonly.append(key)
    print(len(lsMcotonly))
    
    open('list.txt','w').write('\n'.join(lsMcotonly))
    
    st_geneshare = set([ele.split('.')[0] for ele in dc])
    st_mcotonly = set([ele.split('.')[0] for ele in lsMcotonly])
    
    ls_mcotonlyGene = [ele for ele in st_mcotonly if ele not in st_geneshare]
    print(len(st_geneshare), len(st_mcotonly), len(ls_mcotonlyGene))
    
    fout = open('list.txt','w')
    for ele in lsMcotonly:
        if ele.split('.')[0] not in set(ls_mcotonlyGene):
            fout.write(ele+'\t'+'N')
        else:
            fout.write(ele+'\t'+'Y')
        fout.write('\t'+dc_mcot[ele].description+'\n')
    fout.close()
    
    fout = open('20170329MCOTuniqueGene.txt','w')
    for ele in dc_mcot:
        ele = dc_mcot[ele]
        if ele.id.split('.')[0] in ls_mcotonlyGene:
            fout.write('>'+ele.description+'\n'+str(ele.seq)+'\n')
    fout.close()
    
    dc_mcotonlyGeneLength = {}
    for ele in dc_mcot:
        ele = dc_mcot[ele]
        if ele.id.split('.')[0] in ls_mcotonlyGene:
            dc_mcotonlyGeneLength[ele.id] = int(ele.description.split()[1][4:])
    
    for n in [100,200,300,400]:
        print('number longer than %d:'%n,len([ele for ele in dc_mcotonlyGeneLength.values() if ele >n]))

def AgMCOT2OGS45_20170407():
    '''
    compare sequences of AgMCOT 2.0 with AgamP4.5. With the new function by muscle. sequenceComarison.py result
    '''
    import sys
    sys.path.append(r'D:\P\3Language\Xiaolong\python\XCProject\MCuNovoGeneSelector')
    import sequenceComparison
    f_fasta1 = r"D:\Insects\Anopheles_gambiae\MCOT2\MCDpeptides.txt"
    f_fasta2 = r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.5.fa"
    seqs1 = sequenceComparison.openfile2lsFasta(f_fasta1)
    seqs2 = sequenceComparison.openfile2lsFasta(f_fasta2)
    # readin sequences
    #readin comparison result
    import pandas as pd
    df = pd.read_csv(r"D:\Insects\Anopheles_gambiae\MCOT2\compareWIthOGS4.4\20170406MCOT2.0ToOGS4.5_ID_compare.txt",delimiter = '\t')
    df = df.sort(columns = ['seq1_id','matched_length'], ascending = [True, False])
    #keep the best match for each query
    df = df.drop_duplicates(subset=['seq1_id'])
    dfqs = []
    for n in range(df.shape[0]):
        _ls = list(df.iloc[n,:])
        dfqs.append((seqs1[_ls[0]].id, len(seqs1[_ls[0]].seq), seqs2[_ls[1]].id, len(seqs2[_ls[1]].seq),_ls[2]))
    dfqs = pd.DataFrame(dfqs, columns=['query','ql','subject','sl','ml'])
    
    #append those with no match in AgamP4.5
    querys = set(dfqs['query'])
    lsquerys = []
    for seq in seqs1:
        if seq.id not in querys:
            lsquerys.append((seq.id,len(seq.seq),'na',0,0))
    lsquerys = pd.DataFrame(lsquerys,columns=['query','ql','subject','sl','ml'])
    
    dfqs = dfqs.append(lsquerys)
    
    def get_match_comment(s):
        if s[2] == 'na':
            return 'W'
        if s[1] == s[3] and s[1] == s[4]:
            return 'P'
        if s[4] * s[4] /(s[1] * s[3]) > 0.95:
            return 'N'
        if s[4] / (0.7 * s[1]) + s[4]/200 >1:
            return 'O'
        return 'B'
    #add comment column
    match_comment = dfqs.apply(get_match_comment,axis=1)
    dfqs['match_comment'] = match_comment
    dfqs.to_csv(r"D:\Insects\Anopheles_gambiae\MCOT2\compareWIthOGS4.4\20170406MCOT2.0ToOGS4.5_ID_compareResult.csv")


def AgMCuNovo_Trinity20170508():
    '''
    Trinity assemble of RNA-seq reads in one run. 
    '''
    folder = 'F:\\InsectBigData\\AnophelesGambiae\\RNA-seq\\'
    import os
    files = os.listdir(folder)
    files = [e.split('.')[0] for e in files]
    
    #decompress files
    for file in files:
        fout = open('20170508_1_decompress'+file+'.txt', 'w')
        fout.write('''
#!/bin/bash
#PBS -q express
#PBS -N {0}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe
cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq
gzip -d {0}.fa.gz
                   '''.format(file))
        fout.close()
    
    #reads2Transcripts
    for file in files:
        fout = open('20170509_2reads2transcriptsAg_'+file+'.txt','w')
        fout.write('''
#!/bin/bash
#PBS -q batch
#PBS -N {0}
#PBS -l nodes=1:ppn=12
#PBS -l walltime=24:00:00
#PBS -m abe -M atps@outlook.com
#PBS -j oe

module reset
module load python
module load perl
module load gcc-4.9.2
module load jdk/1.8.0_45
module load bowtie2

cd /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/RNA-seq
/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Chrysalis/ReadsToTranscripts -i {0}.fa -f /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/2017Trinity/bundled_iworm_contigs.fasta -o ../reads/{0}.fa  -t 12 -max_mem_reads 50000000 
wc -l {0}.fa >{0}.count
/bin/sort -T . -S 30G -k 1,1n ../reads/{0}.fa > ../reads/{0}.fa.sort
wc -l ../reads/{0}.fa.sort >../reads/{0}.fa.count
rm ../reads/{0}.fa
                   '''.format(file))
        fout.close()
    
    #generate recursive

    
    seqs = range(0,182634)
    fout = open('recursive_trinity.cmds','w')
    for seqn in seqs:
        foldern = seqn//1000
        fout.write('/panfs/panfs.cluster/home/ks2073/p/2016/trinityrnaseq-Trinity-v2.3.2/Trinity --single /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/seqs/{1}/{0} --output /panfs/panfs.cluster/scratch/ks2073/Transcriptome/AnophelesGambiae/seqs/{1}/trinityout{0} --CPU 1 --max_memory 2G --seqType fa --trinity_complete --full_cleanup --KMER_SIZE 27\n'.format(seqn,foldern))
    fout.close()
    
def renameTrinityResult20170515():
    f = r"D:\Insects\Anopheles_gambiae\2017Trinity\list.txt"
    from Bio import SeqIO
    ls = list(SeqIO.parse(open(f),'fasta'))
    gene_id = 1
    i = 0
    j = 1
    fout = open(r"D:\Insects\Anopheles_gambiae\2017Trinity\2017AgTrinity_singleRun.txt",'w')
    while j<len(ls):
        seq1 = ls[i]
        seq2 = ls[j]
        gene1,iso1 = seq1.id.rsplit('_',1)
        gene2,iso2 = seq2.id.rsplit('_',1)
        if gene1 == gene2 and int(iso2[1:]) - int(iso1[1:]) == 1:
            fout.write('>AgTrinity%06d'%gene_id+'.'+iso1[1:]+'\n'+str(seq1.seq) +'\n')
            i+=1
            j+=1
        else:
            fout.write('>AgTrinity%06d'%gene_id+'.'+iso1[1:]+'\n'+str(seq1.seq) +'\n')
            i+=1
            j+=1
            gene_id+=1
    fout.write('>AgTrinity%06d'%gene_id+'.'+iso2[1:]+'\n'+str(seq2.seq) +'\n')
    fout.close()
    
    