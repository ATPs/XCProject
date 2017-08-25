# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:28:14 2015

@author: k
"""

import sys
sys.path.append("D:\\P\\3Language\\Xiaolong\\python\\XCProject\\Shared\\")

import File

def clearfile(filename):
    """
    given the filename of short reads blast, output the filtered result, like
    33-178252_68	gi|685047575|emb|LN596087.1|	95.65	69	1	2	1	68	6210	6143	3e-21	  110
    33-178252_68	gi|685046121|emb|LN594646.1|	95.65	69	1	2	1	68	156230	156297	3e-21	  110
    33-178252_68	gi|685043277|emb|LN591802.1|	95.65	69	1	2	1	68	63566	63633	3e-21	  110
    33-178252_68	gi|685042530|emb|LN591055.1|	95.65	69	1	2	1	68	361333	361266	3e-21	  110
    33-178252_68	gi|685042440|emb|LN590965.1|	95.65	69	1	2	1	68	76679	76612	3e-21	  110
    33-178252_68	gi|685042206|emb|LN590731.1|	95.65	69	1	2	1	68	347611	347678	3e-21	  110
    33-178252_68	gi|685042180|emb|LN590705.1|	95.65	69	1	2	1	68	14955371	14955304	3e-21	  110
    33-178252_68	gi|685042180|emb|LN590705.1|	94.20	69	2	2	1	68	10473029	10472962	1e-19	  104
    33-178252_68	gi|685042180|emb|LN590705.1|	94.20	69	2	2	1	68	13198412	13198345	1e-19	  104
    33-178252_68	gi|685042180|emb|LN590705.1|	94.20	69	2	2	1	68	13923842	13923909	1e-19	  104
    33-178252_68	gi|685042180|emb|LN590705.1|	94.20	69	1	3	1	68	14508261	14508327	5e-19	  102
    33-178252_68	gi|685042180|emb|LN590705.1|	93.85	65	4	0	4	68	2155408	2155472	6e-18	99.0
    33-178252_68	gi|685042180|emb|LN590705.1|	93.75	64	2	2	6	68	2425115	2425177	8e-17	95.3
    33-178252_68	gi|685042175|emb|LN590700.1|	95.65	69	1	2	1	68	14033940	14033873	3e-21	  110
    33-178252_68	gi|685042175|emb|LN590700.1|	94.20	69	2	2	1	68	15778539	15778606	1e-19	  104
    33-178252_68	gi|685042146|emb|LN590671.1|	95.65	69	1	2	1	68	3854573	3854640	3e-21	  110
    33-178252_68	gi|916497858|ref|XM_013494305.1|	94.20	69	2	2	1	68	412	345	1e-19	  104
    filter, identity>95%, matched length > 90% of query length
    only keep the best bitscores, in this case, the first 6 lines
    for the id, only keep the gi numbers. output to the filename+"Out"
    """
    mylist = File.open2list(filename)
    newlist = []
    for ele in mylist:
        eles = ele.split()
        eles[1] = eles[1].split("|")[1]
        if float(eles[2])> 95 and float(eles[3])/float(eles[0].split('_')[1]) > 0.9:
            eles.append(0)
            newlist.append(eles)
    print("original there are ", len(mylist), "lines in the file, ", len(newlist), "lines left after filtering identity and coverage")
    new2list=[]
    newlist[0][-1] = 1
    for n in range(len(newlist)-1):
        if newlist[n][0] == newlist[n+1][0]:
            if newlist[n][11] == newlist[n+1][11]:
                newlist[n+1][12] = newlist[n][12]
        else:
            newlist[n+1][12] = 1
    for ele in newlist:
        if ele[12] == 1:
            new2list.append(ele)
    print("after filter best hits, total lines is", len(new2list))
    fo = open(filename+"filtout","w")
    for ele in new2list:
        elements=""
        for elem in ele[:-1]:
            elements += elem +"\t"
        elements += "\n"
        fo.write(elements)
    fo.close()
    return "Done!"

import os
folder = 'E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\'
files = os.listdir(folder)
for file in files:
    clearfile(folder+file)


folder2 = 'E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\selected\\'
files = os.listdir(folder2)
mygis = set()
for file in files:
    mylist = File.open2list(folder2+file)
    for ele in mylist:
        gi = ele.split()[1]
        mygis.add(gi)
print("total GI number is: ", len(mygis))

fo = open("gi.txt","w")
for ele in mygis:
    fo.write(ele+"\n")
fo.close()


fo = open("E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\SubjectGIs\\gi.txt","r")
mygis = []
for ele in fo.readlines():
    mygis.append(int(ele[:-1]))
fo.close()


mygis =set(mygis)
i = 0
fo2 = open("F:\\gi_taxid_nucl.dmp","r")
fo3 = open("gitaxons.txt","w")
taxon = fo2.readline()
while taxon:
    gi = int(taxon.split()[0])
    if gi in mygis:
        fo3.write(taxon)
    i += 1
    if i % 1000000 ==0:
        print(i, "searched in mygis")
    taxon = fo2.readline()
fo3.close()        
fo2.close()




folder = 'E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\SubjectGIs\\'
file ="noTaxonYet.txt"
mygis = File.open2list(folder+file)
url1 = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=genbank&sort=&val='
url2 = '&from=begin&to=end&extrafeat=976&maxplex=3'
for ele in mygis:
    page = urllib.request.urlopen(url1+ele[:-1]+url2).read()
    fo = open(ele[:-1],"wb")
    fo.write(page)
    fo.close()



folder = 'E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\SubjectGIs\\temp\\'
import os
mylist = os.listdir(folder)
taxons =open("taxons.txt","w")
for file in mylist:
    text = open(folder+file,"r").read()
    n1 = text.find("taxon:")
    if n1>0:
        text2 = text[n1:].split("\"",2)[0]
        text3 = text2.split(":")[1]
        taxons.write(file+"\t"+text3+"\n")
    else:
        print("taxon for this one is not found", file)
taxons.close()


file ="E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\SubjectGIs\\AllGiTaxons.txt"
fo = open(file,"r")
gitaxondic = {}
line = fo.readline()
while line:
    gitaxondic[int(line.split()[0])] = int(line.split()[1])
    line = fo.readline()
fo.close()

import os
folder = "E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\selected\\"
myfiles = os.listdir(folder)
folder2 = "E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\selectedWithTaxon\\"
def addingtaxon(filename, gitaxondic, folder = ""):
    """
    given the selected file of blastn. Add the taxon id to the end of the file
    """
    mylist = open(filename,"r").readlines()
    fout = open(folder+filename.split("\\")[-1]+"taxon", "w")
    for ele in mylist:
        ele = ele[:-1] + str(gitaxondic[int(ele.split()[1])]) +"\n"
        fout.write(ele)
    fout.close()
    print(filename.split("\\")[-1], "done")
for file in myfiles:
    addingtaxon(folder+file, gitaxondic, folder2)
    




folder = "E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\selectedWithTaxon\\"
file = "E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\SubjectGIs\\TaxonIDs.txt"
taxdic ={}
mytaxon = open(file,"r").read().split()
for ele in mytaxon:
    taxdic[int(ele)] =0

import os
myfiles = os.listdir(folder)
for file in myfiles:
    for ele in open(folder+file,"r").readlines():
        taxon = int(ele.split()[-1])
        count = int(ele.split()[0].split('-')[1].split("_")[0])
        taxdic[taxon] += count
    print(file, "done")

sum(taxdic.values())

fo = open("E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\SubjectGIs\\TaxonIDsAndReadsCount.txt", "w")
for ele in taxdic:
    fo.write(str(ele)+"\t"+str(taxdic[ele])+"\n")
fo.close()






folder = "E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\selectedWithTaxon\\"
myfiles = os.listdir(folder)
finalout = open("E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\Blastn1vs1best.txt","w")
numss = 0
for file in myfiles:
    tempdic ={}
    templist = open(folder+file,"r").readlines()
    querynum = set()
    for ele in templist:
        tempdic[ele.split()[0]] =[]
        querynum.add(ele.split()[0])
    for ele in templist:
        tempdic[ele.split()[0]].append(ele)
    for ele in tempdic:
        ttlist = tempdic[ele]
        if ttlist == 1:
            finalout.write(ttlist[0])
        else:
            temptaxvalues = []
            for num in range(len(ttlist)):
                temptaxvalues.append(taxdic[int(ttlist[num].split()[-1])])
            maxvalue = max(temptaxvalues)
            for num in range(len(temptaxvalues)):
                if temptaxvalues[num] == maxvalue:
                    finalout.write(ttlist[num])
                    break
    print(file, "done", len(querynum))
    numss += len(querynum)
finalout.close()






folder = "E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\selectedWithTaxon\\"
myfiles = os.listdir(folder)
myset = set()
for file in myfiles:
    templist = open(folder+file,"r").readlines()
    for ele in templist:
        myset.add(int(ele.split()[0].split('-')[0]))
print(len(myset))


acount = 0
for ele in open("E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\Blastn1vs1best.txt","r").readlines():
    count += int(ele.split()[0].split('-')[1].split('_')[0])
print(acount)


mylist = open("E:\\Study\\store_for_D\\Transcriptome\\unMappedReads\\blastn\\Blastn1vs1best.txt","r").readlines()
taxons = set()
for ele in mylist:
    taxons.add(int(ele.split()[-1]))
print(len(taxons))

fo = open("finaltaxons.txt","w")
for ele in taxons:
    fo.write(str(ele)+'\n')
fo.close()

mydic={}
for ele in taxons:
    mydic[ele] = 0
for ele in mylist:
    elet = int(ele.split()[-1])
    elec = int(ele.split('-')[1].split('_')[0])
    mydic[elet] += elec
fo = open("finaltaxons.txt","w")
for ele in mydic:
    fo.write(str(ele)+'\t'+str(mydic[ele])+'\n')
fo.close()






