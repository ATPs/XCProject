# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 11:26:58 2016

@author: k
"""

import sqlite3

folder = "E:\\store_for_D\\KOBAS2\\DB_using_20160104\\DB\\"
db_ko = "ko.db"
db_bm = "bmor.db"
db_tc = "tca.db"
db_hs = "hsa.db"
db_dm = "dme.db"
db_ag = "aga.db"

conn_ko = sqlite3.connect(folder+db_ko)
conn_bm = sqlite3.connect(folder+db_bm)
conn_tc = sqlite3.connect(folder+db_tc)
conn_hs = sqlite3.connect(folder+db_hs)
conn_dm = sqlite3.connect(folder+db_dm)
conn_ag = sqlite3.connect(folder+db_ag)

c_ko = conn_ko.cursor()
c_bm = conn_bm.cursor()
c_tc = conn_tc.cursor()
c_hs = conn_hs.cursor()
c_dm = conn_dm.cursor()
c_ag = conn_ag.cursor()

conn_use = sqlite3.connect(folder+"20160106using.db")
c_use = conn_use.cursor()



def selectSave(commandline, mycursor, filename = "result.txt"):
    """
    provide a command, save the output to a filename provided. 
    """
    fo = open(filename,"w")
    tempresult = mycursor.execute(commandline)
    mylist = tempresult.fetchall()
    for ele1 in mylist:
        for ele2 in ele1:
            fo.write(ele2+"\t")
        fo.write("\n")
    fo.close()


selectSave("SELECT * from KoGenes WHERE gid LIKE 'hsa:%' ORDER BY koid ASC",c_ko,"ko_hsa.txt")
selectSave("SELECT * from KoGenes WHERE gid LIKE 'bmor:%' ORDER BY koid ASC",c_ko,"ko_bmor.txt")
selectSave("SELECT * from KoGenes WHERE gid LIKE 'tca:%' ORDER BY koid ASC",c_ko,"ko_tca.txt")
selectSave("SELECT * from KoGenes WHERE gid LIKE 'dme:%' ORDER BY koid ASC",c_ko,"ko_dme.txt")
selectSave("SELECT * from KoGenes WHERE gid LIKE 'aga:%' ORDER BY koid ASC",c_ko,"ko_aga.txt")

temp_pathway=c_ko.execute("SELECT DISTINCT pid FROM KoPathways")
pathways = temp_pathway.fetchall()



def generateMapFiles():
    """
    given a KEGG map ko, like ko00010, which represent map00010, glycolysis/gluconeogenesis, return a text file, with name of 
    the ko (ko00010_glycolisis_gluconeogenesis.txt)
    given the Pathway map, return 382 files for each pathway in the file.
    if no manduca genes in the map, do not generate the file
    """
    folder = "E:\\store_for_D\\KOBAS2\\DB_using_20160104\\"
    pathways = open(folder+"Pathways.txt").readlines()
    kopathways = open(folder+"KoPathway.txt").readlines()
    kopathwaypair = []
    koaga = open(folder+"ko_aga.txt").read().split()
    kobmor = open(folder+"ko_bmor.txt").read().split()
    kodme = open(folder+"ko_dme.txt").read().split()
    kohsa = open(folder+"ko_hsa.txt").read().split()
    koMsex = open(folder+"ko_Msex.txt").read().split()
    kotca = open(folder+"ko_tca.txt").read().split()
    kos = open(folder+"Kos.txt").readlines()
    kosdic ={}
    for _ in kos[1:]:
        kosdic[_.split("\t")[0]] = _.split("\t")[1][:-1]
    koec = open(folder+"ko_ec.txt").readlines()
    koecdic ={}
    for _ in koec:
        koecdic[_.split("\t")[0]] = _.split("\t")[1][:-1]
    for _ in kos[1:]:
        if _.split("\t")[0] not in koecdic:
            koecdic[_.split("\t")[0]] = ""
    for _ in kopathways[1:]:
        kopathwaypair.append(tuple(_.split())) #tuple(koid, pid)
    for pathway in pathways[1:]:
        kop,kopname = pathway.split()
        pathwaykolist = []
        for _ in kopathwaypair:
            if _[1] == kop:
                pathwaykolist.append(_[0])
        def pathwaykolistcheck(pathwaykolist,kospeices):
            pathwaykolist_species = []
            for _ in pathwaykolist:
                if _ in kospeices:
                    pathwaykolist_species.append(1)
                else:
                    pathwaykolist_species.append(0)
            return pathwaykolist_species
        pathwaykolist_aga = pathwaykolistcheck(pathwaykolist,koaga)
        pathwaykolist_bmor = pathwaykolistcheck(pathwaykolist,kobmor)
        pathwaykolist_dme = pathwaykolistcheck(pathwaykolist,kodme)
        pathwaykolist_hsa = pathwaykolistcheck(pathwaykolist,kohsa)
        pathwaykolist_Msex = pathwaykolistcheck(pathwaykolist,koMsex)
        pathwaykolist_tca = pathwaykolistcheck(pathwaykolist,kotca)
        if sum(pathwaykolist_Msex) == 0:
            fo = open(folder+"test\\"+kop+"_"+kopname+".txt","w")
            fo.write("ko\tname\tec\taga\tbmor\tdme\tMsex\ttca\t\n")
            for _ in range(len(pathwaykolist)):
                if sum([pathwaykolist_aga[_],pathwaykolist_bmor[_],pathwaykolist_dme[_],\
                pathwaykolist_Msex[_],pathwaykolist_tca[_]]) >0:
                    fo.write(pathwaykolist[_] +"\t" + kosdic[pathwaykolist[_]]+"\t"+\
                    koecdic[pathwaykolist[_]]+"\t"+str(pathwaykolist_aga[_]) +"\t" + \
                    str(pathwaykolist_bmor[_]) +"\t" +str(pathwaykolist_dme[_]) +"\t" +\
                    str(pathwaykolist_Msex[_]) +"\t" +\
                    str(pathwaykolist_tca[_]) +"\t" + "\n")
            fo.close()

generateMapFiles()



folder = "E:\\store_for_D\\KOBAS2\\DB_using_20160104\\"
kos = open(folder+"Kos.txt").readlines()
koecdic ={}
import urllib.request
import re
import time
def getECnumber(konumber):
    with urllib.request.urlopen("http://rest.kegg.jp/find/ko/"+konumber) as response:
        ecs = re.findall("EC:[\d-]*\.[\d-]*\.[\d-]*\.[\d-]*", response.read().decode())
    response.close()
    ec =''
    for ele in ecs:
        ec = ec +ele+","
    return ec[:-1]
fopen = open(folder+"ko_ec.txt","w")
for _ in kos[1:]:
    time.sleep(3)
    tempko = _.split("\t")[0]
    tempec = getECnumber(_.split("\t")[1][:-1])
    print(tempko)
    fopen.write(tempko+"\t"+tempec+"\n")
    koecdic[_.split("\t")[0]] = getECnumber(_.split("\t")[1][:-1])
fopen.close()


"""
the code above does not work well. The website doesnot allow so many requests.
"""


import re
mylist = open("list.txt").readlines()
folder = "E:\\store_for_D\\KOBAS2\\DB_using_20160104\\"
fopen = open(folder+"ko_ec.txt","w")
for ele in mylist:
    tempko = re.findall("K\d{5} ",ele)
    if len(tempko) == 1:
        myko = tempko[0][:-1]
        tempec = re.findall("\[EC:[\d-]*\.[\d-]*\.[\d-]*\.[\d-]*.*\]",ele)
        if len(tempec) == 1:
            myec = tempec[0][1:-1]
            fopen.write(myko+"\t"+myec+"\n")
fopen.close()






"""
I have this table
ko	name	ec	aga	bmor	dme	Msex	tca
K00001	E1.1.1.1, adh	EC:1.1.1.1	0	0	1	0	0
K00002	AKR1A1, adh	EC:1.1.1.2	0	0	1	0	0
K00016	LDH, ldh	EC:1.1.1.27	1	1	1	1	1
K00121	frmA, ADH5, adhC	EC:1.1.1.284 1.1.1.1	1	1	1	1	1

change it to: based on EC numbers.


"""


import os
folder1 = "E:\\store_for_D\\KOBAS2\\DB_using_20160104\\PathwayNohuman\\"
folder2 = "E:\\store_for_D\\KOBAS2\\DB_using_20160104\\PathwayNohumanEC\\"
myfiles = os.listdir(folder1)
def transformFromKO2EC(filename):
    """
    given a file in folder1, output a file of the same name in folder2
    change the first column to EC, instead of KO
    """
    mylist = open(folder1+filename).readlines()
    eclist = []
    for ele in mylist[1:]:
        eclist.extend(ele.split("\t")[2][3:].split(" "))
    ecdic = {}
    for ele in eclist:
        ecdic[ele] =["",0,0,0,0,0]
    for key in ecdic:
        for ele in mylist[1:]:
            eleecs = ele.split("\t")[2][3:].split(" ")
            if key in eleecs:
                elelist = ele.split("\t")
                ecdic[key][0] += elelist[0]+" "
                for i in range(5):
                    ecdic[key][i+1] += int(elelist[i+3])
    fo = open(folder2+filename,"w")
    fo.write("EC\tKO\taga\tbmor\tdme\tMsex\ttca\n")
    for key in ecdic:
        fo.write(key+"\t"+ecdic[key][0][:-1]+"\t"+str(ecdic[key][1])+"\t"+\
        str(ecdic[key][2])+"\t"+str(ecdic[key][3])+"\t"+str(ecdic[key][4])+"\t"+\
        str(ecdic[key][5])+"\n")
    fo.close()

for myfile in myfiles:
    transformFromKO2EC(myfile)
