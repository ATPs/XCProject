# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:01:32 2016

@author: k
"""

#tempfile = "E:\\store_for_D\\Transcriptome\\KEGG_Annotation2016\\KEGG_batchInfo\\KO.txt"
#dicKO2name ={}
#templist = open(tempfile).readlines()
#for tempele in templist:
#    tempelel = tempele.split("\t")
#    dicKO2name[tempelel[0][3:]] = tempelel[1][:-1]

#tempfile = "E:\store_for_D\Transcriptome\KEGG_Annotation2016\KEGG_batchInfo\geneKO2Pathway\\KOpathway.txt"
#dicKO2Pathway = {}
#dicPathway2KO = {}
#templist = open(tempfile).readlines()
#for tempele in templist:
#    tempelel = tempele.split("\t")
#    tempKO = tempelel[0][3:]
#    tempPathway = tempelel[1][5:-1]
#    if tempKO not in dicKO2Pathway:
#        dicKO2Pathway[tempKO] = set()
#    dicKO2Pathway[tempKO].add(tempPathway)
#    if tempPathway not in dicPathway2KO:
#        dicPathway2KO[tempPathway] = set()
#    dicPathway2KO[tempPathway].add(tempKO)

tempfile = "E:\store_for_D\Transcriptome\KEGG_Annotation2016\ManducaKEGG_KO\\20160129ManducaMCOTwithKO.txt"


import sqlite3
conn = sqlite3.connect(r'E:\store_for_D\Transcriptome\KEGG_Annotation2016\KEGG_batchInfo\KEGG_7insects.db') #create Database
c = conn.cursor()



c.execute("CREATE TABLE KO2Annotation (KO TEXT, Annotation TEXT);") #create table KO2Annotation
conn.commit()
tempfile = "E:\\store_for_D\\Transcriptome\\KEGG_Annotation2016\\KEGG_batchInfo\\KO.txt"
templist = open(tempfile).readlines()
for i in range(len(templist)): # add information in the KO.txt file to the table
    ele = templist[i]
    command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
    command = """INSERT INTO KO2Annotation (KO, Annotation) VALUES ('%s','%s');""" %(command2,command3)
    c.execute(command)
conn.commit()



c.execute("CREATE TABLE KO2Pathway (KO TEXT, Pathway TEXT);") #create table KO2Pathway
tempfile = "E:\store_for_D\Transcriptome\KEGG_Annotation2016\KEGG_batchInfo\geneKO2Pathway\\KOpathway.txt"
templist = open(tempfile).readlines()
for i in range(len(templist)): # add information in the KOpathway.txt file to the table
    ele = templist[i]
    command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
    command = """INSERT INTO KO2Pathway (KO, Pathway) VALUES ('%s','%s');""" %(command2,command3)
    c.execute(command)
conn.commit()



c.execute("CREATE TABLE Pathway2Annotation (Pathway TEXT, Annotation TEXT);") #create table Pathway2Annotation
tempfile = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\KEGG_batchInfo\pathway.txt"
templist = open(tempfile).readlines()
for i in range(len(templist)): # add information in the pathway.txt file to the table
    ele = templist[i]
    command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
    command = """INSERT INTO Pathway2Annotation (Pathway, Annotation) VALUES ('%s','%s');""" %(command2,command3)
    c.execute(command)
conn.commit()



c.execute("CREATE TABLE Pathway2Gene (Pathway TEXT, Gene TEXT);") #create table KO2Pathway
tempfolder = "E:\\store_for_D\\Transcriptome\\KEGG_Annotation2016\\KEGG_batchInfo\\geneKO2Pathway\\"
import os
tempfiles = os.listdir(tempfolder)
for tempfile in tempfiles:
    if "KOpathway.txt" not in tempfile:
        templist = open(tempfolder+tempfile).readlines()
        for i in range(len(templist)): # add information in the KOpathway.txt file to the table
            ele = templist[i]
            command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
            command = """INSERT INTO Pathway2Gene (Pathway, Gene) VALUES ('%s','%s');""" %(command2,command3)
            c.execute(command)
conn.commit()



c.execute("CREATE TABLE KO2Gene (KO TEXT, Gene TEXT);") #create table KO2Pathway
tempfolder = "E:\\store_for_D\\Transcriptome\\KEGG_Annotation2016\\KEGG_batchInfo\\KO2gene\\"
import os
tempfiles = os.listdir(tempfolder)
for tempfile in tempfiles:
    templist = open(tempfolder+tempfile).readlines()
    for i in range(len(templist)): # add information in the KOpathway.txt file to the table
        ele = templist[i]
        command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
        command = """INSERT INTO KO2Gene (KO, Gene) VALUES ('%s','%s');""" %(command2,command3)
        c.execute(command)
conn.commit()


c.execute("CREATE TABLE ManducaGene2KO (Gene TEXT, KO TEXT);") #create table Pathway2Annotation
tempfile = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\ManducaKEGG_KO\20160129ManducaMCOTwithKO.txt"
templist = open(tempfile).readlines()
for i in range(len(templist)): # add information in the pathway.txt file to the table
    ele = templist[i]
    command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
    command = """INSERT INTO ManducaGene2KO (Gene, KO) VALUES ('%s','ko:%s');""" %(command2,command3)
    c.execute(command)
conn.commit()



c.execute("CREATE TABLE KO2EC (KO TEXT, EC TEXT);") #create table Pathway2Annotation
tempfile = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\KEGG_batchInfo\ko2ec.txt"
templist = open(tempfile).readlines()
for i in range(len(templist)): # add information in the pathway.txt file to the table
    ele = templist[i]
    command2, command3 = ele.split('\t')[0].replace("'","''"),ele.split('\t')[1].replace("\n","").replace("'","''")
    command = """INSERT INTO KO2EC (KO, EC) VALUES ('%s','ko:%s');""" %(command2,command3)
    c.execute(command)
conn.commit()


conn.close()



"""
get the KOs of Manduca sexta in different maps
"""
conn = sqlite3.connect(r'E:\store_for_D\Transcriptome\KEGG_Annotation2016\KEGG_batchInfo\KEGG_7insects.db') #create Database
c = conn.cursor()

selected = c.execute("""SELECT DISTINCT ManducaGene2KO.KO 
                    FROM ManducaGene2KO JOIN KO2Pathway 
                    ON ManducaGene2KO.KO = KO2Pathway.KO 
                    WHERE KO2Pathway.Pathway = 'path:map01100'; """)
selectedlist = selected.fetchall()

selected2 = c.execute("""SELECT DISTINCT KO2Gene.KO 
                    FROM KO2Gene JOIN KO2Pathway 
                    ON KO2Gene.KO = KO2Pathway.KO 
                    WHERE KO2Pathway.Pathway = 'path:map01100' AND KO2Gene.Gene LIKE 'bmor:%'; """)
fo = open("list2.txt","w")
for ele in selected2.fetchall():
    fo.write(ele[0]+"\n")
fo.close()

print(len(selected2.fetchall()))
selectedlist = selected.fetchall()

selected = c.execute("""SELECT DISTINCT KO2Gene.KO 
                    FROM KO2Gene JOIN KO2Pathway 
                    ON KO2Gene.KO = KO2Pathway.KO 
                    WHERE KO2Pathway.Pathway = 'path:map01100'; """)




"""
deal with KEGG PATHWAY html file
"""



#fh: file, html
fh_ref = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolismKEGG\KEGG PATHWAY_ Metabolic pathways - Reference pathway.html"
fh_api = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolism6Species\KEGG PATHWAY_ Metabolic pathways - Acyrthosiphon pisum (pea aphid).html"
fh_ame = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolism6Species\KEGG PATHWAY_ Metabolic pathways - Apis mellifera (honey bee).html"
fh_aga = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolism6Species\KEGG PATHWAY_ Metabolic pathways - Anopheles gambiae (mosquito).html"
fh_bmor = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolism6Species\KEGG PATHWAY_ Metabolic pathways - Bombyx mori (domestic silkworm).html"
fh_dme = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolism6Species\KEGG PATHWAY_ Metabolic pathways - Drosophila melanogaster (fruit fly).html"
fh_tca = r"E:\store_for_D\Transcriptome\KEGG_Annotation2016\metabolism6Species\KEGG PATHWAY_ Metabolic pathways - Tribolium castaneum (red flour beetle).html"
#dic_file is the file names stored in the dictionary
dic_file={}
dic_file["ref"] = fh_ref
dic_file["api"] = fh_api
dic_file["ame"] = fh_ame
dic_file["aga"] = fh_aga
dic_file["bmor"] = fh_bmor
dic_file["dme"] =fh_dme
dic_file["tca"] =fh_tca























"""
I have got KEGG pictures with 7 insect genes labeled. I want to change the color schemes to make it more significant
"""

DIC_COLOR_OLD = {"ame":(191,255,191),"api":(255,187,204),"aga":(187,204,255),"dme":(207,255,207),"tca":(255,207,239),"bmor":(207,239,255),"msex":(223,239,207)}
DIC_COLOR_NEW = {"ame":(200,100,0),"api":(0,255,255),"aga":(255,0,255),"dme":(0,255,0),"tca":(232,27,0),"bmor":(255,255,0),"msex":(100,100,255)}

def KeggFigureColorChange(filename,outfilename=None):
    from PIL import Image
    im = Image.open(filename)
    width, height = im.size
    pixels = im.load()
    for species in DIC_COLOR_OLD:
        for x in range(width):
            for y in range(height):
                if pixels[x,y] == DIC_COLOR_OLD[species]:
                    pixels[x,y] = DIC_COLOR_NEW[species]
    if outfilename == None:
        outfilename = filename
    im.save(outfilename)

folder = "E:\\store_for_D\\Transcriptome\\KEGG_Annotation2016\\7Species_Insect\\"
import glob
files = glob.glob(folder+"*/*.png") + glob.glob(folder+"*/*/*.png")
for file in files:
    KeggFigureColorChange(file)
