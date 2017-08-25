# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:09:08 2016

@author: k
"""

filename = "E:\\Lab\\works\\Temp\\cornucopia_3.gff" #file name of the original gff file.
mylist = open(filename,"r").readlines()

geneidname = {} 
# create a dictionary, key is like "maker-scaffold_16-exonerate_est2genome-gene-0.0", value like: "Dn_16_1.0-gene" 
for ele in mylist:
    if ele.count("\t") == 8:
        elements = ele.split("\t")
        if ele[:8] == "scaffold" and elements[2] == "gene":
            geneid = elements[-1].split(";")
            geneidname[geneid[0].split("=")[1]] = geneid[1][:-1].split("=")[1]


geneAnnotation = {}
annolist = open("E:\\Lab\\works\\Temp\\gene2name.txt").readlines()#Annotation file, open to list
for ele in annolist:
    geneAnnotation[ele.split("\t")[0]] =  ele.split("\t")[1][:-1]
#mytext = open(filename,"r").read()
#for ele in geneidname:
#    mytext.replace(ele, geneidname[ele])
#fout = open(filename+"out","w")
#fout.write(mytext)
#fout.close()

newlist =[]
for ele in mylist:
    if ele.count("\t") == 8:
        elements = ele.split("\t")
        if ele[:8] == "scaffold" and elements[1] == "maker":
            if elements[2] == "gene":
                geneid = elements[-1].split(";")
                key = geneid[0].split("=")[1] # get the key in geneidname
                elements[-1] = elements[-1].replace(key,geneidname[key])#relace old maker id to id in Name=Dn...
                keyAnno = geneidname[key].split("-")[0]
                elements[-1] = elements[-1].replace("Name="+keyAnno,"Name="+geneAnnotation[keyAnno])
                elen = ele.rsplit("\t",1)[0] + "\t" + elements[-1]
                newlist.append(elen)
            else:
                geneid = elements[-1].split(";",2)
                key = geneid[0][:geneid[0].rfind("-mRNA")].split("=")[1]
                elements[-1] = elements[-1].replace(key,geneidname[key])
                keyAnno = geneidname[key].split("-")[0]
                elements[-1] = elements[-1].replace("Name="+keyAnno,"Name="+geneAnnotation[keyAnno])
                elen = ele.rsplit("\t",1)[0] + "\t" + elements[-1]
                newlist.append(elen)
        else:
            elen = ele
            newlist.append(elen)
    else:
        elen = ele
        newlist.append(elen)

fout = open(filename+"new.gff2","w")#output gff filename.
for ele in newlist:
      fout.write(ele)
fout.close() #save the newlist to file.

#geneIDnamelist =[]
#for ele in mylist:
#    if ele.count("\t") == 8:
#        elements = ele.split("\t")
#        if ele[:8] == "scaffold" and elements[2] == "gene":
#            geneid = elements[-1].split(";")
#            geneIDnamelist.append((geneid[0].split("=")[1], geneid[1][:-1].split("=")[1]))
#
    
    
