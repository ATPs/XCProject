# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 17:40:31 2015

@author: k
"""

def txt_rm_keyword(filename,keyword):
    """
    given a filename, output a newfile with name filename+new
    for line in file, if keyword in that line, do not output it
    """
    fo = open(filename,"r")
    fo2 = open(filename+'new',"w")
    line = fo.readline()
    while line:
        if keyword not in line and line != "\n":
            fo2.write(line)
        line=fo.readline()
    fo2.close()
    fo.close()


def simpleExtractLinesWithKeywords(filename,keywords):
    """
    given a filename of a txt file, extract lines with keywords to filename_keywords
    """
    fout = open(filename+"_keywords","w")
    mylist = open(filename).readlines()
    for ele in mylist:
        if keywords in ele:
            fout.write(ele)
    fout.close()


myfile =r"E:\store_for_D\XuesongMassRawData\txt\band\peptides.txt"
simpleExtractLinesWithKeywords(myfile,"AGAP001657-PA")