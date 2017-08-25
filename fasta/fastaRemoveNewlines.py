# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 12:24:23 2015

@author: k
"""

def baseEncode(number, base=36):
    """
    given a int number, return a string based on base
    10:A, 36:Z, 37:Z1, under default setting
    """
    if base == 10:
        return str(number)
    if not isinstance(number, int):
        raise TypeError('number must be an integer')
    alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    if base > 62 or base <=1:
        print("base should be between 2 and 62")
        return None
    sign = ""
    if number < 0:
        sign = "-"
        number = -number
    alphabet = alphabet[:base+1]
    if 0 <= number and number <base:
        return sign+alphabet[number]
    numberbase=""
    while number != 0:
        number, i = divmod(number, base)
        numberbase = alphabet[i] + numberbase
    return sign+numberbase

def fastaRenameRemoveNewlines(filename, usebaseEncode = False, minlength = 0):
    """
    given a filename of fasta file, output filename_out, remove newlines
    minlength of fasta sequence is 0bp in default.
    use number as the name in default. or use number change to base of 36
    """
    fo = open(filename,"r")
    fout = open(filename+"_out", "w")
    aline = fo.readline()
    faSeq = fo.readline()[:-1]
    readnum = 0
    if usebaseEncode == False:
        base = 10
    else:
        base = 36
    while aline:
        if aline[0] == ">":
            if faSeq != "":
                readnum += 1
                if len(faSeq) >= minlength:
                    fout.write(">"+baseEncode(readnum, base)+"\n"+faSeq+"\n")
            faSeq = ""
        else:
            faSeq = faSeq+aline[:-1]
        aline = fo.readline()
    fout.write(">"+baseEncode(readnum, base)+"\n"+faSeq+"\n")
    fout.close()
    fo.close()

fastaRenameRemoveNewlines("/scratch/ks2073/Transcriptome/AnophelesGambiae/unmapped/single.fa", True, 50)
        