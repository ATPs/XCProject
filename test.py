# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:56:41 2016

@author: k
"""

fn = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/data/20160720sequencingDepth.txt0'
fno = '/panfs/panfs.cluster/scratch/ks2073/Transcriptome/ManducaSexta/sequencingDepth/dataBPKM/20160720sequencingDepth.txt0'
ls = open(fn).readlines()
ls = map(int, ls)
def div(num):
    return int(round(num /2491948975 * 1e9,0))
ls = map(div, ls)
fout = open(fno,'w')
for ele in ls:
    fout.write('%d\n'%ele)
fout.close()