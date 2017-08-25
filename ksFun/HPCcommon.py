# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 13:55:00 2016

@author: k
"""

def compressPigz(files = 'list.txt',qsubnum = 1, qsubname = 'compress.pbs', compress = True, ppn = 12, core = 12):
    """
    files store files to be compressed if compress is True, else, decompress
    output a script for qsub jobs. default, output one job file.
    default qsubname is compress.txt
    ppn, numbers of CPUs to use
    core, numbers of thread for pigz
    """
    
    filenames = open(files).read().split("\n")
    if filenames[-1] == '':
        filenames = filenames[:-1]
    for i in range(qsubnum):
        fo = open(qsubname + str(i),'w')
        fo.write("#!/bin/bash\n")
        fo.write("#PBS -q batch\n")
        fo.write("#PBS -N compress\n")
        fo.write("#PBS -l nodes=1:ppn=%d\n"%ppn)
        fo.write("#PBS -l walltime=120:00:00\n")
        fo.write("#PBS -m abe -M atps@outlook.com\n")
        fo.write("#PBS -j oe\n")
        for j in range(i,len(filenames),qsubnum):
            if compress:
                c = ""
            else:
                c = " -d "
            fo.write("/panfs/panfs.cluster/home/ks2073/p/pigz-2.3.3/pigz -p %d "%core +c + filenames[j] +'\n')
        fo.close()

    
    
    