# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 16:02:46 2017

@author: k
"""

def runInterpro3species20170604():
    '''
    split files in HPC, and generate qsub scripts
    '''
    import os
    from Bio import SeqIO
    num = 20
    folder_full = '/panfs/panfs.cluster/scratch/ks2073/Interpro/Full/'
    files = os.listdir(folder_full)
    folder_split = '/panfs/panfs.cluster/scratch/ks2073/Interpro/Split/'
    def splitfiles(filename, outname, num = 20):
        ls = list(SeqIO.parse(open(filename), 'fasta'))
        fouts = [open(outname+str(i),'w') for i in range(num)]
        for i in range(len(ls)):
            n = i % num
            to_write = '>'+ls[i].id+'\n'+str(ls[i].seq)+'\n'
            fouts[n].write(to_write)
        for n in range(num):
            fouts[n].close()
    for filename in files:
        splitfiles(folder_full + filename, folder_split+filename, num)
    
    #generateQsubScripts
    files_input = os.listdir(folder_split)
    folder_scripts = '/panfs/panfs.cluster/home/ks2073/cuffbuild/interproscan/2017MsTcAp/'
    for filename in files_input:
        fout = open(folder_scripts+filename,'w')
        fout.write('''
#!/bin/bash
#PBS -q batch
#PBS -N {0}
#name you want to give your job 
#PBS -l nodes=1:ppn=12
#request 1 nodes w/8 processors per node
#PBS -l walltime=120:00:00
# request 48 hour of  walltime.  your job will be killed if it goes over.  
#PBS -m abe -M atps@outlook.com
#PBS -j oe

module load jdk
module load python
module load perl
module load gcc-4.9.2
module load signalp/4.1
module load  blast+/2.6.0


cd /panfs/panfs.cluster/scratch/ks2073/Interpro/Result/

/panfs/panfs.cluster/opt/interproscan/5.17-56.0-64-bit/prebuilt/interproscan.sh -d ./  -goterms -iprlookup -pa -f TSV,XML,HTML -i {1}{0}
                   '''.format(filename,folder_split))
        fout.close()
    