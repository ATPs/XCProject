# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 09:59:16 2016

@author: k
"""

def prepareForSubmission20160610():
    folder = 'E:\\Lab\\works\\20121023_manduca_genome_new\\TranscriptAssembleNew\\201510NewAssembly\\'
    fn_trinity = 'Trinity_asPairNewName.fasta'
    from Bio import SeqIO
    ls_seqs = list(SeqIO.parse(open(folder+fn_trinity),'fasta'))
    
    fn_to_delete = 'NCBI_to_delete.txt'
    _to_deletelines = open(folder+fn_to_delete).readlines()
    ls_to_delete =[]
    for ele in _to_deletelines:
        ls_to_delete.append(ele.split('\t')[0])
    print(len(ls_to_delete))
    
    fn_to_trim = 'NCBI_to_trim.txt'
    _to_trimlines = open(folder+fn_to_trim).readlines()
    dc_to_trim = {}
    for ele in _to_trimlines:
        eles = ele.split('\t')
        elename = eles[0]
        elelength = int(eles[1])
        if len(eles[2].split('..')) >2:
            ls_to_delete.append(elename)
        else:
            [ele_start, ele_end] = eles[2].split('..')
            ele_start = int(ele_start)
            ele_end = int(ele_end)
            if ele_start == 1:
                dc_to_trim[elename] = [ele_end-1, elelength]
            elif ele_end == elelength:
                dc_to_trim[elename] = [0,ele_start]
            else:
                ls_to_delete.append(elename)
    ls_to_delete.append('TRINITY_DN84339_c0_g2_i3')
    ls_to_delete.append('TRINITY_DN82415_c1_g6_i1')
    ls_to_delete.append('TRINITY_DN85500_c6_g5_i1')
    print(len(ls_to_delete))
    
    ls_keep = []
    for ele in ls_seqs:
        if ele.id in dc_to_trim:
            ele.seq = ele.seq[dc_to_trim[ele.id][0]:dc_to_trim[ele.id][1]]
        if ele.id not in ls_to_delete:
            if len(ele.seq)>200:
                ls_keep.append(ele)
    
    fo_new = open('20160613manducaNCBI_corrected.fa','w')
    for ele in ls_keep:
        fo_new.write('>'+ele.id+'\n'+str(ele.seq)+'\n')
    fo_new.close()
    
