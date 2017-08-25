# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 17:57:43 2017

@author: k
"""

def SRP1_HP5_match20170716():
    '''
    match peaks to possible fragments of SRP1
    '''
    #new peaks after mix HP5 and proSRP1
    peaks = [1636.25830186, 1840.32080662, 2095.42465747, 2594.30975289, 2814.71737553, 2906.50214659, 3053.5291964, 4109.11402151, 4672.04049909, 4682.80815472, 5154.23860096]
    
    #get all possible fragments from proSRP1. cut after histidine, lysine or argine
    SRP1seq = 'MGSSHHHHHHDYDIPTTENLYFQGHMKTKEFPLPAEGIVMEFLVKNKSDMAHPISTLKMDIDQLTKKIIVTQSEDKMPVVRHLTPKPDVPSFENRFGVRVGTCPSGYVRRGTFCFPDDDYLE'
    aa_cut = "HKR"
    seqlen = len(SRP1seq)
    cut_positon = [-1]
    for i in range(seqlen):
        if SRP1seq[i] in aa_cut:
            cut_positon.append(i)
    cut_positon.append(seqlen-1)
    from itertools import combinations
    cut_choices = len(cut_positon)
    fragments =[SRP1seq[i+1:j+1] for (i,j) in combinations(cut_positon,2)]
    fragments_pos = [(i+1,j+1) for (i,j) in combinations(cut_positon,2)]
    #fragments that can be connected by cysteines
    cys_comball = list(combinations(range(len(fragments)),2))
    cys_comb = []
    for e in cys_comball:
        i = fragments_pos[e[0]]
        j = fragments_pos[e[1]]
        if i[1] <j[0]:
            if "C" in fragments[e[0]] and "C" in fragments[e[1]]:
                cys_comb.append(fragments[e[0]]+'-'+fragments[e[1]])
    
    fragments_all = fragments + cys_comb
#    fragments_random = [SRP1seq[i:j] for (i,j) in combinations(range(seqlen+1),2)]
#    fragments_all = fragments_all + fragments_random

    #calculate all mz of all fragments, charge 1
    from pyteomics import mass
    dc_mass = {}
    for seq in fragments_all:
        dc_mass[seq] = mass.calculate_mass(seq.replace('-',''),charge = 1)
    
    dc_peak = {}
    ppm = 300
    for peak in peaks:
        dc_peak[peak] = []
        for seq in fragments_all:
            if peak * (1-ppm/100000) < dc_mass[seq] and dc_mass[seq] < peak * (1+ppm/100000):
                dc_peak[peak].append((seq,dc_mass[seq]))
    