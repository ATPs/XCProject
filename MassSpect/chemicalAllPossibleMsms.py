# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 17:53:16 2016

@author: k
play with mgf file, calculate all possible m/z's for a given molecule
"""

"""
examples
from pyteomics import mass
molecule = 'C2H5OH'
mass.calculate_mass(formula = molecule,charge = 1)

mass.calculate_mass(composition={'C':32,'H':55,'N':9,'O':15}, charge = 1)
mass.calculate_mass(sequence='PEPTIDE')
"""
import numpy as np
from pyteomics import mass
AA_comp = dict(mass.std_aa_comp)
AA_comp['d'] = mass.Composition({'H':-1}) #dyhydrogen
AA_comp['c'] = mass.Composition({'H':3, 'C':2, 'N':1, 'O':1}) #carbamidomethyl
AA_comp['o'] = mass.Composition({'O':1}) #oxidation
AA_comp['a'] = mass.Composition({'H':2, 'C':2, 'O':1}) #acetylation
DCaaModi = {'C':'c', 'M':'o', 'Nterminal':'a'}
DCaaModi2 = {'C':'d', 'M':'o', 'Nterminal':'a'}
neutron_w = 0.998
ISOTOP = {'H':['H[1]', 'H[2]'], 'C':['C[12]','C[13]'], 'N':['N[14]','N[15]'],
           'O':['O[16]','O[17]','O[18]'], 'S':['S[32]','S[33]','S[34]']}

def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    for ms/ms data
    The function generates all possible m/z for fragments of types 
    `types` and of charges from 1 to `maxharge`.
    """
    for i in range(1, len(peptide)-1):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)

def getPepmassWithModi(peptide, charge, ion_type = 'M', neutron_num = None):
    """
    function the same as mass.calculate_mass, add several common modifications.
    extra neutrons can be added
    d: dehydrogen, -H
    c: Carbamidomethyl, +H(3) C(2) N O
    o: oxidation, +O
    a: acetylation, +H(2) C(2) O
    if neutron num is not None, return mz and its relative abundance to monoisotopic one
    """
    m = mass.calculate_mass(peptide, aa_comp = AA_comp, charge = charge, ion_type = ion_type)
    if neutron_num is None:
        return m
    comp = mass.Composition(peptide,aa_comp = AA_comp)
    com_ori = {}
    for key in comp:
        com_ori[ISOTOP[key][0]] = comp[key]
    if neutron_num == 0:
        return m, 1
    keys = list(comp.keys())
    p_ori = mass.isotopic_composition_abundance(composition = com_ori)
    
#    if neutron_num == 1:
#        com_dc = {}
#        for key in keys:
#            com_dc[key] = dict(com_ori)
#            com_dc[key][ISOTOP[key][0]] -= 1
#            com_dc[key][ISOTOP[key][1]] = 1
#        mzs = [mass.calculate_mass(com_dc[key], aa_comp = AA_comp, charge = charge, ion_type = ion_type) for key in keys]
#        ps = [mass.isotopic_composition_abundance(composition = com_dc[key]) for key in keys]
#        return sum(np.array(mzs) * np.array(ps)) / sum(ps), sum(ps) /p_ori
    
    #neutron_num >= 1
    keys1more = keys
    keys2more = [ele for ele in keys1more if ele in ['O','S']]
    onetwo = [(neutron_num -n*2, n) for n in range(5+1) if neutron_num - n *2 >= 0]
    import itertools
    from collections import Counter
    ls_atoms = [[key] * com_ori[key] for key in com_ori]
    ls_atoms = list(itertools.chain.from_iterable(ls_atoms))
    changes = []
    for n in onetwo:
        change1 = list(itertools.product(keys1more, repeat = n[0]))
        change2 = list(itertools.product(keys2more, repeat = n[1]))
        for m in range(len(change1)):
            change1[m] = [(ISOTOP[ele][0], ISOTOP[ele][1]) for ele in change1[m]]
        for m in range(len(change2)):
            change2[m] = [(ISOTOP[ele][0], ISOTOP[ele][2]) for ele in change2[m]]
        changes = changes + list(itertools.product(change1, change2))
    changes = [ele[0] + ele[1] for ele in changes]
    changesTemp = set()
    for ele in changes:
        ele.sort()
        changesTemp.add(tuple(ele))
    changes = changesTemp
    mzs = []
    ps = []
    comps = []
    for ele in changes:
        molenew = ';'.join(ls_atoms)
        for pair in ele:
            if pair[0] not in molenew:
                break
            molenew = molenew.replace(pair[0],pair[1],1)
        if pair[0] not in molenew:
            continue
        com_new = dict(Counter(molenew.split(';')))
        if com_new not in comps:
            comps.append(com_new)
            mz = mass.calculate_mass(com_new, aa_comp = AA_comp, charge = charge, ion_type = ion_type)
            p = mass.isotopic_composition_abundance(composition = com_new)
            mzs.append(mz)
            ps.append(p)
    return sum(np.array(mzs) * np.array(ps)) / sum(ps), sum(ps) /p_ori
        
    
        
            
            




def fragmentsModi(peptide, types=('b', 'y'), maxcharge=1, full_include = False, modification = {'C':'cC'}):
    """
    for ms/ms data. for peptides with modification
    The function generates all possible m/z for fragments of types 
    `types` and of charges from 1 to `maxharge`.
    """
    if full_include:
        start = 0
    else:
        start = 1
    for i in range(start, len(peptide)-1):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    pep = peptide[:i]
                    for key in modification:
                        pep = pep.replace(key, modification[key])
                    yield mass.calculate_mass(pep, ion_type=ion_type, charge=charge, aa_comp = AA_comp)
                else:
                    pep = peptide[i:]
                    for key in modification:
                        pep = pep.replace(key, modification[key])
                    yield mass.calculate_mass(pep, ion_type=ion_type, charge=charge, aa_comp = AA_comp)


def modiPeptide(peptide, dcaaModi = DCaaModi, random = True):
    """
    given some peptide sequences, modifiy amino acid .
    return all possible sequences.
    if random is True, modification may happen or not for each aa.
    modifications can be 'acdo' as defined in getPepmassWithModi.
    for dcaaModi, key can be signle AA, or multiple letters like "NM" to represent more AA.
    Nterminal modification is NOT included currently.
    for d, dehydrogen of C, must come in pairs
    Only one modification for one aa.
    >>>modiPeptide('FCAC',{'C','c'}, True)
    ['FCAC', 'FcCAC', 'FcCAcC', 'FCAcC']
    >>>modiPeptide('FCAC',{'C','c'}, False)
    ['FcCAcC']    
    """
    import itertools
    dcaaWithModi = {}
    lsaaWithModi = []
    lspep = list(peptide)
    for ele in dcaaModi:
        if ele != 'Nterminal':
            for aa in ele:
                dcaaWithModi[aa]= dcaaModi[ele] + aa
                lsaaWithModi.append(aa)
    for num in range(len(lspep)):
        if lspep[num] in lsaaWithModi:
            if random:
                lspep[num] = [lspep[num], dcaaWithModi[lspep[num]]]
            else:
                lspep[num] = [dcaaWithModi[lspep[num]]]
    to_return = [''.join(ele) for ele in itertools.product(*lspep)]
    return [ele for ele in to_return if ele.count('d') % 2 ==0]
    
