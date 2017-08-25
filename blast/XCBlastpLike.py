# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 13:20:28 2016
a blast-like program for proteins
@author: k
"""

AA2VALUE = {} #dictionary. 'A': 1, 'B':2, ..., 'Z': 26
for num in range(65,65+26):
    AA2VALUE[chr(num)] = num - 64

def aas2value(seq):
    """
    given a amino acid sequence, return a value by converting each character to int based on AA2VALUE
    "A" -> 1
    "AA" -> 1*27 +1 = 28
    """
    if type(seq) != str:
        print("input is not a string!")
        return None
    value = 0
    for ele in seq:
        if ele not in AA2VALUE:
            print("invalid letter exists!")
            return None
        else:
            value = value*27 + AA2VALUE[ele]
    return value

def value2aas(value):
    """
    given a value, return amino acid sequence, based on definition above
    """
    if type(value) != int:
        print("input value not int!")
        return None
    if value <= 0:
        print("too small value of ", value)
        return None
    aas = ""
    while value != 0:
        value, letternum = divmod(value,27)
        aas = chr(letternum +64) + aas
    return aas

def seq2kmer(seq,kmerlen=5):
    """
    given a seq, return a set of kmers with length of kmerlen
    """
    kmerset = set()
    for i in range(len(seq)+1-kmerlen):
        kmernum = seq[i:i+kmerlen]
        kmerset.add(kmernum)
    return kmerset

def seqs2kmerdic(seqs, kmerlen=5):
    """
    given many seqs, return a dictionary with keys of all kmers from seqs as key, and tuple of index of seqs with the kmer
    >>>seqs = ['AAA','AAB']
    >>>seqs2kmerdic(seqs,2)
    >>>{'AA':(0,1),'AB':(1)}
    """
    kmerdic ={}
    for num in range(len(seqs)):
        seq = seqs[num]
        kmers = seq2kmer(seq,kmerlen)
        for kmer in kmers:
            if kmer not in kmerdic:
                kmerdic[kmer] = set()
            kmerdic[kmer].add(num)
    for ele in kmerdic:
        kmerdic[ele] = tuple(kmerdic[ele])
    return kmerdic

def seqFindtargets(seq, kmerdic):
    """
    given a seq and a kmerdic generated from seqs2kmerdic, return a Counter object,
    with target id in seqs and number of common kmers.
    Counter object is sorted according to the number of common kmers.
    """
    for kmer in kmerdic: #get kmer length in kmerdic
        kmerlen = len(kmer)
        break
    kmers_fromseq = seq2kmer(seq, kmerlen)
    seqtargets = ()
    for kmer in kmers_fromseq:
        if kmer in kmerdic:
            seqtargets = seqtargets + kmerdic[kmer]
    from collections import Counter
    seqtargets = Counter(seqtargets)
    seqtargets = seqtargets.most_common()
    return seqtargets

def proteinPairwiseAlignGlobal(seq1,seq2,gap_open = -22, gap_extend = 0, matrix = 62):
    """
    given two protein sequence, return two sequence of global alignment.
    seq1 = "MPKSSSNDLP"
    seq2 = "MPKASSNADLP"
    return:
    MPKSSSN-DLP
    MPKASSNADLP
    """
    from Bio.SubsMat import MatrixInfo
    if matrix == 62:
        matrix = MatrixInfo.blosum62
    else:
        print("you need to change source code of function proteinPairwiseAlignGlobal")
        return None
    from Bio import pairwise2
    alns = pairwise2.align.globalds(seq1,seq2,matrix,gap_open,gap_extend)
    top_aln = alns[0]
    return top_aln[0], top_aln[1]


def proteinAlignLength(seq1_aln, seq2_aln,mincommon = 10,error_rate = 0.02):
    """
    return alignment length based on the output of proteinPairwiseAlignGlobal
    seq1_aln = "MPKSSSN-DLP"
    seq2_aln = "MPRASSNADLP"
    proteinAlignLength(seq1_aln,seq2_aln,1,0), return 8
    proteinAlignLength(seq1_aln,seq2_aln,3,0), return 6
    proteinAlignLength(seq1_aln,seq2_aln,3,0.5), return 10
    """
    seqlen = len(seq1_aln)
    if seqlen != len(seq2_aln):
        print("input two seqs are not the same length!")
        return None
    seqcommon = "" #find common elements. if two aa not the same, use #. "MP##SSN-AP"
    for num in range(seqlen):
        if seq1_aln[num] == seq2_aln[num]:
            seqcommon = seqcommon + seq1_aln[num]
        else:
            if seq1_aln[num] == "-" or seq2_aln[num] == "-":
                seqcommon = seqcommon + "-"
            else:
                seqcommon = seqcommon + "#"
    import re
    seqs = re.split('-+|#{2,}',seqcommon) #split at gap region or if two amino acids are not the same
    common = 0
    for seq in seqs:
        seqlen = len(seq)
        maxerror = int(seqlen * error_rate)
        error = seq.count("#") 
        if error > maxerror:
            seqlen = seqlen - error +maxerror
        if seqlen < mincommon:
            seqlen = 0
        common += seqlen
#        print(seq,seqlen)
    return common
        
