# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 10:52:39 2016

@author: k
"""


def patternCounter(text, pattern):
    """
    return number of pattern in text
    """
    count = 0
    lp = len(pattern)
    lt = len(text)
    for i in range(lt - lp+1):
        if pattern == text[i:i+lp]:
            count += 1
    return count

def frequentWords(text, k):
    """
    word length is k
    return the most frequent kmer
    """
    kmers = [text[i:i+k] for i in range(0,len(text) - k+1)]
    from collections import Counter
    kcount = Counter(kmers)
    kcount = kcount.most_common()
    if len(kcount) >0:
        m = kcount[0][1]
    mkers = []
    for ele in kcount:
        if ele[1] == m:
            mkers.append(ele[0])
        else:
            break
    return mkers

def complementDNA(text):
    """
    return a complement sequence of a DNA/RNA
    """
    ntDNA = {'A': 'T','T':'A','C':'G','G':'C'}
    reverse = ''
    for ele in text:
        reverse = ntDNA[ele] + reverse
    return reverse

def complementDNAFromFile(filename):
    """
    output a file with reverse complement
    """
    fo = open(filename)
    text = fo.read()
    fo.close()
    fout = open(filename + '.complement','w')
    reverse = ''
    ntDNA = {'A': 'T','T':'A','C':'G','G':'C'}
    for ele in text:
        if ele in ntDNA:
            reverse = ntDNA[ele] + reverse
        else:
            reverse = ele + reverse
    fout.write(reverse)
    fout.close()

def patternPosition(pattern, genome):
    """
    return a list of position of the pattern
    """
    import re
    pattern = re.compile(r'(?=('+pattern+'))')
    m = re.finditer(pattern, genome)
    starts = []
    for ele in m:
        starts.append(ele.start())
    return starts

def readFile(filename):
    "return lines in filename"
    fo = open(filename)
    ls = fo.readlines()
    ls = [ele.replace('\n','') for ele in ls]
    return ls

def clumpFinder(genome, k,L,t):
    """
    given a genome, and kmer length of k, clump length of L, and k repeat t times in L
    return the kmer seq
    """
    dckmer = {}
    for i in range(0,len(genome) +1 -k):
        kmer = genome[i:i+k]
        if kmer not in dckmer:
            dckmer[kmer] = [[],[]]
        dckmer[kmer][0].append(i)
        dckmer[kmer][1].append(i+k)
    goodkmer = []
    for kmer in dckmer:
        kmergood = False
        kmercount = len(dckmer[kmer][0])
        if kmercount >= t:
            for i in range(0, kmercount +1 - t):
                boundaries = [dckmer[kmer][0][i:i+t],dckmer[kmer][1][i:i+t]]
                if boundaries[1][-1] - boundaries[0][0] <= L:
                    kmergood = True
                    break
            if kmergood:
                goodkmer.append(kmer)
    return goodkmer

def pattern2number(pattern):
    '''
    A:0,C:1,G:2,T:3
    given pattern, return number
    '''
    dc = {'A':0,'C':1,'G':2,'T':3}
    num = 0
    for ele in pattern:
        n = dc[ele]
        num = num *4 + n
    return num

def number2pattern(number, k):
    '''
    given a number, return kmer. length of kmer is k
    '''
    dc = {0:'A',1:'C',2:'G',3:'T'}
    pattern = ''
    while number != 0:
        n = number % 4
        number = number //4
        pattern = dc[n] + pattern
    if len(pattern) <= k:
        pattern = 'A'*k+pattern
        pattern = pattern[-k:]
    return pattern

def computingFreq(seq, k):
    '''
    k is the kmer length.
    return a list of count for different kmers.
    kmers ordered from small to big, defined by pattern2number function
    '''
    count = [0 for _i in range(4**k)]
    for i in range(len(seq)+1 - k):
        kmer = seq[i:i+k]
        n = pattern2number(kmer)
        count[n] += 1
    return count

def pattern2number_cursive(pattern):
    '''
    A:0,C:1,G:2,T:3
    given pattern, return number
    call recursive
    '''
    dc = {'A':0,'C':1,'G':2,'T':3}
    if pattern == '':
        return 0
    else:
        symbol = pattern[-1]
        prefix = pattern[:-1]
        return 4 * pattern2number_cursive(prefix) + dc[symbol]



#20161216 week2
def skewCount(seq):
    '''
    given a seq, return a list of skew
    skew: occurrrences of G minus occurences of C
    CGA
    return [0,-1,0,0]
    '''
    ls = [0,]
    c = 0
    g = 0
    for ele in seq:
        if ele == 'C':
            c += 1
        if ele == 'G':
            g += 1
        ls.append(g - c)
    return ls

def findSkewMin(seq):
    '''
    given a seq, return the positon of minimun skew number
    '''
    l = skewCount(seq)
    m = min(l)
    p = []
    for n in range(len(l)):
        if l[n] == m:
            p.append(n)
#    print(' '.join(str(ele) for ele in p))
    return p

def hammingdistance(p,q):
    '''
    p and q are two sequence with same length.
    return the number of non-identical letters
    '''
    l = len(p)
    count = 0
    for n in range(l):
        if p[n] != q[n]:
            count += 1
    return count

def approxiMatch(kmer, seq, error):
    '''
    return position where kmer can match with seq with maximun error mismatches
    Sample Input:
    ATTCTGGA
    CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
    3
    Sample Output:
    6 7 26 27
    '''
    l = len(kmer)
    ps = []
    for n in range(len(seq) + 1 - l):
        s = seq[n:n+l]
        e = hammingdistance(kmer,s)
        if e <= error:
            ps.append(n)
#    print(' '.join(str(ele) for ele in ps))
    return ps

def approxiMatchCount(kmer, seq, error):
    '''
    return position where kmer can match with seq with maximun error mismatches
    Sample Input:
    ATTCTGGA
    CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
    3
    Sample Output:
    6 7 26 27
    '''
    l = len(kmer)
    ps = []
    for n in range(len(seq) + 1 - l):
        s = seq[n:n+l]
        e = hammingdistance(kmer,s)
        if e <= error:
            ps.append(n)
    return len(ps)


def mostFreqKmerFinder(seq,klen,error):
    '''
    find the most frequent kmer with klen as length, with maximum error number of mismatch
    '''
    kmers = [seq[i:i+klen] for i in range(len(seq))]
    from collections import Counter
    kcount = Counter(kmers)
    kcount = dict(kcount)
    kmers = list(kcount.keys())
    def kerror(kmer, error):
        '''
        given a kmer, return a list of kmers with less than error differences
        kerror('AT',1)
        return: ['AT','TT','GT','CT','AC','AG,'AA']
        '''
        import itertools
        changes = itertools.combinations(range(len(kmer)), error)
        ks = []
        for change in changes:
            k_ls = []
            for n in range(len(kmer)):
                if n not in change:
                    k_ls.append(kmer[n])
                else:
                    k_ls.append('ATCG')
            ks += list(itertools.product(*k_ls))
        ks = ks =[''.join(ele) for ele in ks]
        return set(ks)
    dc_kmerkerror = {}
    for kmer in kmers:
        dc_kmerkerror[kmer] = kerror(kmer,error)
    dc_kerrorCount = {}
    for kmer in dc_kmerkerror:
        for k in dc_kmerkerror[kmer]:
            if k not in dc_kerrorCount:
                dc_kerrorCount[k] = 0
            dc_kerrorCount[k] += kcount[kmer]
    m = max(dc_kerrorCount.values())
    ls = []
    for k in dc_kerrorCount:
        if dc_kerrorCount[k] == m:
            ls.append(k)
#    print(m)
#    print(' '.join(ls))
    return ls


def mostFreqKmerFinderBothStrand(seq,klen,error):
    '''
    find the most frequent kmer with klen as length, with maximum error number of mismatch
    '''
    kmers = [seq[i:i+klen] for i in range(len(seq))]
    seq2 = complementDNA(seq)
    kmers2 = [seq2[i:i+klen] for i in range(len(seq2))]
    kmers = kmers + kmers2
    from collections import Counter
    kcount = Counter(kmers)
    kcount = dict(kcount)
    kmers = list(kcount.keys())
    def kerror(kmer, error):
        '''
        given a kmer, return a list of kmers with less than error differences
        kerror('AT',1)
        return: ['AT','TT','GT','CT','AC','AG,'AA']
        '''
        import itertools
        changes = itertools.combinations(range(len(kmer)), error)
        ks = []
        for change in changes:
            k_ls = []
            for n in range(len(kmer)):
                if n not in change:
                    k_ls.append(kmer[n])
                else:
                    k_ls.append('ATCG')
            ks += list(itertools.product(*k_ls))
        ks = ks =[''.join(ele) for ele in ks]
        return set(ks)
    dc_kmerkerror = {}
    for kmer in kmers:
        dc_kmerkerror[kmer] = kerror(kmer,error)
    dc_kerrorCount = {}
    for kmer in dc_kmerkerror:
        for k in dc_kmerkerror[kmer]:
            if k not in dc_kerrorCount:
                dc_kerrorCount[k] = 0
            dc_kerrorCount[k] += kcount[kmer]
    m = max(dc_kerrorCount.values())
    ls = []
    for k in dc_kerrorCount:
        if dc_kerrorCount[k] == m:
            ls.append(k)
#    print(m)
#    print(' '.join(ls))
    return ls

def find_DnaA_Salmonella_enterica():
    filename = r'F://Salmonella_enterica.txt'
    from Bio import SeqIO
    seq = str(SeqIO.read(open(filename),'fasta').seq)
    skew = findSkewMin(seq)
    seq2 = seq[min(skew) - 500: max(skew) +500]
    dnaA = mostFreqKmerFinderBothStrand(seq2,9,1)
    return dnaA
#    print(dnaA)



#20161219 week3
def kerror(kmer, error):
    '''
    given a kmer, return a list of kmers with less than error differences
    kerror('AT',1)
    return: ['AT','TT','GT','CT','AC','AG,'AA']
    '''
    import itertools
    changes = itertools.combinations(range(len(kmer)), error)
    ks = []
    for change in changes:
        k_ls = []
        for n in range(len(kmer)):
            if n not in change:
                k_ls.append(kmer[n])
            else:
                k_ls.append('ATCG')
        ks += list(itertools.product(*k_ls))
    ks = ks =[''.join(ele) for ele in ks]
    return set(ks)

def MotifEnumeration(klen,error,seqs):
    '''
    return kmers if it appears in every string with at most 'error' errors
    '''
    kmers = [set([seq[i:i+klen] for i in range(len(seq)+1-klen)]) for seq in seqs]
    kerrors = [set() for i in range(len(seqs))]
    for i in range(len(kmers)):
        for k in kmers[i]:
            kes = kerror(k,error)
            kerrors[i] = kerrors[i].union(kes)
    kcommon = kerrors[0]
    for i in range(len(kmers)):
        kcommon = kcommon & kerrors[i]
#    print(' '.join(kcommon))
    return kcommon

def entropySinglesite(ps):
    '''
    A: 0.2; C: 0.1; G: 0; T: 0.7
    single site, entropy H = -sum(p*log2(p))
    '''
    import numpy as np
    s = 0
    def e(p):
        '''
        enttropy for single value
        '''
        if p == 0:
            return 0
        else:
            return -p * np.log2(p)
    for p in ps:
        s += e(p)
#    print(s)
    return s


def hammingdistanceKwithSeq(kmer,seq):
    '''
    p and q are two sequence with same length. p is kmer, q is the longer seq. return the minimus score between p and q
    return the number of non-identical letters
    '''
    klen = len(kmer)
    seqlen = len(seq)
    if seqlen < klen:
        return None
    distances = []
    for i in range(seqlen +1 - klen):
        s = seq[i:i+klen]
        d = hammingdistance(kmer, s)
        distances.append(d)
    return min(distances)

def medianString(klen, seqs):
    '''
    Input: An integer k, followed by a collection of strings Dna.
    Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern. (If there are
    multiple such strings Pattern, then you may return any one.)
    Sample Input:
    3
    AAATTGACGCAT
    GACGACCACGTT
    CGTCAGCGCCTG
    GCTGAGCACCGG
    AGTTCGGGACAG
    Sample Output:
    GAC
    '''
    distance = float('inf')
    median = []
    import itertools
    for kk in itertools.product('ATCG',repeat = klen):
        kmer = ''.join(kk)
        d = 0
        for seq in seqs:
            d += hammingdistanceKwithSeq(kmer,seq)
        if d < distance:
            median = []
            median.append(kmer)
            distance = d
        elif d == distance:
            median.append(kmer)
#        break
    return median


def kmerProbability(kmer,matrix):
    '''
    serve for the function profileMostProbableKmer. given a kmer and matrix, return the probability
    '''
    dc = {'A': 0,'C': 1,'G':2,'T':3}
    product = 1
    for i in range(len(kmer)):
        base = kmer[i]
        product = product * matrix[dc[base]][i]
    return product

def profileMostProbableKmer(seq,klen,matrix):
    '''
    Given a profile matrix Profile, we can evaluate the probability of every k-mer in a string Text and find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated by Profile among all k-mers in Text. For example, ACGGGGATTACC is the Profile-most probable 12-mer in GGTACGGGGATTACCT. Indeed, every other 12-mer in this string has probability 0. In general, if there are multiple Profile-most probable k-mers in Text, then we select the first such k-mer occurring in Text.
    
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    Sample Input:
    ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
    5
    0.2 0.2 0.3 0.2 0.3
    0.4 0.3 0.1 0.5 0.1
    0.3 0.3 0.5 0.2 0.4
    0.1 0.2 0.1 0.1 0.2
    Sample Output:
    CCGAG
    '''
    kmers = list([seq[i:i+klen] for i in range(len(seq)+1-klen)])
    ps = []
    for kmer in kmers:
        ps.append((kmer, kmerProbability(kmer,matrix)))
    pmax = max([i[1] for i in ps])
    for ele in ps:
        if ele[1] == pmax:
            return ele[0]



def seq2matrix(seqs):
    '''
    given seqs of the same length, build a matrix shown in profileMostProbableKmer. 
    the rows are 'A','C','G','T'
    '''
    nseq = len(seqs)
    mseq = len(seqs[0])
    matrix = [[0 for i in range(mseq)] for j in range(4)]
    for col in range(mseq):
        s = ''
        for row in range(nseq):
            s += seqs[row][col]
        matrix[0][col] = s.count('A') /nseq
        matrix[1][col] = s.count('C') /nseq
        matrix[2][col] = s.count('G') /nseq
        matrix[3][col] = s.count('T') /nseq
    return matrix

def seq2matrixWithPseudocounts(seqs):
    '''
    given seqs of the same length, build a matrix shown in profileMostProbableKmer. 
    the rows are 'A','C','G','T'
    '''
    nseq = len(seqs)
    mseq = len(seqs[0])
    matrix = [[0 for i in range(mseq)] for j in range(4)]
    for col in range(mseq):
        s = ''
        for row in range(nseq):
            s += seqs[row][col]
        matrix[0][col] = (s.count('A')+1) /(nseq + 4)
        matrix[1][col] = (s.count('C')+1) /(nseq + 4)
        matrix[2][col] = (s.count('G')+1) /(nseq + 4)
        matrix[3][col] = (s.count('T')+1) /(nseq + 4)
    return matrix
    

def scoringMatrix(matrix):
    '''
    entropySinglesite calculates entropy for each site.
    after have a matrix profile, calculate entropy for the whole matrix
    '''
    entropy = 0
    for ps in matrix:
        entropy += entropySinglesite(ps)
    return entropy

def greedyMotifSearch(seqs,klen,t):
    '''
    '''
    bestMotifs = [seq[:klen] for seq in seqs]
    bestScore = scoringMatrix(seq2matrix(bestMotifs))
    for kmer in [seqs[0][i:i+klen] for i in range(len(seqs[0]) + 1 - klen)]:
        motifs = []
        motifs.append(kmer)
        for i in range(1,len(seqs)):
            matrix = seq2matrix(motifs)
            m = profileMostProbableKmer(seqs[i],klen,matrix)
            motifs.append(m)
        score = scoringMatrix(seq2matrix(motifs))
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
    return bestMotifs

def greedyMotifSearchWithPseudocounts(seqs,klen,t):
    '''
    '''
    bestMotifs = [seq[:klen] for seq in seqs]
    bestScore = scoringMatrix(seq2matrixWithPseudocounts(bestMotifs))
    for kmer in [seqs[0][i:i+klen] for i in range(len(seqs[0]) + 1 - klen)]:
        motifs = []
        motifs.append(kmer)
        for i in range(1,len(seqs)):
            matrix = seq2matrixWithPseudocounts(motifs)
            m = profileMostProbableKmer(seqs[i],klen,matrix)
            motifs.append(m)
        score = scoringMatrix(seq2matrixWithPseudocounts(motifs))
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
    return bestMotifs


#20161220 week4

def seq2kmer(seq,klen):
    '''
    given a seq, return a list of kmers with length klen
    '''
    l = len(seq)
    return [seq[i:i+klen] for i in range(l+1-klen)]


def scoreMotifs(motifs):
    '''
    given a group of motifs, return a score.
    first calculate consensus string, then use sonsensus string to compute the sum of hamming distances
    '''
    n = len(motifs)
    m = len(motifs[0])
    score = 0
    for i in range(m):
        s = ''.join(motif[i] for motif in motifs)
        score = score + n - max(s.count(w) for w in 'ATCG')
    return score
        

def randomizedMotifSearch(seqs, klen, t):
    '''
    In general, we can begin from a collection of randomly chosen k-mers Motifs in Dna, 
    construct Profile(Motifs), and use this profile to generate a new collection of k-mers:
    Motifs(Profile(Motifs), Dna).
    Why would we do this? Because our hope is that Motifs(Profile(Motifs), 
    Dna) has a better score than the original collection of k-mers Motifs. 
    We can then form the profile matrix of these k-mers,
    Profile(Motifs(Profile(Motifs), Dna))
    and use it to form the most probable k-mers,
    Motifs(Profile(Motifs(Profile(Motifs), Dna)), Dna).
    We can continue to iterate. . .
    Profile(Motifs(Profile(Motifs(Profile(Motifs), Dna)), Dna))...
    for as long as the score of the constructed motifs keeps improving, 
    which is exactly what RandomizedMotifSearch does. 
    To implement this algorithm, you will need to randomly select the initial collection 
    of k-mers that form the motif matrix Motifs. To do so, you will need a random number 
    generator (denoted Random(N)) that is equally likely to return any integer from 1 to N. 
    You might like to think about this random number generator as an unbiased N-sided die.
    '''
    if t != len(seqs):
        print('something wrong with t')
        return 0
    kmers = [seq2kmer(seq,klen) for seq in seqs]
    import random
    bestMotifs = [random.choice(k) for k in kmers]
    bestprofile = seq2matrixWithPseudocounts(bestMotifs)
    bestScore = scoreMotifs(bestMotifs)
    for n in range(100000):
        motifs = [profileMostProbableKmer(seq,klen,bestprofile) for seq in seqs]
        profile = seq2matrixWithPseudocounts(motifs)
        score = scoreMotifs(motifs)
        if score < bestScore:
            bestScore = score
            bestprofile = profile
            bestMotifs = motifs
        else:
            break
#    print(bestMotifs,bestScore)
    return bestMotifs, bestScore
    
def randomizedMotifSearchMany(seqs,klen, t, n):
    '''
    repeat function randomizedMotifSearch n times, return the best motifs
    '''
    bestMotifs, bestScore = randomizedMotifSearch(seqs,klen,t)
    for i in range(1,n):
        motifs, score = randomizedMotifSearch(seqs,klen,t)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
    return bestMotifs, bestScore



def randomizedMotifSearchGibbs(seqs, klen, t, n):
    '''
    the same funtion as randomizedMotifSearchMany. Change only one motif in motifs with gibbs sampling
    '''
    kmers = [seq2kmer(seq,klen) for seq in seqs]
    seql = len(seqs[0])
    import numpy as np
    bestMotifs = [np.random.choice(k) for k in kmers]
    bestprofile = seq2matrixWithPseudocounts(bestMotifs)
    bestScore = scoreMotifs(bestMotifs)
    for j in range(n):
        i = np.random.randint(t)
        motifs = bestMotifs.copy()
        ps = [kmerProbability(kmer,bestprofile) for kmer in kmers[i]]
        ps = [p/sum(ps) for p in ps]
        motifs[i] = kmers[i][np.random.choice(range(seql + 1 - klen), p = ps)]
        profile = seq2matrixWithPseudocounts(motifs)
        score = scoreMotifs(motifs)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
            bestprofile = profile
#        else:
#            bestMotifs, bestScore
    return bestMotifs, bestScore
    
def randomizedMotifSearchManyGibbs(seqs, klen, t, n, N):
    '''
    run randomizedMotifSearchGibbs N times. return the best
    '''
    bestMotifs, bestScore = randomizedMotifSearchGibbs(seqs,klen,t,n)
    for i in range(1,N):
        motifs, score = randomizedMotifSearchGibbs(seqs,klen,t,n)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
    print('\n'.join(bestMotifs))
    return bestMotifs
