# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 17:39:27 2016

@author: k
"""
def errorMatch(seq1, seq2, errors=2):
    """
    Given seq1 and seq2, len(seq1) <= len(seq2), return whether they match each other with allowed error.
    """
    if len(seq1) > len(seq2):
        return False
    step = len(seq1)//(errors+1)
    if step == 0:
        return True
    if errors == 0:
        return seq1 in seq2
    parts = [seq1[i:i+step] for i in range(0,len(seq1),step)] #separate seq1 to error+1 parts
    if len(parts[-1]) < step:
        parts[-2] = parts[-2]+parts[-1]
        parts.pop()
    similar = False
    sameslist =[]
    for i in range(errors+1):
        findsame = seq2.find(parts[i])
        if  findsame >= 0:
            similar = True
            sameslist.append((i,findsame))
    if similar == False:
        return False
    for i,j in sameslist:
        if j-step*(i)>=0 and j-step*(i)+len(seq1) <= len(seq2):
            seq2n = seq2[j-step*i:j-step*i +len(seq1)]
            missmatched = 0
            for k in range(len(seq1)):
                if seq1[k] != seq2n[k]:
                    missmatched += 1
                if missmatched > errors:
                    break
            if missmatched <= errors:
                return True
    return False


fasta = r"E:\store_for_D\Insects Info\Anopheles gambiae\AllTrinity.fasta99.transdecoder.pepselected"
myfasta = list(SeqIO.parse(fasta,"fasta"))
import time
time1 = time.time()
dickmer3 = {}
for dummyi in range(len(myfasta)):
    seq = myfasta[dummyi]
    for i in range(len(seq.seq)+1-3):
        kmer3 = str(seq.seq[i:i+3])
        if kmer3 not in dickmer3:
            dickmer3[kmer3] = set()
        dickmer3[kmer3].add(dummyi)
print(time.time()-time1)

time1 = time.time()
for kmer3 in dickmer3:
    dickmer3[kmer3] = list(dickmer3[kmer3])
print(time.time()-time1)

#count=0
#for ele in dickmer3:
#    if 0 in dickmer3[ele]:
#        count+=1
#print(count)

time1 = time.time()
toremove = set()
from collections import Counter
error_rate = 0.02
for num1 in range(2):
    seq1 = myfasta[num1]
    seq1kmers = set()
    for i in range(len(seq1.seq)+1-3):
        seq1kmers.add(str(seq1.seq[i:i+3]))
#    print(time.time()-time1)
    seq1targets = []
    for kmer3 in seq1kmers:
        seq1targets += dickmer3[kmer3]
    seq1targets = Counter(seq1targets)
    seq1targets = seq1targets.most_common()
#    print(time.time()-time1)
    errors = int(len(seq1.seq)*error_rate)
    for seq2id, seq2_counts in seq1targets:
        seq2 = myfasta[seq2id]
        if seq2_counts >= len(seq1kmers)-errors*3 and seq2id != num1 and seq2id not in toremove and len(seq1.seq) <= len(seq2.seq) \
        and seq2_counts >=10:
#            print(time.time()-time1)
            if errorMatch(str(seq1.seq),str(seq2.seq),errors):
                toremove.add(num1)
                break
#    print(num1,time)

print(time.time()-time1)
print(toremove)

time1 = time.time()
toremove = set()
from collections import Counter
error_rate = 0.02
for num1 in range(20):
    seq1 = str(myfasta[num1].seq)
    errors = int(len(seq)*error_rate)
    for num2 in range(len(myfasta)):
        if num1!=num2:
            seq2 = str(myfasta[num2].seq)
            if seq2 not in toremove:
                if errorMatch(seq1,seq2,errors):
                    toremove.add(num1)
print(time.time()-time1)
print(toremove)




from Bio import SeqIO
fasta = r"E:\store_for_D\Insects Info\Anopheles gambiae\AllTrinity.fasta99.transdecoder.pepselected"
myfasta = list(SeqIO.parse(fasta,"fasta"))
import time
time1 = time.time()
kmerlen=5
dickmernum = {}
for dummyi in range(len(myfasta)):
    seq = myfasta[dummyi]
    for i in range(len(seq.seq)+1-kmerlen):
        kmernum = str(seq.seq[i:i+kmerlen])
        if kmernum not in dickmernum:
            dickmernum[kmernum] = set()
        dickmernum[kmernum].add(dummyi)
print(time.time()-time1)


time1 = time.time()
for kmernum in dickmernum:
    dickmernum[kmernum] = list(dickmernum[kmernum])
print(time.time()-time1)


time1 = time.time()
toremove = set()
from collections import Counter
error_rate = 0.02
for num1 in range(len(myfasta)):
    seq1 = myfasta[num1]
    seq1kmers = set()
    for i in range(len(seq1.seq)+1-kmerlen):
        seq1kmers.add(str(seq1.seq[i:i+kmerlen]))
#    print(time.time()-time1)
    seq1targets = []
    for kmernum in seq1kmers:
        seq1targets += dickmernum[kmernum]
    seq1targets = Counter(seq1targets)
    seq1targets = seq1targets.most_common()
#    print(time.time()-time1)
    errors = int(len(seq1.seq)*error_rate)
    for seq2id, seq2_counts in seq1targets:
        seq2 = myfasta[seq2id]
        if seq2_counts >= len(seq1kmers)-errors*kmerlen and seq2id != num1 and seq2id not in toremove and len(seq1.seq) <= len(seq2.seq) \
        and seq2_counts >=10:
#            print(time.time()-time1)
            if errorMatch(str(seq1.seq),str(seq2.seq),errors):
                toremove.add(num1)
                break
    if num1 %10000 ==0:
        print(num1,time.time()-time1)

print(time.time()-time1)
print(len(toremove))

fout = open(fasta + "less0.02error.txt","w")
for i in range(len(myfasta)):
    if i not in toremove:
        SeqIO.write(myfasta[i],fout,"fasta")
fout.close()



from Bio import SeqIO
fasta = r"E:\store_for_D\Insects Info\Anopheles gambiae\AllTrinity.fasta99.transdecoder.pepselected"
myfasta = list(SeqIO.parse(fasta,"fasta"))
import time
time1 = time.time()
aa = set()
for i in range(len(myfasta)):
    seq = str(myfasta[i].seq)
    for a in seq:
        aa.add(a)
print(len(aa))
print(time.time()-time1)
