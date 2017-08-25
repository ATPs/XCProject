
from Bio import SeqIO

def reads_lengths(filename, fmt = "fasta"):
    """
    given fasta filename
    return a dictionary, with reads length as key, counts as value
    """
    myfa = SeqIO.parse(open(filename,"r"),fmt)
    mylens = []
    for aseq in myfa:
        mylens.append(len(aseq.seq))
    from collections import Counter
    mydic = dict(Counter(mylens))
    print("totally", len(mydic), "different lengths, total reads", len(mylens))
    return mydic

def reads_lengths_special_unmappedreads(filename,fmt = "fasta"):
    """
    special for reads that with counts in their name
    >1-1000
    AATTGG
    length is 6, counts is 1000, instead of 1
    """
    myfa = SeqIO.parse(open(filename,"r"),fmt)
    mydic ={}
    for aseq in myfa:
        aseqlen = len(aseq.seq)
        if aseqlen not in mydic.keys():
            mydic[aseqlen] = 0
        aseqcount = int(aseq.id.split("-")[1])
        mydic[aseqlen] += aseqcount
    return mydic




fopen = open("/scratch/ks2073/Transcriptome/20151002Manduca/Oases/out5/_25/Sequences","r")
fout = open("/scratch/ks2073/Transcriptome/20151002Manduca/Oases/out5/_25/Sequences.out","w")
aline = fopen.readline()
faSeq = ""
faHead = ""
readnum = 0
while aline:
    if aline[0] == ">":
        if faSeq != "":
            facounts = int(faHead.split("-")[1].split("_")[0])
            if faHead.split()[0].split("_")[1] == "0":
                for temp_i in range(facounts):
                    readnum += 1
                    fout.write(">"+str(readnum)+"\n"+faSeq)
        faHead = aline[1:]
        faSeq = ""
    else:
        faSeq += aline
    aline = fopen.readline()
fout.close()
    