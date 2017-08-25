# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 13:57:07 2016

@author: k
"""

def seqRegExtract(filename,re_str,outname = "out.fasta",fmt = None,length_min = 0, length_max = 1e10):
    '''
    given a fastafilename, output a file in defined by fmt.
    if fmt is None, write fasta format, with no newline symbols for AA sequence
    the sequence length is in the defined length range
    only keep sequences
    '''
    from Bio import SeqIO
    import re
    f = open(filename,'r')
    fout = open(outname,'w')
    seqs = SeqIO.parse(f,'fasta')
    count = 0
    for seq in seqs:
        seqaa = str(seq.seq)
        if len(seqaa) >= length_min and len(seqaa) < length_max:
            if re.findall(re_str,seqaa):
                count += 1
                if fmt is None:
                    fout.write('>'+seq.id+'\n'+seqaa+'\n')
                else:
                    SeqIO.write(seq,fout,fmt)
    print("%d sequences fit the parameters"%(count))

def downloadFastaFromNCBI(lsID, database = 'nucleotide', outFilename = 'sequences.txt', fmt = 'fasta', nameOnly = True):
    """
    given a list of GI numbers or accession numbers, download the sequence from NCBI
    default database is nucleotide
    the default format is fasta, and default output name of sequences.txt
    """
    from Bio import Entrez, SeqIO
    fout = open(outFilename,'w')
    for num in range(0,len(lsID),100):
        sub_ls = lsID[num:min(num+100, len(lsID))]
        sub_ls = [str(ele) for ele in sub_ls]
        sub_str = ','.join(sub_ls)
        handle = Entrez.efetch(db = database, id = sub_str, rettype = fmt, retmode = 'text')
        if not nameOnly:
            fout.write(handle.read())
            handle.close()
        else:
            record = SeqIO.parse(handle, fmt)
            for ele in record:
                fout.write(ele.description+'\n')
            
        print(min(num+100, len(lsID)), 'finished')
    fout.close()
        
            
    
        