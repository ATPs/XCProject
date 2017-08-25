'''
functions to analyze alignment related problems
'''
import pandas as pd
from Bio import SeqIO

def readMSA_fromfile(filename,fmt = 'fasta'):
    '''
    return a dataframe, given a multiple sequence alignment filename
    protein name as the index, columns is the site of residues, start from 0
    '''
    seqlist = list(SeqIO.parse(open(filename),fmt))
    index = [e.id for e in seqlist]
    columns = range(len(seqlist[0].seq))
    data = [list(str(e.seq)) for e in seqlist]
    df = pd.DataFrame(data = data, index = index, columns = columns)
    return df

def msa_common_position(msa,ratio = 1,print_len = True):
    '''
    return locations where most sequence have value instead of gap
    '''
    sites = []
    for i in range(msa.shape[1]):
        ls = list(msa.iloc[:,i])
        if (1 - ls.count('-')/msa.shape[0]) >= ratio:
            sites.append(i)
    if print_len:
        print(len(sites))
    return sites

def save_msa_fasta(msa,filename = 'seqs.fa'):
    '''
    save msa to fasta format
    '''
    fout = open(filename,'w')
    for head in msa.index:
        seq = list(msa.loc[head,:])
        fout.write('>'+head+'\n'+''.join(seq)+'\n')
    fout.close()