# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 15:42:21 2015

@author: k
"""

def removeSeqsInFile(filename,listname):
    """
    given a file of fasta format, and a list file is names to be removed
    output a filename_out file, to remove those un-wanted sequences
    """
    from Bio import SeqIO
    mylist = open(listname).read().split()
    myfas = list(SeqIO.parse(filename,"fasta"))
    fo = open(filename+"_out","w")
    for ele in myfas:
        if ele.id not in mylist:
            fo.write(">"+ele.description+"\n"+str(ele.seq)+"\n")
    fo.close()
    return None


filename = r"E:\store_for_D\Transcriptome\20150105MCOT31166proteinSequenceNewIDfasta.txt"
listname = "list.txt"
removeSeqsInFile(filename, listname)


folder = "E:\store_for_D\Transcriptome\\"
filename = "20150105MCOT30303proteinSequenceNewIDfastaRemoveOtherSpecies.txt"
from Bio import SeqIO
myfasta = list(SeqIO.parse(open(folder+filename),"fasta"))
for i in range(0,len(myfasta),5000):
    fo = open(folder+"part"+str(i)+".txt","w")
    for j in range(i,min(i+5000,len(myfasta))):
        SeqIO.write(myfasta[j],fo,"fasta")
    fo.close()

def signalP20190215():
    '''
    get signal peptide annotation for OGS2 and MCOT
    based on signalP program
    '''
    import requests
    
    payload = {
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.90 Safari/537.36',
        'Sentence': 'Once upon a time, there was a little red hat and a wolf.',
        'Constituents': 'on',
        'NullLinks': 'on',
        'AllLinkages': 'on',
        'LinkDisplay': 'on',
        'ShortLegth': '6',
        'PageFile': '/docs/submit-sentence-4.html',
        'InputFile': "/scripts/input-to-parser",
        'Maintainer': "sleator@cs.cmu.edu"
    }
    
    r = requests.post("http://www.link.cs.cmu.edu/cgi-bin/link/construct-page-4.cgi#submit", 
                      data=payload)
    
    print(r.text)
    
    pay = {
           'R2':'OEuk',
           'S1':'MKMASSLAFLLLNFHVSLFLVQLLTPCSAQFSVLGPSGPILAMVGEDADLPCHLFPTM SAETMELRWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSIL',
           
           }
    
    rr = requests.post("http://www.csbio.sjtu.edu.cn/bioinf/Signal-3Lv1/", 
                      data=pay)
    
    import requests
    
    species = 'OEuk'
    seq = '>1\nMKQHHRHVHLKLRPKLRSCCHCGVKVAGHLRAAHLEQAHGVPAPACGACGKKFNYPNQVLRHQKHFHMGEKSFKCQSCDMMFSSQANLTQHAIKHSTQREHKCEYCGKAFKWRKNLTTHIMMHLNNRRHVCSACDERFVQQTSLKYHITKRHPEMV*'
    # species and seq may be read from file. Here you may do similar checking of the seq to remove space, new lines, or the existence of invalid amino acid letter, and to make sure its length is above 50 aa or give warning and break.
    
    url = 'http://www.csbio.sjtu.edu.cn/cgi-bin/Signal3L.cgi'
    form_data = {'mode':'string', 'R2':species, 'S1':seq}
    
    response = requests.post(url, data = form_data)
    
    #signal peptide by signalP
    filename = r"E:\Insects\ManducaSexta\signalPeptide\ms_signalP"
    outfile = r"E:\Insects\ManducaSexta\signalPeptide\ms_signalP.tab"
    fout = open(outfile,'w')
    fout.write("ms_id\tStart_nonSignalAA\n")
    for line in open(filename):
        if len(line)<1:
            continue
        if line[0] == '#':
            continue
        es = line.split()
        if es[9] == 'Y':
            fout.write(es[0] + '\t' + es[2] + '\n')
    fout.close()
    
    import requests
    import time
    import re
    import random
    from Bio import SeqIO
    
    def getSignalPeptide(seq):
        '''
        seq is a amino acid string
        return None if there is no signal peptide, else ,return position of signal peptide, like 
        '''
        species = 'OEuk'
        # species and seq may be read from file. Here you may do similar checking of the seq to remove space, new lines, or the existence of invalid amino acid letter, and to make sure its length is above 50 aa or give warning and break.
        
        url = 'http://www.csbio.sjtu.edu.cn/cgi-bin/Signal3L.cgi'
        form_data = {'mode':'string', 'R2':species, 'S1':seq}
        response = requests.post(url, data = form_data)
        text = response.text
        if "According to Signal-3L engine for your selected species, your input sequence <font color='red' face='times new roman'>does not</font> include a signal peptide." in text:
            return None
        text_signal = re.findall("According to Signal-3L engine for your selected species, the signal peptide is:<font color=red face='times new roman'>\d*-\d*",text)[0]
        start_nonsing = re.sub("According to Signal-3L engine for your selected species, the signal peptide is:<font color=red face='times new roman'>\d*-","",text_signal)
        return int(start_nonsing) + 1
    
    file1 = r"E:\Insects\ManducaSexta\MCOT_OGS\20150105MCOT30303proteinSequenceNewIDfastaRemoveOtherSpecies.fasta"
    file2 = "E:\Insects\ManducaSexta\MCOT_OGS\ms_ogs_proteins.fa"
    file_signalP = r"E:\Insects\ManducaSexta\signalPeptide\ms_signalP.tab"
    outfile = r"E:\Insects\ManducaSexta\signalPeptide\ms_signal3L.tab"
    file_missing = r"E:\Insects\ManducaSexta\signalPeptide\ms_signal3L.missing"
    finished_ids = open(outfile).read().split()
    finished_ids = [e for e in finished_ids if e.startswith('Msex2') or e.startswith('MCOT')]
    fout = open(outfile,'w')
    fmissing = open(file_missing,'w')
    fout.write("ms_id\tStart_nonSignalAA\n")
    seqs_exclude = open(file_signalP).read().split()
    seqs_exclude = set([e for e in seqs_exclude if e.startswith('Msex2') or e.startswith('MCOT')])
    seqs_torun = list(SeqIO.parse(file1,'fasta')) + list(SeqIO.parse(file2,'fasta'))
    seqs_torun = [e for e in seqs_torun if e.id not in seqs_exclude]
    seqlen = len(seqs_torun)
    for n,s in enumerate(seqs_torun):
        if s.id in finished_ids:
            continue
        seq = str(s.seq).strip('*')
        try:
            start = getSignalPeptide(seq)
            print(n+1,'out of',seqlen,s.id,'finished, result:',start)
        except:
            fmissing.write('>'+s.id+'\n'+str(s.seq)+'\n')
            print(n+1,'out of',seqlen,s.id,'something wrong')
            continue
        #time.sleep(random.random())
        fout.write(s.id+'\t'+str(start)+'\n')
    #    if n==10:
    #        break
    fout.close()
    fmissing.close()

def summarizeSignalPresult20190303():
    '''
    summarize signalP signal3L result
    '''
    file_signalP = r"E:\Insects\ManducaSexta\signalPeptide\ms_signalP"
    
    
    
    