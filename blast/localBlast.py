# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 17:00:16 2016

@author: k

doing local blast
"""

BLASTFOLDER = "C:\\P\\blast-2.6.0+\\bin\\"
BLASTP = BLASTFOLDER+"blastp.exe"
BLASTN = BLASTFOLDER+"blastn.exe"
BLASTx = BLASTFOLDER+"blastx.exe"
TBLASTN = BLASTFOLDER+"tblastn.exe"
WINDOWMASKER = BLASTFOLDER+"windowmasker.exe"
MAKEBLASTDB = BLASTFOLDER+"makeblastdb.exe"
SEGMASKER = BLASTFOLDER + "segmasker.exe"


def blastDBbuild(fastafile, outname = None, infmt = "fasta", title = None, dbtype = "prot"):
    """
    make sure there is no space in the names
    given a fastafile, build a blast search database.
    fastafile is the location of fasta file. outname is the database name.
    title is the title for the database
    dbtype ="prot", default database is protein.
    """
    if title == None:
        title = fastafile
    if outname == None:
        outname = "blastDB"
    import subprocess
    command1 = SEGMASKER + " -in " + fastafile + " -infmt "+ infmt + " -parse_seqids  -outfmt maskinfo_asn1_bin  -out " + outname +".asnb"
    command2 = MAKEBLASTDB + " -in " + fastafile + " -input_type " + infmt + " -dbtype " + dbtype + " -parse_seqids -mask_data " + \
    outname +".asnb -out " + outname +" -title " + title
#    print(command1)
#    print(command2)
    subprocess.call(command1)
    subprocess.call(command2)


def localblastp(query,db,outname = None, outfmt = 0,evalue = 10, threads =3, matrix = "BLOSUM62",max_target_seqs = 20):
    """
    local blast. query is a filename of input. output blastp_result.txt by default. default format is 0
    """
    if outname == None:
        outname = "blastp_result.txt"
    import subprocess
    command = BLASTP + " -threshold 24 -task blastp-fast -db " + db + " -query " + query+ " -out " + outname + " -matrix " + matrix +\
    " -evalue " + str(evalue) + " -num_threads " + str(threads) + " -outfmt " + str(outfmt) +\
    " -max_target_seqs " + str(max_target_seqs)
    print(command)
    subprocess.call(command)

#fastafile = """E:\\Lab\\fastaDB\\Ag\\AgCufflinks20160301.fa.transdecoder.pepnonredun"""
#outname = "E:\\Lab\\fastaDB\\Ag\\CuffPrNonredunDB"
#blastDBbuild(fastafile,outname,"fasta","2015AgCuffNonredunDB")
#
#localblastp("D:\\P\\3Language\\Xiaolong\\python\\list.txt","E:\\Lab\\fastaDB\\Ag\\CuffPrNonredunDB","D:\\P\\3Language\\Xiaolong\\python\\blastp_result.txt")

