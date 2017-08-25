# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 11:42:01 2016

@author: k
"""

import sys
sys.path.append("D:\\P\\3Language\\Xiaolong\\python\\XCProject\\fasta")
import largeFastaFile as lf


folder = "E:\\store_for_D\\DrFeng\\Assemblies\\TrinityTranslation\\"
filename = "Trinity.fasta.transdecoder.pep"

myfasta = lf.open_fasta_to_list(folder+filename)
mynewfasta = lf.fasta_within_seq_big_faster(myfasta)
lf.SaveFastaListToFile(mynewfasta, folder + "TrinityPepNoneWithin.fasta")

fo = open("TrinityPepNoneWithin2.fasta","w")
for myseq in mynewfasta:
    fo.write('>'+myseq.id+'\n'+str(myseq.seq)+'\n')
fo.close()