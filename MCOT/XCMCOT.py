# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 16:38:04 2015

@author: k
"""
import sys
sys.path.append("D:\\P\\3Language\\Xiaolong\\python\\XCProject\\MCOT")

from Bio import SeqIO
from XCBlast import *
from largeFastaFile import *



class XCGene(object):
    """
    store the protein sequence, CDS sequence and transcript sequence in the same class
    for easier usage in the future
    the protein, cds and transcript sequences are stored in the seq format of Biopython
    in default, protein_seq, cds_seq and transcript_seq were set to None.
    """
    def __init__(self, name = None, protein_seq = None, cds_seq = None, transcript_seq = None):
        self.name = name
        self.protein_seq = protein_seq
        self.cds_seq = cds_seq
        self.transcript_seq = transcript_seq
    
    def protein_length(self):
        if self.protein_seq != None:
            return len(self.protein_seq.seq)
        else:
            return None
    
    def cds_length(self):
        if self.cds_seq != None:
            return len(self.cds_seq.seq)
        else:
            return None
    
    def transcript_length(self):
        if self.transcript_seq != None:
            return len(self.transcript_seq.seq)
        else:
            return None
        
    def __str__(self):
        """
        a string with the description of the gene, including details of protein, cds, transcripts and their sequence
        """
        line0 = self.name        
        if self.protein_seq == None:
            line1 = "protein sequence is missing!"
        else:
            line1 = "length of protein is "+ str(len(self.protein_seq.seq))
        if self.cds_seq == None:
            line2 = "cds sequence is missing!"
        else:
            line2 = "length of cds is" + str(len(self.cds_seq.seq))
        if self.transcript_seq == None:
            line3 = "transcript sequence is missing!"
        else:
            line3 = "length of transcript is " + str(len(self.transcript_seq.seq))
        return line0 + "\n" + line1 + "\n" + line2 + "\n" + line3

    def description(self):
        """
        print out descriptions of XCgene.
        >>> mygene.description()
        length of protein is 500
        length of cds is 1500
        length of transcript is 2000
        """
        print(self)
        
    def details(self):
        """
        return protein, cds and transcript sequences
        """
        return (self.protein_seq, self.cds_seq, self.transcript_seq)



    
    

def transdecoder_to_dic_XCGene(file_names):
    """
    given the file name of output of transdecoder, return a dictionary in XCgene format
    the name should be like: "20150709Diaphorina_citri_Transcripts.fasta"
    and the corresponding file names are:
    20150709Diaphorina_citri_Transcripts.fasta.transdecoder.cds
    20150709Diaphorina_citri_Transcripts.fasta.transdecoder.mRNA
    20150709Diaphorina_citri_Transcripts.fasta.transdecoder.pep
    """
    fo_protein = open(file_names + ".transdecoder.pep", "r")
    fo_cds = open(file_names + ".transdecoder.cds", "r")
    fo_mRNA = open(file_names + ".transdecoder.mRNA", "r")
    dic_protein = SeqIO.to_dict(SeqIO.parse(fo_protein, "fasta"))
    dic_cds = SeqIO.to_dict(SeqIO.parse(fo_cds, "fasta"))
    dic_mRNA =SeqIO.to_dict(SeqIO.parse(fo_mRNA, "fasta"))
    dic_XCGene = {}
    for name in dic_protein:
        dic_XCGene[name] = XCGene(name, dic_protein[name], dic_cds[name], dic_mRNA[name])
    return dic_XCGene

def maker_to_cuff(dic_maker2cuff_ml, dicGene_maker, dicGene_cuff, coverage = 0.7, min_length = 200):
    """
    given the dictionary of maker2cufflinks matched length, maker and cufflink sequences in the XCGene format    
    return a list of maker genes that will be kept as part of MCOT
    """
    keepmaker=[]
    makerlist=[]
    for elements in dic_maker2cuff_ml:
        makerlist.append(elements[0])
        if dic_maker2cuff_ml[elements] * (1/ (coverage * dicGene_maker[elements[0]].protein_length()) + 1/min_length) < 1:
            keepmaker.append((elements[0], "Maker, bad match in cufflinks"))
    makerset = set(makerlist)
    for elements in dicGene_maker:
        if elements not in makerset:
            keepmaker.append((elements,"Maker, no match in cufflinks"))
    return keepmaker


def cuff_to_denovo(file_maker, file_cuff, file_Trinity, file_Oases, \
file_maker2cuff, file_cuff2maker, file_maker2Trinity, file_Trinity2maker, file_maker2Oases, file_Oases2maker, \
file_cuff2Trinity, file_Trinity2cuff, file_cuff2Oases, file_Oases2cuff,\
file_maker2Uni, file_cuff2Uni, file_Trinity2Uni, file_Oases2Uni, file_uniprot, fout_mcot, coverage = 0.7, min_length =200,\
 file_cuffUniqueAndNoneWithin = None, file_TrinityUN = None, file_OasesUN = None):
    """
    given transdecoder file of cufflinks, Trinity, and Oases
    their blast results to each other and to uniprot
    return a list of seleted model with description
    save them to file_mcot
    """
    mcot=[]
    mcot_to_chose = {}
    dicGene_maker = transdecoder_to_dic_XCGene(file_maker)    
    dicGene_cuff = transdecoder_to_dic_XCGene(file_cuff)
    dicGene_Trinity = transdecoder_to_dic_XCGene(file_Trinity)
    dicGene_Oases = transdecoder_to_dic_XCGene(file_Oases)
    dicGene={}
    dicGene.update(dicGene_maker)
    dicGene.update(dicGene_cuff)
    dicGene.update(dicGene_Trinity)
    dicGene.update(dicGene_Oases)
    print("read in all gene sequences of maker, cufflinks, Trinity and Oases")
    
    diclen_maker = fasta_length(open_fasta_to_list(file_maker+".transdecoder.pep"))
    diclen_cuff = fasta_length(open_fasta_to_list(file_cuff+".transdecoder.pep"))
    diclen_Trinity = fasta_length(open_fasta_to_list(file_Trinity+".transdecoder.pep"))
    diclen_Oases = fasta_length(open_fasta_to_list(file_Oases+".transdecoder.pep"))
    diclen_all = {}
    diclen_all.update(diclen_maker)
    diclen_all.update(diclen_cuff)
    diclen_all.update(diclen_Trinity)
    diclen_all.update(diclen_Oases)
    diclen_Uniprot = fasta_length(open_fasta_to_list(file_uniprot))
    print("finishing generating dictionary of protein length of maker, cufflinks, Trinity, Oases, and Uniprot")
    
    dic_ml_m2c_ori = file_Blast6_to_best_match_length_list_refine(file_maker2cuff, file_cuff2maker)
    dic_ml_c2m = dic_keyoftupple_trans(file_Blast6_to_best_match_length_list_refine(file_cuff2maker, file_maker2cuff))
    dic_ml_m2c = dic_keyoftupple_trans(file_Blast6_to_best_match_length_list_refine(file_maker2cuff, file_cuff2maker))
    dic_ml_m2t = dic_keyoftupple_trans(file_Blast6_to_best_match_length_list_refine(file_maker2Trinity, file_Trinity2maker))
    dic_ml_m2o = dic_keyoftupple_trans(file_Blast6_to_best_match_length_list_refine(file_maker2Oases, file_Oases2maker))
    dic_ml_c2t = dic_keyoftupple_trans(file_Blast6_to_best_match_length_list_refine(file_cuff2Trinity, file_Trinity2cuff))
    dic_ml_c2o = dic_keyoftupple_trans(file_Blast6_to_best_match_length_list_refine(file_cuff2Oases, file_Oases2cuff))
    print("finishing generating matching dictionary of maker to cufflinks, trinity and oases, cufflinks to Trinity/ Oases")
    
    mcot_maker_ori = maker_to_cuff(dic_ml_m2c_ori, dicGene_maker, dicGene_cuff)
    for elements in mcot_maker_ori:
        mcot_to_chose[elements[0]] = []
        if elements[0] not in dic_ml_m2t:
            mtml = 0
        else:
            mtml = dic_ml_m2t[elements[0]][1]
        if mtml / dicGene_maker[elements[0]].protein_length() >= coverage:
            mcot_to_chose[elements[0]].append(dic_ml_m2t[elements[0]][0])
        if elements[0] not in dic_ml_m2o:
            moml = 0
        else:
            moml = dic_ml_m2o[elements[0]][1]
        if moml / dicGene_maker[elements[0]].protein_length() >= coverage:
            mcot_to_chose[elements[0]].append(dic_ml_m2o[elements[0]][0])
    print("finishing connecting maker-unique with trintiy and oases, matched length / query length > coverage (default: 0.7)")
    
    for elements in dicGene_cuff:
        mcot_to_chose[elements] = []
        if elements not in dic_ml_c2t:
            ctml = 0
        else:
            ctml = dic_ml_c2t[elements][1]
        if ctml / dicGene_cuff[elements].protein_length() >= coverage:
            mcot_to_chose[elements].append(dic_ml_c2t[elements][0])
        if elements not in dic_ml_c2o:
            coml = 0
        else:
            coml = dic_ml_c2o[elements][1]
        if coml / dicGene_cuff[elements].protein_length() >= coverage:
            mcot_to_chose[elements].append(dic_ml_c2o[elements][0])
    print("finish connecting cufflinks with trinity and oases, matched length / query length > coverage (default: 0.7)")
    
    dic_matchscore_maker = file_Blast6_to_Uniprot_matching_score(file_maker2Uni, diclen_maker, diclen_Uniprot)
    dic_matchscore_cuff = file_Blast6_to_Uniprot_matching_score(file_cuff2Uni, diclen_cuff, diclen_Uniprot)
    dic_matchscore_trinity = file_Blast6_to_Uniprot_matching_score(file_Trinity2Uni, diclen_Trinity, diclen_Uniprot)
    dic_matchscore_Oases = file_Blast6_to_Uniprot_matching_score(file_Oases2Uni, diclen_Oases, diclen_Uniprot)
    dic_matchscore={}
    dic_matchscore.update(dic_matchscore_maker)
    dic_matchscore.update(dic_matchscore_cuff)
    dic_matchscore.update(dic_matchscore_Oases)
    dic_matchscore.update(dic_matchscore_trinity)
    print("finishing calcualte matching score with Uniprot")
    
    def dicscore(dic_matchscore, elements):
        """
        if element in the dic_matchscore, return the matching score. else, return 0
        """
        if elements in dic_matchscore:
            return dic_matchscore[elements][2]
        else:
            return 0
        
    for elements in mcot_to_chose:
        if mcot_to_chose[elements] == []:
            mcot.append((elements, elements))
        elif len(mcot_to_chose[elements]) == 1:
            lcuff = dicGene[elements].protein_length()
            uniscore_cuff = dicscore(dic_matchscore, elements)
            if "X" in str(dicGene[elements].protein_seq.seq):
                lcuff = lcuff * 0.7
                uniscore_cuff = uniscore_cuff * 0.7
            ldenovo = dicGene[mcot_to_chose[elements][0]].protein_length()
            uniscore_denovo = dicscore(dic_matchscore,mcot_to_chose[elements][0])
            if "X" in str(dicGene[mcot_to_chose[elements][0]].protein_seq.seq):
                ldenovo = ldenovo * 0.7
                uniscore_denovo = uniscore_denovo * 0.7
            if lcuff >= ldenovo:
                if max(uniscore_denovo, uniscore_cuff) > 0.5 and uniscore_cuff - uniscore_denovo < -0.3:
                    mcot.append((elements, mcot_to_chose[elements][0]))
                else:
                    mcot.append((elements, elements))
            else:
                if max(uniscore_denovo, uniscore_cuff) > 0.5 and uniscore_cuff - uniscore_denovo > 0.3:
                    mcot.append((elements, elements))
                else:
                    mcot.append((elements, mcot_to_chose[elements][0]))
        else:
            lcuff = dicGene[elements].protein_length()
            uniscore_cuff = dicscore(dic_matchscore,elements)
            if "X" in str(dicGene[elements].protein_seq.seq):
                lcuff = lcuff * 0.7
                uniscore_cuff = uniscore_cuff * 0.7
            lTrinity = dicGene[mcot_to_chose[elements][0]].protein_length()
            uniscore_Trinity = dicscore(dic_matchscore,mcot_to_chose[elements][0])
            if "X" in str(dicGene[mcot_to_chose[elements][0]].protein_seq.seq):
                lTrinity = lTrinity * 0.7
                uniscore_Trinity = uniscore_Trinity * 0.7
            lOases = dicGene[mcot_to_chose[elements][1]].protein_length()
            uniscore_Oases = dicscore(dic_matchscore,mcot_to_chose[elements][1])
            if "X" in str(dicGene[mcot_to_chose[elements][1]].protein_seq.seq):
                lOases = lOases * 0.7
                uniscore_Oases = uniscore_Oases * 0.7
            if lcuff >= lTrinity and lTrinity >= lOases:
                (ss_cuff, ss_Trinity, ss_Oases) = (3,2,1) # ss, short for selection score. to simplify the code
            if lcuff >= lOases and lOases > lTrinity:
                (ss_cuff, ss_Trinity, ss_Oases) = (3,1,2)
            if lTrinity > lcuff and lcuff >= lOases:
                (ss_cuff, ss_Trinity, ss_Oases) = (2,3,1)
            if lTrinity >= lOases and lOases > lcuff:
                (ss_cuff, ss_Trinity, ss_Oases) = (1,3,2)
            if lOases > lcuff and lcuff >= lTrinity:
                (ss_cuff, ss_Trinity, ss_Oases) = (2,1,3)
            if lOases > lTrinity and lTrinity > lcuff:
                (ss_cuff, ss_Trinity, ss_Oases) = (1,2,3)
            if max(uniscore_cuff, uniscore_Trinity, uniscore_Oases) > 0.5:
                ss_cuff += 10/3 * uniscore_cuff
                ss_Trinity += 10/3 * uniscore_Trinity
                ss_Oases += 10/3 *uniscore_Oases
            ss_max = max(ss_cuff, ss_Trinity, ss_Oases)
            if ss_cuff == ss_max:
                mcot.append((elements, elements))
            elif ss_Trinity == ss_max:
                mcot.append((elements, mcot_to_chose[elements][0]))
            else:
                mcot.append((elements, mcot_to_chose[elements][1]))
    print("mcot choose the best finished!")
    
    mcot_nonredundant =[]
    if file_cuffUniqueAndNoneWithin == None:
        listGene_cuffUandNoneWithin = fasta_within_seq_big_fast(open_fasta_to_list(file_cuff+".transdecoder.pep"))
        dicGene_cuffUandNoneWithin = {}
        for cuffelements in listGene_cuffUandNoneWithin:
            dicGene_cuffUandNoneWithin[cuffelements.id] = cuffelements
    else:
        dicGene_cuffUandNoneWithin = open_fasta_to_dic(file_cuffUniqueAndNoneWithin)
    for elements in mcot:
        if elements[0] in dicGene_cuffUandNoneWithin:
            mcot_nonredundant.append(elements)
    print("finished removing redundant elements. gene number before and after in mcot is")
    print(len(mcot), " ", len(mcot_nonredundant))
    
    if file_TrinityUN == None:
        list_GeneTrinityUN = fasta_within_seq_big_fast(open_fasta_to_list(file_Trinity+".transdecoder.pep"))
    else:
        list_GeneTrinityUN = open_fasta_to_list(file_TrinityUN)
    if file_OasesUN == None:
        list_GeneOasesUN = fasta_within_seq_big_fast(open_fasta_to_list(file_Oases+".transdecoder.pep"))
    else:
        list_GeneOasesUN = open_fasta_to_list(file_OasesUN)
    print("finished reading/ converting Trinity / Oases Unique and NonWithin")
    
    list_denovoUN = list_GeneTrinityUN + list_GeneOasesUN
    dic_ml_t2c = file_Blast6_to_best_match_length_list_refine(file_Trinity2cuff, file_cuff2Trinity)
    dic_ml_o2c = file_Blast6_to_best_match_length_list_refine(file_Oases2cuff, file_cuff2Oases)
    dic_ml_t2m = file_Blast6_to_best_match_length_list_refine(file_Trinity2maker, file_maker2Trinity)
    dic_ml_o2m = file_Blast6_to_best_match_length_list_refine(file_Oases2maker, file_maker2Oases)
    dic_ml_denovo = {}
    dic_ml_denovo.update(dic_ml_t2c)
    dic_ml_denovo.update(dic_ml_o2c)
    dic_ml_denovo.update(dic_ml_t2m)
    dic_ml_denovo.update(dic_ml_o2m)
    print("finish matching length comparing of denovo to cufflinks and maker")
    
    denovo_goodmatch = []
    for elements in dic_ml_denovo:
        if dic_ml_denovo[elements] * (1/ (coverage * dicGene[elements[0]].protein_length()) + 1/min_length) >= 1:
            denovo_goodmatch.append(elements[0])
    denovo_goodmatch = set(denovo_goodmatch)

    dic_matchscoreFilter = {}
    for elements in dic_matchscore:
        if dic_matchscore[elements][2] > 0.5:
            dic_matchscoreFilter[elements] = dic_matchscore[elements]
    denovo_keep = []
    for elements in list_denovoUN:
        if elements.id in dic_matchscoreFilter and elements.id not in denovo_goodmatch and "type:complete" in elements.description:
            denovo_keep.append(elements.id)
    print("finishing find denovo unique genes. There are ", len(denovo_keep), "genes.")
    
    for elements in denovo_keep:
        mcot.append((elements, elements))
    print("raw number of genes in mcot is ", len(mcot))
    
    list_mcot_protein = []    
    for elements in mcot:    
        list_mcot_protein.append(dicGene[elements[1]].protein_seq)
    list_mcot_protein_UandNwithin3 = fasta_within_seq_big_faster(list_mcot_protein, 3, 3)
    dic_mcot_protein_UandNwithin3={}
    for elements in list_mcot_protein_UandNwithin3:
        dic_mcot_protein_UandNwithin3[elements.id] = elements
    mcotUN3 =[]
    for elements in mcot:
        if elements[1] in dic_mcot_protein_UandNwithin3:
            mcotUN3.append(elements)
    print("After removing similar sequences, the number of genes in mcot is ", len(list_mcot_protein_UandNwithin3))
        
    saveFastaListToFile(list_mcot_protein_UandNwithin3, fout_mcot+"mcot_before_rename.protein")
    fout_temp1 = open(fout_mcot+"list_of_mcot_ori_selected_for_analysis.txt", "w")
    for elements in mcotUN3:
        fout_temp1.write(str(elements) + "\n")
    fout_temp1.close()
    
    
    dic_protein2Genename ={} # the content is like {transcript_id: gene_id; ...}
    list_cuff = open_fasta_to_list(file_cuff)
    for elements in list_cuff:
        dic_protein2Genename[elements.id] = elements.description.split()[1]
    
    dic_mcotUN3={} # {mcot:[ori_id, gene_id]}, mcot, selected id; ori_id, original id before selection
    for elements in mcotUN3:
        if elements[0].split('|')[0] in dic_protein2Genename:
            dic_mcotUN3[elements[1]] = [elements[0], dic_protein2Genename[elements[0].split('|')[0]]]
        else:
            dic_mcotUN3[elements[1]] = [elements[0], None]
    print("finishing assign origin gene names")
    
    dic_mcotUN3gene = {} # {gene_id:[[mcot, ori_id], [mcot2, ori_id2],...]; ...}
    for elements in dic_mcotUN3:
        dic_mcotUN3gene[dic_mcotUN3[elements][1]] = []
    for elements in dic_mcotUN3:
        dic_mcotUN3gene[dic_mcotUN3[elements][1]].append([elements, dic_mcotUN3[elements][0]])
    
    genenum = 1
    mcot_final = [] #[('MCOT23245.2', ['TCONS_00043156|m.30467', 'TCONS_00043156|m.30467']),...]
    for elements in dic_mcotUN3gene:
        if elements == None:
            for item in dic_mcotUN3gene[elements]:
                mcot_final.append(("MCOT"+"{0:05}".format(genenum)+".0",item))
                genenum += 1
        else:
            transcriptnum = 1
            if len(dic_mcotUN3gene[elements]) == 1:
                transcriptnum = 0
            for item in dic_mcotUN3gene[elements]:
                mcot_final.append(("MCOT"+"{0:05}".format(genenum)+"." + str(transcriptnum),item))
                transcriptnum += 1
            genenum +=1
    
    def return_des(elements, dicofml):
        """
        if elements in dicofml, return subject_id, subject_length, and matched_length
        else, return "NA, NA, NA"
        """
        if elements in dicofml:
            return (dicofml[elements][0], str(diclen_all[dicofml[elements][0]]), str(dicofml[elements][1]))
        else:
            return ("NA", "NA", "NA")
    dic_uniprot = open_fasta_to_dic(file_uniprot)    
    
    fout_mcot_pr = open(fout_mcot+"protein.fa", "w")
    fout_mcot_cds = open(fout_mcot+"cds.fa", "w")
    fout_mcot_transcript = open(fout_mcot+"transcript.fa", "w")
    
    for elements in mcot_final:
        if elements[1][0] in dic_matchscore:
            des_uni_id = dic_matchscore[elements[1][0]][0]
            des_uni_len = str(dic_matchscore[elements[1][0]][4])
            des_uni_ml = str(dic_matchscore[elements[1][0]][1])
            des_uni_des = str(dic_uniprot[des_uni_id].description)
            des_uni = "uni: " + des_uni_id + " ul: " + des_uni_len + " ml_u: " + des_uni_ml + " "+ des_uni_des
        else:
            des_uni = "no good match in uniprot"
        short_description = elements[0]
        if elements[1][1] in dicGene_cuff:
            des_cuffid = elements[1][1]
            des_cufflen = str(diclen_all[elements[1][1]])
            (des_makerid, des_makerlen, des_ml_cm) = return_des(elements[1][1], dic_ml_c2m)
            (des_Trinityid, des_Trinitylen, des_ml_ct) = return_des(elements[1][1], dic_ml_c2t)
            (des_Oasesid, des_Oaseslen, des_ml_co) = return_des(elements[1][1], dic_ml_c2o)
            long_description = "ori: cuff " + des_cuffid + " clen: " + des_cufflen + \
            " maker: " + des_makerid + " ml: " + des_makerlen + " ml_cm: " + des_ml_cm +\
            " Trinity: " + des_Trinityid + " tl: " + des_Trinitylen + " ml_ct: " + des_ml_ct +\
            " Oases: " + des_Oasesid + " ol: " + des_Oaseslen + " ml_co: " + des_ml_co
            if elements[1][0] in dicGene_Trinity:
                short_description += ".CT"
            elif elements[1][0] in dicGene_Oases:
                short_description += ".CO"
            elif elements[1][0] in dicGene_cuff:
                short_description += ".CC"
            else:
                short_description +=".immpossible"
        elif elements[1][1] in dicGene_maker:
            des_makerid = elements[1][1]
            des_makerlen = str(diclen_all[elements[1][1]])
            (des_cuffid, des_cufflen, des_ml_mc) = return_des(elements[1][1], dic_ml_m2c)
            (des_Trinityid, des_Trinitylen, des_ml_mt) = return_des(elements[1][1], dic_ml_m2t)
            (des_Oasesid, des_Oaseslen, des_ml_mo) = return_des(elements[1][1], dic_ml_m2o)
            long_description = "ori: maker " + des_makerid + " mlen: " + des_makerlen + \
            " cuff: " + des_cuffid + " cl: " + des_cufflen + " ml_mc: " + des_ml_mc +\
            " Trinity: " + des_Trinityid + " tl: " + des_Trinitylen + " ml_mt: " + des_ml_mt +\
            " Oases: " + des_Oasesid + " ol: " + des_Oaseslen + " ml_mo: " + des_ml_mo
            if elements[1][0] in dicGene_Trinity:
                short_description += ".MT"
            elif elements[1][0] in dicGene_Oases:
                short_description += ".MO"
            elif elements[1][0] in dicGene_maker:
                short_description += ".MM"
            else:
                short_description +=".immpossible"
        elif elements[1][1] in dicGene_Trinity:
            short_description += ".TT"
            long_description = "ori: Trinity " + elements[1][1] + " tlen: " + str(diclen_all[elements[1][1]])
        elif elements[1][1] in dicGene_Oases:
            short_description += ".OO"
            long_description = "ori: Oases " + elements[1][1] + " olen: " + str(diclen_all[elements[1][1]])
        fout_mcot_pr.write(">"+short_description + " " + long_description + " " + des_uni + "\n" +\
        str(dicGene[elements[1][0]].protein_seq.seq) + "\n")
        fout_mcot_cds.write(">"+short_description + " " + long_description + " " + des_uni + "\n" +\
        str(dicGene[elements[1][0]].cds_seq.seq) + "\n")
        fout_mcot_transcript.write(">"+short_description + " " + long_description + " " + des_uni + "\n" +\
        str(dicGene[elements[1][0]].transcript_seq.seq) + "\n")
    fout_mcot_pr.close()
    fout_mcot_cds.close()
    fout_mcot_transcript.close()

