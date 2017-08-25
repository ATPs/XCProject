from Bio import SeqIO

import sys
sys.path.append("D:\\P\\3Language\\Xiaolong\\python\\XCProject\\MCOT")

#from XCMCOT import *
#from NucleotideSequence20140707 import *


"""
deal with blast result. First develop those simple ones deal with data in the outfmt 6, 
e.g. "TCONS_00002781|m.70	TR:T1I2R8_RHOPR	46.74	92	48	1	1	92	130	220	6e-22	94.0"
other functions will be added later.
"""

class Blast6(object):
    """
    to outfmt 6 of blast result in a class. 
    elements inculde "query_id, subject_id, identity, alignment_length, mismatches, gap_opens, q_start, q_end, \
    s_start, s_end, e_value, bit_score"
    """
    def __init__(self, query_id, subject_id, identity, alignment_length, mismatches, gap_opens, q_start, q_end, \
    s_start, s_end, e_value, bit_score, query_length = None, subject_length = None):
        self.query_id = query_id
        self.subject_id = subject_id
        self.identity = identity
        self.alignment_length = alignment_length
        self.mismatches = mismatches
        self.gap_opens = gap_opens
        self.q_start = q_start
        self.q_end = q_end
        self.s_start = s_start
        self.s_end = s_end
        self.e_value = e_value
        self.bit_score = bit_score
        self.query_length = query_length
        self.subject_length = subject_length
    def details(self):
        """
        return a line with all information of this blast result line
        the format is the same as in outfmt 6
        also print out this result
        """
        line_details = str(self.query_id)\
        + "\t" + str(self.subject_id) \
        + "\t" + str(self.identity) \
        + "\t" + str(self.alignment_length) \
        + "\t" + str(self.mismatches) \
        + "\t" + str(self.gap_opens) \
        + "\t" + str(self.q_start) \
        + "\t" + str(self.q_end) \
        + "\t" + str(self.s_start) \
        + "\t" + str(self.s_end) \
        + "\t" + str(self.e_value) \
        + "\t" + str(self.bit_score)
        return line_details
    def print_details(self):
        """
        return a line with all information of this blast result line
        the format is the same as in outfmt 6
        also print out this result
        """
        line_details = self.details()
        print("query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")
        print(line_details)
        return line_details
    def __str__(self):
        """
        print the result in orignal format
        """
        return self.details()





def file_to_Blast6(file_name, dic_XCGene_query = None, dic_XCGene_subject = None):
    """
    given file name including the path. open it and return a list, with each element in Blast6 format
    """
    lines = open(file_name).readlines()
    list_Blast6=[]
    for line in lines:
        element_list=line.split()
        element_Blast6 = Blast6(element_list[0],element_list[1],float(element_list[2]),int(element_list[3]),\
        int(element_list[4]),int(element_list[5]),int(element_list[6]),int(element_list[7]),\
        int(element_list[8]),int(element_list[9]),float(element_list[10]),float(element_list[11]))
        if dic_XCGene_query != None:
            element_Blast6.query_length = dic_XCGene_query[element_Blast6.query_id].protein_length()
        if dic_XCGene_subject != None:
            element_Blast6.subject_length = dic_XCGene_subject[element_Blast6.subject_id].protein_length()
        list_Blast6.append(element_Blast6)
    return list_Blast6

def filter_identity_Blast6(list_Blast6, min_identity=0, max_identity=100):
    """
    given a list of element in Blast6 format, return a list with identity greater or equal to the min_identity, 
    less or equal to the max_identity
    """
    new_list_Blast6=[]
    for element in list_Blast6:
        if element.identity >= min_identity and element.identity <= max_identity:
            new_list_Blast6.append(element)
    return new_list_Blast6


def filter_alignment_length_Blast6(list_Blast6, min_alignment_length=0, max_alignment_length=10e10):
    """
    given a list of element in Blast6 format, return a list with alignment length greater or equal to the min_alignment_length, 
    less or equal to the max_alignment_length
    default setting: min_alignment_length=0, max_alignment_length=10e10
    """
    new_list_Blast6=[]
    for element in list_Blast6:
        if element.alignment_length >= min_alignment_length and element.alignment_length <= max_alignment_length:
            new_list_Blast6.append(element)
    return new_list_Blast6


def filter_e_value_Blast6(list_Blast6, min_e_value=0, max_e_value=100):
    """
    given a list of element in Blast6 format, return a list with alignment length greater or equal to the min_e_value, 
    less or equal to the max_e_value
    default setting: min_e_value=0, max_e_value=100
    """
    new_list_Blast6=[]
    for element in list_Blast6:
        if element.e_value >= min_e_value and element.e_value <= max_e_value:
            new_list_Blast6.append(element)
    return new_list_Blast6


def filter_bit_score_Blast6(list_Blast6, min_bit_score=0,max_bit_score=10e10):
    """
    given a list of element in Blast6 format, return a list with alignment length greater or equal to the min_bit_score, 
    less or equal to the max_bit_score
    default setting: min_bit_score=0,max_bit_score=10e10
    """
    new_list_Blast6=[]
    for element in list_Blast6:
        if element.bit_score >= min_bit_score and element.bit_score <= max_bit_score:
            new_list_Blast6.append(element)
    return new_list_Blast6

def best_match_list_Blast6(list_Blast6, simple=1):
    """
    given a list of elements in Blast6 format, return a list.
    for each query_id, only return those with subject_id. 
    Under default setting, the first suject_id have best match with the query_id
    """
    if simple != 1:
        print("not a simple case. Function will be developed later!")
        return None
    list_Blast6_best=[]
    query_list=[]
    for element in list_Blast6:
        if element.query_id not in query_list:
            query_list.append(element.query_id)
            list_Blast6_best.append(element)
        elif element.query_id == list_Blast6_best[-1].query_id and element.subject_id == list_Blast6_best[-1].subject_id:
            list_Blast6_best.append(element)
    return list_Blast6_best
        


def list_Blast6_to_dic(list_Blast6):
    """
    return a dictonary, with a tuple of (query_id, subject_id) as key, and all members in the list with the same query_id
    and subject_id into a list as the value.
    >>>myBlast6=[element1, element2, ...] 
    >>>list_Blast6_to_dic(myBlast6)
    {(query_id, subject_id):[element1, element...]; ...}
    """
    if list_Blast6 ==[]:
        print("empty list")
        return None
    dic_Blast6={}
    for element in list_Blast6:
        dic_key=(element.query_id, element.subject_id)
        if dic_key not in dic_Blast6:
            dic_Blast6[dic_key]=[]
        dic_Blast6[dic_key].append(element)
    return dic_Blast6

def sum_match_length_list_same_query_subject(list_Blast6):
    """
    given a list in of Blast6 class with the same query_id and subject_id, calculate the matched length.
    Here the match length is only based on the query sequence.
    query from (30, 60) and (70,90), the matched length should be 52
    query from (30, 90) and (60,70), the matched length should be 61
    query from (30, 70) and (60,90), the matched length should be 61
    """
    if len(list_Blast6) > 1:
        for i in range(1, len(list_Blast6)):
            if list_Blast6[0].query_id != list_Blast6[i].query_id or list_Blast6[0].subject_id != list_Blast6[i].subject_id:
                print("Error! not the same query and subject!")
                return None
    matched_sites=set()
    for element_Blast6 in list_Blast6:
        matched_sites.update(range(element_Blast6.q_start,element_Blast6.q_end+1))
    matched_sites2=set()
    for element_Blast6 in list_Blast6:
        matched_sites2.update(range(element_Blast6.s_start,element_Blast6.s_end+1))
    return min(len(matched_sites),len(matched_sites2))
#    return len(matched_sites)


def sum_match_length_list(list_Blast6):
    """
    given a list of Blast6 class, return a dictionary, with a tuple of (query_id, subject_id) as key,
    and the matched length based on the query as the value.
    >>>myBlast6=[element1, element2, ...] 
    >>>list_Blast6_to_dic(myBlast6)
    {(query_id, subject_id): matched_length; ...}
    """
    dic_Blast6=list_Blast6_to_dic(list_Blast6)
    dic_matched_length={}
    for key in dic_Blast6:
        dic_matched_length[key] = sum_match_length_list_same_query_subject(dic_Blast6[key])
    return dic_matched_length

def file_Blast6_to_best_match_length_list(file_name, min_identity = 95):
    """
    given a file name, return a dictionary, with a tuple of (query_id, subject_id) as key,
    and the matched length based on the query as the value.
    filter identity >95% first
    keep the only those with best subject_id
    then calculate the matched_length
    >>>myBlast6=[element1, element2, ...] 
    >>>list_Blast6_to_dic(myBlast6)
    {(query_id, subject_id): matched_length; ...}
    """
    list_Blast6 = file_to_Blast6(file_name)
    list_Blast6_quality_filter = filter_identity_Blast6(list_Blast6, min_identity)
    list_Blast6_quality_filter_best = best_match_list_Blast6(list_Blast6_quality_filter)
    return sum_match_length_list(list_Blast6_quality_filter_best)

def file_Blast6_to_best_match_length_list_refine(file_name, file_name_reference):
    """
    given a file name and a reference file name. very similar to function file_Blast6_to_best_match_length_list.
    use file_name reference to help improve the matched length result. keep the smaller matched length in file_name or file_name_reference
    for example, maker model blast to cufflinks in file_name, cufflinks blast to maker in file_name_reference.
    based on maker model, matched length is 500, while based on cufflinks, matched length is 400. 
    this difference is because some part of maker model matched to cufflinks more than one time. use 400 as the matched length.
    return a dictionary, with a tuple of (query_id, subject_id) as key,
    and the matched length based on the query as the value.
    """
    match_length_list = file_Blast6_to_best_match_length_list(file_name)
    match_length_list_reference = file_Blast6_to_best_match_length_list(file_name_reference)
    for element_key in match_length_list:
        if (element_key[1],element_key[0]) in match_length_list_reference:
            match_length_list[element_key] = min(match_length_list[element_key], match_length_list_reference[(element_key[1],element_key[0])])
    return match_length_list

def dic_keyoftupple_trans(dic):
    """
    dic like {(cuff_id, trinity_id): matched_length; ...}
    change to {cuff_id: [trinity_id, matched_length]; ...}
    """
    dicnew = {}
    for elements in dic:
        dicnew[elements[0]] = [elements[1], dic[elements]]
    return dicnew


def file_Blast6_to_match_length_list(file_name, min_identity = 35):
    """
    given a file name, return a dictionary, with a tuple of (query_id, subject_id) as key,
    and the matched length based on the query as the value.
    >>>myBlast6=[element1, element2, ...] 
    >>>list_Blast6_to_dic(myBlast6)
    {(query_id, subject_id): matched_length; ...}

    """
    list_Blast6 = file_to_Blast6(file_name)
    list_Blast6_quality_filter = filter_identity_Blast6(list_Blast6, min_identity)
    return sum_match_length_list(list_Blast6_quality_filter)


def file_Blast6_to_Uniprot_matching_score(file_Blast6, diclen_fasta_query, diclen_fasta_subject):
    """
    given a file of Blast6 format Here, it is some query sequence blast against uniprot seq
    return the matching score. 
    for example, query sequence length 100, and matched length and subject length can be
    (100, 1000), (50, 50), matching score is defined as ml*ml/(ql*sl)
    here, the best matching score, we choose the 100*100/(100*1000) = 0.1, not the 50*50/(100*50) = 0.5 one
    return dictionary like {query_id: (subject_id, matched length, matching score, query_length, subject_length)}
    """
    dic_matches = file_Blast6_to_match_length_list(file_Blast6)
    dic_query_matches = {}
    for elements in dic_matches:
        dic_query_matches[elements[0]]=[]
    for elements in dic_matches:
        dic_query_matches[elements[0]].append((elements[1], \
        dic_matches[elements],dic_matches[elements] ** 2 /(diclen_fasta_query[elements[0]] * diclen_fasta_subject[elements[1]]), \
        diclen_fasta_query[elements[0]], diclen_fasta_subject[elements[1]]))
    for elements in dic_query_matches:
        dic_query_matches[elements] = max(dic_query_matches[elements], key=lambda x:x[1])
    return dic_query_matches

#20160314 use normal scoring matrix, get the output in blast6 format, and calculate matched length using global alignment
def mcotMLwithblast6(dcFa1,dcFa2,file_blast6):
    """
    dcFa1, dcFa2 is two dictionary of fasta sequences.  file_blast6 is blast6 output of blast dcFa1 to dcFa2.
    return a dictionary, with key in dcFa1 as key, best match in dcFa2 and querylen, subjectlen, their matched length as value.
    calculate based on file_blast6, so some sequence in dcFa1 may not exist in the output
    * symbol in the seq is not counted in sequence length
    {"seq1":("seq2",100,100,100),...}
    """
    lsBlast6 = open(file_blast6).readlines()
    for num in range(len(lsBlast6)):
        lsBlast6[num] = lsBlast6[num].split()
    import XCBlastpLike
    dcML ={}
    count = 0
    import time
    time0 = time.time()
    for blast6 in lsBlast6:
        count += 1
        if count % 100 == 0:
            print(count, "finished, time used: ", time.time() - time0)
#            break
        if float(blast6[2]) > 70:
            query = blast6[0]
            subject = blast6[1]
            query_seq = dcFa1[query].seq
            if query_seq[-1] == "*":
                query_seq = query_seq[:-1]
            subject_seq = dcFa2[subject].seq
            if subject_seq[-1] == "*":
                subject_seq = subject_seq[:-1]
            querylen = len(query_seq)
            subjectlen = len(subject_seq)
            
            if (int(blast6[3]) - int(blast6[4])) / querylen > 0.6:
                qeury_aln, subject_aln = XCBlastpLike.proteinPairwiseAlignGlobal(query_seq,subject_seq)
                ml = XCBlastpLike.proteinAlignLength(qeury_aln, subject_aln)
                if query not in dcML:
                    dcML[query] = (subject,querylen,subjectlen,ml)
                else:
                    if ml > dcML[query][3]:
                        dcML[query] = (subject,querylen,subjectlen,ml)
    return dcML


     