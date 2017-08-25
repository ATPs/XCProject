# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:48:03 2017
data processing for AgSPSPH paper
@author: k
"""

def AgOGS43_44compare():
    '''
    compare AgOGS4.3 and 4.4
    '''
    f_43 = r"D:\Insects\Anopheles_gambiae\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.3.fa"
    f_44 = r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.4.fa"
    import largeFastaFile
    dc43 = largeFastaFile.open_fasta_to_dic(f_43)
    dc44 = largeFastaFile.open_fasta_to_dic(f_44)
    
    
    ls_both = []
    ls_43only = []
    ls_44only = []
    keys = set()
    keys.update(dc43.keys())
    keys.update(dc44.keys())
    for key in keys:
        if key in dc43:
            if key in dc44:
                ls_both.append(key)
            else:
                ls_43only.append(key)
        else:
            if key in dc44:
                ls_44only.append(key)
    
    ls_both_identical = []
    ls_both_nonIdentical = []
    for key in ls_both:
        if str(dc43[key].seq) == str(dc44[key].seq):
            ls_both_identical.append(key)
        else:
            ls_both_nonIdentical.append(key)
    print('''AgP4.3, {len43} proteins
    AgP4.4, {len44} proteins
    
    4.3 only: {len43only}
    4.4 only: {len44only}
    Both: {lenboth43_44}, and {len_both_identical} are identical, {len_both_nonIdentical} are not identical.
    '''.format(len43=len(dc43),len44=len(dc44),len43only = len(ls_43only), len44only=len(ls_44only),lenboth43_44=len(ls_both),len_both_identical = len(ls_both_identical), len_both_nonIdentical = len(ls_both_nonIdentical)))
    
    #20170329 compare AgOGS4.4 and 4.5
    f_45 = r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.5.fa"
    f_44 = r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.4.fa"
    import largeFastaFile
    dc45 = largeFastaFile.open_fasta_to_dic(f_45)
    dc44 = largeFastaFile.open_fasta_to_dic(f_44)
    
    
    ls_both = []
    ls_45only = []
    ls_44only = []
    keys = set()
    keys.update(dc45.keys())
    keys.update(dc44.keys())
    for key in keys:
        if key in dc45:
            if key in dc44:
                ls_both.append(key)
            else:
                ls_45only.append(key)
        else:
            if key in dc44:
                ls_44only.append(key)
    
    ls_both_identical = []
    ls_both_nonIdentical = []
    for key in ls_both:
        if str(dc45[key].seq) == str(dc44[key].seq):
            ls_both_identical.append(key)
        else:
            ls_both_nonIdentical.append(key)
    print('''AgP4.5, {len45} proteins
    AgP4.4, {len44} proteins
    
    4.5 only: {len45only}
    4.4 only: {len44only}
    Both: {lenboth45_44}, and {len_both_identical} are identical, {len_both_nonIdentical} are not identical.
    '''.format(len45=len(dc45),len44=len(dc44),len45only = len(ls_45only), len44only=len(ls_44only),lenboth45_44=len(ls_both),len_both_identical = len(ls_both_identical), len_both_nonIdentical = len(ls_both_nonIdentical)))
    #protein sequences in 4.4 and 4.5 are identical. Check the differet lines
    ls44 = open(f_44).readlines()
    ls45 = open(f_45).readlines()
    diff = []
    for i in range(len(ls44)):
        if ls44[i] != ls45[i]:
            diff.append(ls44[i]+ls45[i])
    print(len(diff))
    

def compare340AgSPSPHproteinsWithOGS44():
    '''
    finally there are 340 proteins of SPSPH from Ag. compare final version of proteins with OGS4.4
    '''
    #read in 340 proteins
    f_340 = r"D:\mine\OneDrive\Lab\works\AgSPSPH\Jiang\340goodProteinSeq.txt"
    import largeFastaFile
    ls340 = largeFastaFile.open_fasta_to_list(f_340)
    f_44 = r"F:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.4.fa"
    dc44 = largeFastaFile.open_fasta_to_dic(f_44)
    dc340to44 = {}
    txt = '''CLIPE8	NA
    SP101	NA
    SP58	AGAP010798-PA
    SP57	NA
    CLIPE14	NA
    CLIPE20	NA
    CLIPA14	AGAP011788-PA
    CLIPA30	NA
    CLIPE30	AGAP010628-PB
    CLIPE31	NA
    SP208	AGAP001365-PA
    SP207	AGAP001366-PA
    CLIPD21	AGAP003686-PA
    CLIPC7	AGAP003689-PA
    CLIPD22	AGAP008996-PA
    SP129	AGAP012777-PA
    SP116	AGAP012817-PA
    SP127	AGAP028710-PA
    SP117	AGAP028711-PA
    SPH106	AGAP028660-PA
    SP144	AGAP028662-PA
    SP109	AGAP028707-PA
    CLIPE33	AGAP012504-PA
    CLIPE9	AGAP028075-PA
    CLIPD9	AGAP008997-PA
    SPH17	AGAP013254-PA
    SPH34	AGAP013164-PA
    CLIPE10	AGAP010545-PA
    CLIPE5	AGAP010546-PA
    SP111	AGAP011427-PA
    SP120	AGAP011429-PA
    SP139	AGAP011430-PA
    SP68	AGAP011431-PA
    CLIPA10	AGAP006954-PA
    CLIPE1	AGAP008091-PA
    CLIPB3b	AGAP013487-PA
    CLIPE11	AGAP028010-PA
    CLIPC4	AGAP000573-PB
    CLIPB7	AGAP002270-PA
    SP213	AGAP005625-PA
    CLIPB9a	AGAP013442-PA
    CLIPB9b	AGAP013442-PD
    SP142	NA
    SP113	NA
    SPH12	NA
    CLIPE32	NA
    SP135	NA
    SPH42	NA
    SPH132	NA
    GPH68	NA
    CLIPC15	NA
    SPH124	NA
    SP44	NA
    SP49	NA
    SPH51	AGAP000182-PA
    CLIPA27	AGAP000290-PA
    SPH216	AGAP000303-PA
    CLIPC6	AGAP000315-PA
    SP125	AGAP000411-PA
    CLIPC5	AGAP000571-PA
    CLIPC10	AGAP000572-PA
    GP20	AGAP001190-PA
    GP13	AGAP001198-PA
    GP12	AGAP001199-PA
    GP44	AGAP001241-PA
    GP63	AGAP001244-PA
    GPH57	AGAP001245-PA
    GP72	AGAP001246-PA
    GP83	AGAP001247-PA
    GP70	AGAP001248-PA
    GP52	AGAP001249-PA
    SP43	AGAP001250-PA
    GP98	AGAP001251-PA
    GP88	AGAP001252-PA
    SP1	AGAP001395-PA
    CLIPD3	AGAP001433-PA
    CLIPB17	AGAP001648-PA
    SP210	AGAP001707-PA
    SPH217	AGAP001798-PA
    SP128	AGAP001924-PA
    CLIPA26	AGAP001964-PA
    CLIPD1	AGAP002422-PA
    GP39	AGAP002432-PA
    SP71	AGAP002543-PA
    CLIPD8	AGAP002784-PA
    CLIPD4	AGAP002811-PA
    CLIPD6	AGAP002813-PA
    CLIPA15	AGAP002815-PA
    SP56	AGAP002842-PA
    CLIPB8	AGAP003057-PA
    CLIPA19	AGAP003245-PA
    CLIPB2	AGAP003246-PA
    CLIPB19	AGAP003247-PA
    SPH6	AGAP003248-PA
    CLIPB3	AGAP003249-PA
    CLIPB4	AGAP003250-PA
    CLIPB1	AGAP003251-PA
    CLIPB6	AGAP003252-PA
    GPH36	AGAP003626-PA
    SP19	AGAP003627-PA
    CLIPE12	AGAP003691-PA
    CLIPE13	AGAP003748-PA
    GP75	AGAP003807-PA
    SP122	AGAP003960-PA
    GPH87	AGAP003971-PA
    SPH52	AGAP003977-PA
    CLIPB5	AGAP004148-PA
    CLIPD10	AGAP004149-PA
    CLIPC11	AGAP004317-PA
    CLIPC3	AGAP004318-PA
    SP26	AGAP004552-PA
    SP25	AGAP004566-PA
    SP18	AGAP004567-PA
    SP24	AGAP004568-PA
    SP64	AGAP004569-PA
    SP22	AGAP004570-PA
    SP21	AGAP004571-PA
    SP47	AGAP004638-PA
    SPH46	AGAP004639-PA
    SPH23	AGAP004644-PA
    GPH3	AGAP004700-PA
    CLIPC9	AGAP004719-PA
    GPH47	AGAP004740-PA
    GP45	AGAP004741-PA
    GP27	AGAP004770-PA
    CLIPB13	AGAP004855-PA
    CLIPD11	AGAP004858-PA
    SPH5	AGAP004859-PA
    GP1	AGAP004900-PA
    GP17	AGAP005065-PA
    SP219	AGAP005072-PA
    SP3	AGAP005194-PA
    SP4	AGAP005195-PA
    SP2	AGAP005196-PA
    SP70	AGAP005303-PA
    SP141	AGAP005304-PA
    GPH18	AGAP005310-PA
    GPH35	AGAP005587-PA
    GPH38	AGAP005591-PA
    GPH37	AGAP005592-PA
    SPH38	AGAP005593-PA
    SPH65	AGAP005594-PA
    GPH28	AGAP005596-PA
    GPH29	AGAP005597-PA
    SPH60	AGAP005598-PA
    GPH46	AGAP005642-PA
    GP71	AGAP005663-PA
    GP73	AGAP005664-PA
    GP82	AGAP005665-PA
    GP58	AGAP005669-PA
    GP49	AGAP005670-PA
    GP65	AGAP005671-PA
    SP31	AGAP005686-PA
    SP33	AGAP005687-PA
    SP30	AGAP005688-PA
    GP43	AGAP005689-PA
    GP42	AGAP005690-PA
    GP41	AGAP005691-PA
    GPH91	AGAP005702-PA
    GPH61	AGAP005703-PA
    GPH89	AGAP005704-PA
    GPH94	AGAP005705-PA
    GPH95	AGAP005706-PA
    GPH77	AGAP005707-PA
    GPH69	AGAP005708-PA
    GPH90	AGAP005709-PA
    SPH79	AGAP005790-PA
    SPH80	AGAP005791-PA
    SPH78	AGAP005792-PA
    SPH8	AGAP005793-PA
    SP73	AGAP006087-PA
    SP218	AGAP006120-PA
    SP54	AGAP006192-PA
    GP19	AGAP006385-PA
    GP5	AGAP006416-PA
    GPH55	AGAP006485-PA
    GPH53	AGAP006486-PA
    GPH54	AGAP006487-PA
    GPH84	AGAP006488-PA
    GPH78	AGAP006489-PA
    SP63	AGAP006539-PA
    SP214	AGAP006631-PA
    GP85	AGAP006672-PA
    GP81	AGAP006673-PA
    GP66	AGAP006674-PA
    GP60	AGAP006675-PA
    GPH50	AGAP006676-PA
    GPH67	AGAP006677-PA
    GP15	AGAP006707-PA
    GP22	AGAP006709-PA
    GP23	AGAP006710-PA
    GP21	AGAP006711-PA
    GP56	AGAP006869-PA
    SP112	AGAP007043-PB
    SP62	AGAP007141-PA
    GP25	AGAP007142-PA
    GPH16	AGAP007165-PA
    GP86	AGAP007251-PA
    GP51	AGAP007252-PA
    GP40	AGAP007262-PA
    SP212	AGAP007280-PA
    SP138	AGAP007795-PA
    CLIPD2	AGAP008183-PA
    SP75	AGAP008276-PA
    SP77	AGAP008277-PA
    GP6	AGAP008290-PA
    GP7	AGAP008291-PA
    GP10	AGAP008293-PA
    GP9	AGAP008294-PA
    GP24	AGAP008295-PA
    GP26	AGAP008296-PA
    CLIPE15	AGAP008403-PA
    SP137	AGAP008558-PA
    SP61	AGAP008649-PA
    SP130	AGAP008808-PA
    CLIPC1	AGAP008835-PA
    GP11	AGAP008861-PA
    SP143	AGAP008891-PA
    SPH29	AGAP008911-PA
    CLIPD12	AGAP008995-PA
    CLIPD7	AGAP008998-PA
    CLIPD13	AGAP009000-PA
    CLIPD14	AGAP009006-PA
    GP92	AGAP009121-PA
    GPH93	AGAP009122-PA
    CLIPD15	AGAP009211-PA
    CLIPB11	AGAP009214-PA
    CLIPB18	AGAP009215-PA
    CLIPD16	AGAP009216-PA
    CLIPB12	AGAP009217-PA
    SP67	AGAP009218-PA
    SP66	AGAP009219-PA
    CLIPD17	AGAP009220-PA
    SP136	AGAP009249-PA
    CLIPE16	AGAP009251-PB
    CLIPE17	AGAP009252-PA
    CLIPB16	AGAP009263-PA
    CLIPE18	AGAP009273-PA
    GPH4	AGAP009828-PA
    CLIPB15	AGAP009844-PA
    CLIPD18	AGAP009849-PA
    SP72	AGAP009966-PA
    SPH10	AGAP010015-PA
    GP97	AGAP010240-PA
    GP14	AGAP010243-PA
    SP105	AGAP010415-PA
    CLIPE4	AGAP010530-PA
    GP32	AGAP010614-PA
    GP33	AGAP010615-PA
    SPH14	AGAP010617-PA
    SP20	AGAP010618-PA
    SPH36	AGAP010619-PA
    SPH35	AGAP010620-PA
    SP114	AGAP010635-PA
    SP133	AGAP010659-PA
    GPH74	AGAP010661-PA
    SP126	AGAP010662-PA
    GPH79	AGAP010663-PA
    CLIPA28	AGAP010730-PA
    CLIPA8	AGAP010731-PA
    CLIPB14	AGAP010833-PA
    CLIPA9	AGAP010968-PA
    CLIPE19	AGAP011040-PA
    CLIPD19	AGAP011325-PA
    SPH27	AGAP011432-PA
    SP123	AGAP011432-PA
    GP96	AGAP011477-PA
    GP2	AGAP011590-PA
    GPH48	AGAP011608-PA
    SP140	AGAP011669-PA
    CLIPE21	AGAP011719-PA
    CLIPA4	AGAP011780-PA
    CLIPA12	AGAP011781-PA
    CLIPE2	AGAP011782-PA
    CLIPA13	AGAP011783-PA
    CLIPE6	AGAP011785-PA
    CLIPE7	AGAP011786-PA
    CLIPA5	AGAP011787-PA
    CLIPA6	AGAP011789-PA
    CLIPA2	AGAP011790-PB
    CLIPA1	AGAP011791-PA
    CLIPA7	AGAP011792-PA
    CLIPA31	AGAP011793-PA
    CLIPA32	AGAP011794-PA
    SP205	AGAP011908-PA
    SP201	AGAP011909-PA
    SP202	AGAP011910-PA
    SP206	AGAP011912-PA
    SP203	AGAP011913-PA
    SP204	AGAP011914-PA
    GPH64	AGAP011917-PA
    GP62	AGAP011918-PA
    GPH76	AGAP011919-PA
    GP59	AGAP011920-PA
    CLIPE22	AGAP012020-PA
    CLIPE23	AGAP012021-PA
    CLIPE24	AGAP012022-PA
    CLIPC12	AGAP012034-PA
    CLIPB20	AGAP012037-PA
    SP81	AGAP012315-PA
    SP45	AGAP012328-PA
    SPH37	AGAP012470-PA
    SP40	AGAP012473-PA
    SP13	AGAP012492-PA
    CLIPE25	AGAP012502-PA
    SPH15	AGAP012566-PA
    CLIPA3	AGAP012591-PA
    GP31	AGAP012670-PA
    SP41	AGAP012671-PA
    GP80	AGAP012692-PA
    SP118	AGAP012749-PA
    SP11	AGAP012778-PA
    GP8	AGAP012842-PA
    SP53	AGAP012946-PA
    SP209	AGAP013020-PA
    CLIPD20	AGAP013089-PA
    SPH48	AGAP013117-PA
    CLIPB36	AGAP013184-PA
    SP7	AGAP013221-PA
    SP134	AGAP013252-PA
    SP28	AGAP013452-PA
    SP74	AGAP013716-PA
    CLIPC13	AGAP028007-PA
    GPH99	AGAP028031-PA
    CLIPE26	AGAP028069-PA
    SP119	AGAP028071-PA
    CLIPE34	AGAP028102-PA
    SPH50	AGAP028136-PA
    CLIPC14	AGAP028167-PA
    SP55	AGAP028179-PA
    CLIPE27	AGAP028183-PA
    SPH16	AGAP028216-PA
    SP76	AGAP028218-PA
    CLIPE28	AGAP028229-PA
    CLIPE29	AGAP028641-PA
    SP108	AGAP028661-PA
    SP39	AGAP028704-PA
    GP30	AGAP028705-PA
    SPH107	AGAP028706-PA
    SP115	AGAP012270-PA'''
    for line in txt.split('\n'):
        _sp,_ogs = line.split()
        dc340to44[_sp] = _ogs
    
    comments = {}
    fout = open('comment.txt','w')
    for s in ls340:
        _sp = s.id
        if dc340to44[_sp] == 'NA':
            comments[_sp] = 'no AgID'
            fout.write(_sp+'\t'+comments[_sp]+'\n')
        else:
            _ogs = dc340to44[_sp]
            seqsp = str(s.seq)
            seqsp = seqsp.replace('*','')
            if _ogs in dc44:
                seqogs = str(dc44[_ogs].seq)
                if seqsp == seqogs:
                    comments[_sp] = 'same'
                    fout.write(_sp+'\t'+comments[_sp]+'\n')
                elif seqsp[:len(seqogs)] == seqogs:
                    comments[_sp] = 'C_incomplete'
                    fout.write(_sp+'\t'+comments[_sp]+'\n')
                elif seqsp[len(seqsp)-len(seqogs):] == seqogs:
                    comments[_sp] = 'N_incomplete'
                    fout.write(_sp+'\t'+comments[_sp]+'\n')
                elif seqogs in seqsp:
                    comments[_sp] = 'NC_incomplete'
                    fout.write(_sp+'\t'+comments[_sp]+'\n')
                else:
                    comments[_sp] = 'other'
                    fout.write(_sp+'\t'+comments[_sp]+'\n')
            else:
                comments[_sp] = 'OGS not in 4.4'
                fout.write(_sp+'\t'+comments[_sp]+'\n')
    fout.close()
    
    ogsIDs = list(dc340to44.values())
    for _ogs in ogsIDs:
        if _ogs != 'NA' and _ogs not in ls_both:
            print(_ogs)

def getGenomeLocationAgSPSPH():
    '''
    get genome location of AgSPSPHs, with OGS P4.4 ids
    '''
    f_ogsIDs = r"D:\mine\OneDrive\Lab\works\AgSPSPH\GenomeLocationOfSPSPHAg340\AgIDs.txt"
    ls_ogsIDs = open(f_ogsIDs).read().split()
    ls_ogsIDs = [ele.split('-')[0] for ele in ls_ogsIDs]
    
    dc_geneLocation = {}
    f_gff3 = r'F:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3'
    for ele in open(f_gff3):
        if '\tgene\t' in ele:
            eles = ele.split('\t')
            e_chr = eles[0]
            e_start = eles[3]
            e_end = eles[4]
            e_strand = eles[6]
            e_gene = eles[-1].split('=')[1].split(';')[0]
            dc_geneLocation[e_gene] = [e_chr, e_start, e_end, e_strand]
    
    fout = open(r"D:\mine\OneDrive\Lab\works\AgSPSPH\GenomeLocationOfSPSPHAg340\AgIDsLocation.txt",'w')
    for ele in ls_ogsIDs:
        fout.write(ele + '\t'+'\t'.join(dc_geneLocation[ele])+'\n')
    fout.close()
    
    f_ogsIDs = r"D:\mine\OneDrive\Lab\works\AgSPSPH\GenomeLocationOfSPSPHAg340\AgIDs.txt"
    ls_mRNAIDs = open(f_ogsIDs).read().split()
    ls_mRNAIDs = [ele.replace('-P','-R') for ele in ls_mRNAIDs]
    
    dc_mRNALocation = {}
    f_gff3 = r'F:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.4.gff3'
    for ele in open(f_gff3):
        if '\tmRNA\t' in ele:
            eles = ele.split('\t')
            e_chr = eles[0]
            e_start = eles[3]
            e_end = eles[4]
            e_strand = eles[6]
            e_mRNA = eles[-1].split('=')[1].split(';')[0]
            e_exon = 0
            dc_mRNALocation[e_mRNA] = [e_chr, e_start, e_end, e_strand, e_exon]
        if '\texon\t' in ele:
            eles = ele.split('\t')
            e_mRNA = eles[-1].split('=')[1].split(';')[0]
            if e_mRNA in dc_mRNALocation:
                dc_mRNALocation[e_mRNA][-1] += 1
    
    fout = open(r"D:\mine\OneDrive\Lab\works\AgSPSPH\GenomeLocationOfSPSPHAg340\AgIDsmRNALocation.txt",'w')
    for ele in ls_mRNAIDs:
        fout.write(ele + '\t'+'\t'.join(str(e) for e in dc_mRNALocation[ele])+'\n')
    fout.close()

def findSPSPHcdsFromMCOTidsOrOtherInfo20170302():
    '''
    with the 340 good protein sequences, try to find the corresponding cds sequences.
    '''
    #read in the sequences
    import largeFastaFile
    ls_SPs = largeFastaFile.open_fasta_to_list(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\Jiang\340goodProteinSeq.txt")
    #remove star symbol,convert seq to str
    for _s in ls_SPs:
        if _s.seq[-1] == '*':
            _s.seq = str(_s.seq[:-1])
        else:
            _s.seq = str(_s.seq)
    
    # read in mcot sequences
    ls_mcot = largeFastaFile.open_fasta_to_list(r"D:\Insects\Anopheles_gambiae\mcot\20160316AgMcotproteinNoStar.fa")
    for _s in ls_mcot:
        if _s.seq[-1] == '*':
            _s.seq = str(_s.seq[:-1])
        else:
            _s.seq = str(_s.seq)
    
    #use seq as key to get identical sequences in SPSPH and mcot
    dcSeq2idSP = {}
    for _s in ls_SPs:
        dcSeq2idSP[_s.seq] = _s.id
    dcSeq2idmcot = {}
    for _s in ls_mcot:
        dcSeq2idmcot[_s.seq] = _s.id
    ls_SPMCOTpair = []
    ls_SPnotfound = []
    for _e in dcSeq2idSP:
        if _e in dcSeq2idmcot:
            ls_SPMCOTpair.append((dcSeq2idSP[_e],dcSeq2idmcot[_e]))
        else:
            ls_SPnotfound.append(dcSeq2idSP[_e])
    print('identical: ',len(ls_SPMCOTpair),'not found yet:', len(ls_SPnotfound))
    
    #read in AgOGS4.4
    ls_ogs = largeFastaFile.open_fasta_to_list(r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_PEPTIDES_AgamP4.4.fa")
    for _s in ls_ogs:
        if _s.seq[-1] == '*':
            _s.seq = str(_s.seq[:-1])
        else:
            _s.seq = str(_s.seq)
    dcSeq2idogs = {}
    for _s in ls_ogs:
        dcSeq2idogs[_s.seq] = _s.id
    ls_SPogsPair = []
    ls_SPogsnotfound = []
    for _e in dcSeq2idSP:
        if _e in dcSeq2idogs:
            ls_SPogsPair.append((dcSeq2idSP[_e],dcSeq2idogs[_e]))
        else:
            ls_SPogsnotfound.append(dcSeq2idSP[_e])
    print('identical: ',len(ls_SPogsPair),'not found yet:', len(ls_SPogsnotfound))
    
    ls_SPnotIdenticalOGSmcot = [_e for _e in ls_SPnotfound if _e in ls_SPogsnotfound]
    print(len(ls_SPnotIdenticalOGSmcot))
    
    fout = open('list1.txt','w')
    for _e in ls_SPMCOTpair:
        fout.write('%s\t%s\n'%(_e[0],_e[1]))
    fout.close()
    fout = open('list2.txt','w')
    for _e in ls_SPogsPair:
        fout.write('%s\t%s\n'%(_e[0],_e[1]))
    fout.close()
    
    #get a dc with mcot and ogs4.4 cds or transcript
    dc_cds = {}
    dc_cds.update(largeFastaFile.open_fasta_to_dic(r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.4.fa"))
    dc_cds.update(largeFastaFile.open_fasta_to_dic(r"D:\Insects\Anopheles_gambiae\mcot\20160316AgMcotcds.fa"))
    
    ls_cdsid = open('list.txt').readlines()
    fout = open('cds.txt','w')
    for ele in ls_cdsid:
        _gene, _key = ele.split()
        if 'AGAP' in _key:
            _key = _key.replace('-P','-R')
        fout.write('>'+_gene+'\n'+str(dc_cds[_key].seq) +'\n')
    fout.close()
    
    #make recombinant cds, remove ids in list.txt
    ls_toRemove = open('list.txt').read().split()
    ls_toRemove = set(ele.replace('-P','-R') for ele in ls_toRemove)
    
    dc_ogs = largeFastaFile.open_fasta_to_dic(r"D:\Insects\Anopheles_gambiae\2017OGS\Anopheles-gambiae-PEST_TRANSCRIPTS_AgamP4.4.fa")
    dc_SPSPH = largeFastaFile.open_fasta_to_dic(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\Jiang\340goodProteinCDS.txt")
    fout = open('AgamP4.4RecombinantCDS.fa','w')
    for ele in dc_SPSPH:
        fout.write('>'+dc_SPSPH[ele].id+'\n'+str(dc_SPSPH[ele].seq)+'\n')
    for ele in dc_ogs:
        if ele not in ls_toRemove:
            fout.write('>'+dc_ogs[ele].id+'\n'+str(dc_ogs[ele].seq)+'\n')
    fout.close()

def getFPKMvaluesForAGAP008292_20170403():
    '''
    it's hard to read my previous code about how to get FPKM values for AgSPSPHs.
    decide to write new function
    '''
    ls_srr = '

def plotFPKMfigs20170306():
    '''
    '''
    import HeatmapCluster
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170305AgSPSPH_RSEM_FPKM_expression340SimplifiedLibs.csv", index_col = 0)
    print(df.shape)
    df.head(2)
    dfGP = df.loc[[ele for ele in df.index if ele[:2] == 'GP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfGP,'20170306AgGP_')
    
    dfSP = df.loc[[ele for ele in df.index if ele[:2] == 'SP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfSP,'20170306AgSP_')
    
    dfCLIP = df.loc[[ele for ele in df.index if ele[:2] == 'CL'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfCLIP,'20170306AgCLIP_')

def plotFPKMfigs20170314():
    '''
    '''
    import HeatmapCluster
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170314AgSPSPH_RSEM_FPKM_expression335SimplifiedLibs.csv", index_col = 0)
    print(df.shape)
    df.head(2)
    dfGP = df.loc[[ele for ele in df.index if ele[:2] == 'GP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfGP,'20170314AgGP_')
    
    dfSP = df.loc[[ele for ele in df.index if ele[:2] == 'SP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfSP,'20170314AgSP_')
    
    dfCLIP = df.loc[[ele for ele in df.index if ele[:2] == 'CL'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfCLIP,'20170314AgCLIP_')

def plotFPKMfigs20170320():
    '''
    correct one error. keep the order of result in 20170314
    '''
    import HeatmapCluster
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170320AgSPSPH_RSEM_FPKM_expression335SimplifiedLibs.csv", index_col = 0)
    print(df.shape)
    df.head(2)
    dfGP = df.loc[[ele for ele in df.index if ele[:2] == 'GP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfGP,'20170320AgGP_',col_cluster=False, row_cluster=False)
    
    dfSP = df.loc[[ele for ele in df.index if ele[:2] == 'SP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfSP,'20170320AgSP_',col_cluster=False, row_cluster=False)
    
    dfCLIP = df.loc[[ele for ele in df.index if ele[:2] == 'CL'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfCLIP,'20170320AgCLIP_',col_cluster=False, row_cluster=False)
    
    #replot GP figure
    import HeatmapCluster
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170320AgSPSPH_RSEM_FPKM_expression335SimplifiedLibs.csv", index_col = 0)
    print(df.shape)
    df.head(2)
    dfGP = df.loc[[ele for ele in df.index if ele[:2] == 'GP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfGP,'20170404AgGP_',col_cluster=False, row_cluster=True)
    
    #plotSerpinManduca
    import HeatmapCluster
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\20121023_manduca_genome_new\PaperWriting\2017Serpin\2017RSEM_Serpin_ForFigureFPKM.csv", index_col = 0)
    print(df.shape)
    df.head(2)
    HeatmapCluster.fpkmHierarchicalCluster(df,'20170404MsSerpin_',col_cluster=False, row_cluster=True)
    
    #plot SP 20170430
    import HeatmapCluster
    import pandas as pd
    df = pd.read_csv(r"C:\Users\k\OneDrive\Lab\works\AgSPSPH\RSEM_expression\20170320AgSPSPH_RSEM_FPKM_expression335SimplifiedLibs.csv", index_col = 0)
    print(df.shape)
    dfSP = df.loc[[ele for ele in df.index if ele[:2] == 'SP'],:]
    HeatmapCluster.fpkmHierarchicalCluster(dfSP,'20170430AgSP_',col_cluster=False, row_cluster=True)


def compareOriginalSPSPFileWithExcelCompleteList20170327():
    '''
    some proteins may be missing from TableS1 complete sheet
    '''
    #read in original file
    txt = open('temp.txt').read()
    import re
    lsori1 = re.findall('AGAP\d\d\d\d\d\d-\w\w',txt)
    lsori2 = re.findall('MCOT\d*\.\d*\.\w\w',txt)
    lsori = lsori1 + lsori2
    
    #read in names in tableS1 complete
    tablenames = open('list.txt').read()
    
    for ele in lsori:
        if ele not in tablenames:
            print(ele)