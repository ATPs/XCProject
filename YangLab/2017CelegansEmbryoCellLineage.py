# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 19:51:26 2017

@author: ATPs
"""

def checkfiles20170928():
    '''
    test files from Prof. Yang.
    
    '''
    import pandas as pd
    folder = r'C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans'+'\\'
    f_embryo = folder + 'CelegansEmbryonicCellLineage.alm'
    df_embryo = pd.read_csv(f_embryo,sep='\t',dtype =str)
    df_all = pd.read_csv(folder + 'C. elegans Cell List - WormAtlas.tsv',sep='\t')
    '''df_all looks like
            Cell   Lineage Name             Description
        0    AB           P0 a  Embryonic founder cell
        1  ADAL  AB plapaaaapp       Ring interneurons
        2  ADAR  AB prapaaaapp       Ring interneurons
    '''
    #check if all Name in df_embryo in df_all
    cellnamesAll = set(df_all.loc[:,'Cell'])
    cellnamesEmb = set(df_embryo.loc[:,'Name'])
    nameOnlyInEmb = [e for e in cellnamesEmb if e not in cellnamesAll]
    print(len(nameOnlyInEmb))
    #count "Cell" in df_all if Cell is more than one
    from collections import Counter
    cellMoreThan1 = [e for e in Counter(df_all.loc[:,'Cell']).most_common() if e[1]>1]
    #some "Cell" have more than 1 lineage name, find them
    df_complexLineage = df_all[df_all.loc[:,'Lineage Name'].apply(lambda x: ',' in x or '/' in x)]
    #duplicated "Cell Lineages"
    df_duplicateLineage = df_all[df_all.loc[:,'Lineage Name'].apply(lambda x: list(df_all.loc[:,'Lineage Name']).count(x)>1)]
    #find living embryonic cells
    print(len(cellMoreThan1), df_complexLineage.shape, df_duplicateLineage.shape)
    import re
    df_livingEmbryo = df_all[df_all.loc[:,'Lineage Name'].apply(lambda x: bool(re.match('\w+ \w+',x)))]
    df_livingPostEmbryoNormal = df_all[df_all.loc[:,'Lineage Name'].apply(lambda x: bool(re.match('^\w+[\.]\w+$',x)))]
    df_abnormalCell = df_all[df_all.loc[:,'Lineage Name'].apply(lambda x: not bool(re.match('^\w+[ \.]\w+$',x)))]
    print(df_livingEmbryo.shape, df_livingPostEmbryoNormal.shape, df_abnormalCell.shape)
    #find headers of normal cell (in the form of 'XX.aabb' or 'XX aabb')
    df_normalCell = df_all[df_all.loc[:,'Lineage Name'].apply(lambda x: bool(re.match('^\w+[ \.]\w+$',x)))]
    lineageName = df_normalCell.loc[:,'Lineage Name']
    lineageHeader = set([re.split('[ \.]',e)[0] for e in lineageName])
    for i in lineageHeader:
        print(i)
    #there should be only aplrdv letters discribing position of the cells. take a look
    lineagePosition = set([f for e in lineageName for f in list(re.split('[ \.]',e)[1]) ])
    print(lineagePosition) #It's true that there are only 6 letters describing the cell positions
    
    #20170930
    #20171002 the meaning of P1 to P12 is different from exprected. correct the errors
    PREFIX = 'Ab:Za;AB:Za;B:Zaprppppapa;C:Zppa;D:Zpppa;E:Zpap;F:Zaplppppapp;G1:Zaprpaaaapa;G2:Zaplapaapa;H1L:Zaplaaappp;H1R:Zaarpapppp;H2L:Zaarppaaap;H2R:Zaarpppaap;In:;K:Zaplpapppaa;M:Zpaaapaapp;MS:Zpaa;P0:Z;P1:Zaplapaapp;P2:Zaprapaapp;P3:Zaplappaaa;P4:Zaprappaaa;P5:Zaplappaap;P6:Zaprappaap;P7:Zaplappapp;P8:Zaprappapp;P9:Zaplapapap;P10:Zaprapapap;P11:Zaplapappa;P12:Zaprapappa;QL:Zaplapapaaa;QR:Zaprapapaaa;TL:Zaplappppp;TR:Zaprappppp;U:Zaplppppapa;V1L:Zaarppapaa;V1R:Zaarppppaa;V2L:Zaarppapap;V2R:Zaarppppap;V3L:Zaplappapa;V3R:Zaprappapa;V4L:Zaarppappa;V4R:Zaarpppppa;V5L:Zaplapapaap;V5R:Zaprapapaap;V6L:Zaarppappp;V6R:Zaarpppppp;W:Zaprapaapa;Y:Zaprpppaaaa;Z1:Zpaapppaap;Z4:Zpaaappaap;Z:Z'
    dcLineagePrefix = dict([e.split(':') for e in PREFIX.split(';')])
    df_all_editing = df_all.copy()
    def lineagePrefixReplace(lineage,dcLineagePrefix=dcLineagePrefix):
        '''
        for lineage name like 'AB plapaaaapp' or 'H2L.aa'
        replace the prefix based on dcLineagePrefix
        if prefix is not in dcLineagePrefix, return empty string
        '''
        dc1 = dcLineagePrefix
        dc2 = dcLineagePrefix.copy()
        dc2['P4'] = 'Zpppp'
        import re
        if re.match('^\w+[ \.]\w+$',lineage):
            prefix, subfix = re.split('[ \.]',lineage)
        elif re.match('^\w+[ \.]\w+, .*$',lineage):
            prefix, subfix = re.split('[ \.]',lineage.split(',')[0])
        else:
            return ''
        if prefix in dcLineagePrefix:
            if re.match('^\w+[ ]\w+$',lineage):
                return dc2[prefix]+subfix
            else:
                return dc1[prefix]+subfix
        else:
            return ''
    df_all_editing.loc[:,'TrueLineage'] = df_all_editing.loc[:,'Lineage Name'].apply(lineagePrefixReplace,args=(dcLineagePrefix,))
    
    #20171002 all cells with a lineage starts with Z after manually annotation. 
    #read in the modifed 'C. elegans Cell List - WormAtlas.tsv' file
    df_allNew = pd.read_csv(folder + 'C. elegans Cell List - WormAtlas.tsv',sep='\t',index_col = 0)
    #a -> 0; p -> 1; a Anterior; p Posterior. Z ->1
    #l -> 0; r -> 1; l left; r right
    #d -> 0; v -> 1; d dorsal; v ventral
    df_allNew.loc[:,'Lineage (all start with 1)'] = df_allNew.loc[:,'TrueLineage'].apply(lambda x:re.sub('[ald]','0',re.sub('[Zprv]','1',x)))
    df_allNew.to_csv('list.csv')
    #check if there is duplicated 'TrueLineage'
    trueLineage = []
    for _e in list(df_allNew.loc[:,'TrueLineage']):
        if ';' not in _e:
            trueLineage.append(_e)
        else:
            trueLineage = trueLineage + _e.split(';')
    #found lineages exist more than once
    _c = Counter(trueLineage).most_common()
    trueLineageDup = [e for e in _c if e[1]>1]
    #save the duplicated lines and take a look
    trueLineageDup2 = [e[0] for e in _c if e[1]>1]
    def _f(x):
        if ';' not in x:
            return x in trueLineageDup2
        x1,x2 = x.split(';')
        if x1 in trueLineageDup2:
            return True
        if x2 in trueLineageDup2:
            return True
        return False
    df__temp = df_allNew[df_allNew.loc[:,'TrueLineage'].apply(_f)]
    #all lines were manually checked. due to naming and male/female duplicates. and some errors.
    
    #20171005
    #read in "Embryo file" and find full lineage for Dead cells.
    df_Em = pd.read_csv(folder+'CelegansEmbryonicCellLineage.alm',index_col=0,sep='\t',dtype =str)
    df_EmTrue = df_Em[df_Em.loc[:,'TrueLineage'].apply(lambda x:type(x) is str)] #store lines with known "TrueLineage"
    dc_Lineage = {}#dictionary store relationship of TrueLineage and "Lineage (all start with 1)"
    _dc =  {'Z':'Z','a':'p','l':'r','d':'v','p':'a','r':'l','v':'d'}
    _dc2 = {'Z':'1','a':'0','l':'0','d':'0','p':'1','r':'1','v':'1'}
    for _n in df_EmTrue.index:
        _v = df_EmTrue.loc[_n,'TrueLineage']
        _k = df_EmTrue.loc[_n,'Lineage(allstartwith1)']
        if len(_k) != len(_v):
            print(_k,_v,'error')
            break
        for _i in range(1,len(_k)+1):
            _kk = _k[:_i]
            _vv = _v[:_i]
            _vv2 = _vv[:-1] + _dc[_vv[-1:]]#symmetric TrueLineage only the end is different
            _kk2 = ''.join(_dc2[i] for i in _vv2)
#            print(_kk,_vv,_kk2,_vv2)
            for _kk,_vv in [(_kk,_vv),(_kk2,_vv2)]:
                if _kk not in dc_Lineage:
                    dc_Lineage[_kk] = []
                dc_Lineage[_kk].append(_vv)
    dc_Lineage={_k:set(_v) for (_k,_v) in dc_Lineage.items()}
    print([e for e in dc_Lineage.values() if len(e) >1])
    dc_Lineage={_k:list(_v)[0] for (_k,_v) in dc_Lineage.items()}
    
    #change 'TrueLineage' values for df_Em
    for _n in df_Em.index:
        if df_Em.loc[_n,'TrueLineage'] is not str:
            df_Em.loc[_n,'TrueLineage'] = dc_Lineage[df_Em.loc[_n,'Lineage(allstartwith1)']]
    df_Em.to_csv('list.csv')

def checkEpicData20170929():
    '''
    check EPIC data
    '''
    folder = r'X:\Insects\C_elegans\201709Cell_lineage\EPIC'+'\\'
    import os
    files = [folder + e for e in os.listdir(folder)]
    ls_all = []
    for _f in files:
        _ls = open(_f,encoding = 'utf-8').readlines()
        ls_all += [e.split(',')[1] for e in _ls]
    print(len(ls_all))
    for _f in files:
        _ls = open(_f,encoding = 'utf-8').readlines()
        for e in _ls:
            if e.split(',')[1] == 'Nuc1' or e.split(',')[1] == 'Nuc4':
                print(_f)
                break
    #conclusion: this file X:\Insects\C_elegans\201709Cell_lineage\EPIC\CD20070403_his72_D1A.csv
    #have two lines. one with Nuc1 and one with Nuc4. no idea what these two "cells" are.

def buildTreeForCelegansEmbryo20171006():
    '''
    buid tree with ETE with the embryo data
    '''
    #read in the embryo data
    import pandas as pd
    df_E = pd.read_csv(r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\CelegansEmbryonicCellLineage.alm",sep='\t',dtype=str,index_col=0)
    
    #initiate the tree
    from ete3 import Tree
    t = Tree()
    #use TrueLineage as name. for leaves, use name as name
    t.name = 'Z'
    for _i in df_E.index:
        _lineage = df_E.loc[_i,'Lineage(allstartwith1)']
        _TrueLineage = df_E.loc[_i,'TrueLineage']
        _name = df_E.loc[_i,'Name']
        _class = df_E.loc[_i,'Class']
        _leaf = Tree()
        _leaf.name = _name
        _leaf.add_features(Lineage = _lineage,TrueLineage = _TrueLineage, Class=_class)
        for _n in range(1,len(_TrueLineage)):
            _subTreeName = _TrueLineage[:_n]
            _upperTree = t.search_nodes(name = _subTreeName[:-1])
            if _upperTree:
                if not t.search_nodes(name = _subTreeName):
                    _subnode = Tree(name=_subTreeName)
                    if len(_upperTree) == 1:
                        _upperTree = _upperTree[0]
                        _upperTree.add_child(_subnode)
                    else:
                        print('more than one upper tree')
        _upperTree = t.search_nodes(name = _leaf.TrueLineage[:-1])[0]
        _upperTree.add_child(_leaf)
    
    #add "Lineage" feature to non-leaves
    for _leaf in t.traverse():
        if not _leaf.is_leaf():
            _leaf.add_features(Lineage = f(_leaf.name))
    #save newick tree
    t.write(features=['Lineage','TrueLineage','Class'],outfile=r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\CelegansEmbryonicCellLineageFormat8.nwk",format=8)
    ##change name of leaves to "Class"
    for _l in t.iter_leaves():
        _l.Class, _l.name = _l.name, _l.Class
    t.write(outfile=r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\CelegansEmbryonicCellLineageSimple.nwk")
    #Format8 file stored all information. including 'Lineage(allstartwith1)','TrueLineage','Class'
    
def plotTreeForCelegansEmbryo20171008():
    #plot embryo tree. align cell names. color by class
    #    Class	Count	description
    #    Neu	222	neuron
    #    Dea	113	dead
    #    Str	46	neural structural
    #    Epi	93	epithelial
    #    Mus	123	muscle
    #    Bla	39	blast cell
    #    Gla	13	gland
    #    Int	20	intestine
    #    Ger	2	germline
    from ete3 import Tree, TreeStyle, NodeStyle, AttrFace
    from matplotlib import colors
    import numpy as np
    colordc = {'Neu': 'Red', 'Dea': 'grey', 'Str' : 'Magenta' , 'Epi' : 'Green' , 'Mus' : 'MediumBlue' , 'Bla' : 'LightGreen' , 'Gla' : 'Brown' , 'Int' : 'Cyan' , 'Ger' : 'Orange' }
    #read in the tree
    t = Tree(r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\CelegansEmbryonicCellLineageFormat8.nwk",format=8)
    ## add AttrFace
    #iterate all nodes. set node to 0. color leaves
    for _l in t.traverse():#set all leaves
        ns = NodeStyle()
        ns['size'] = 0
        if _l.is_leaf():
            _n = _l.Class
            ns['hz_line_width'] = 1
            ns['vt_line_width'] = 1
            ns['hz_line_color'] = colordc[_n]
            ns['vt_line_color'] = colordc[_n]
            _l.img_style = ns
    for _l in t.traverse():
        ns = NodeStyle()
        ns['size'] = 0
        if not _l.is_leaf(): #average color from leaves for internal nodes
            ns['hz_line_width'] = 1
            ns['vt_line_width'] = 1
            ns['hz_line_type'] = 1
            _r,_g,_b = [],[],[]
            for _ll in _l.iter_leaves():
                _color = colors.ColorConverter.to_rgb(_ll.img_style['hz_line_color'])
                _r.append(_color[0])
                _g.append(_color[1])
                _b.append(_color[2])
            _r = np.mean(_r)
            _g = np.mean(_g)
            _b = np.mean(_b)
#            print(_r,_g,_b,colors.rgb2hex((_r,_g,_b)))
#            break
            ns['hz_line_color'] = colors.rgb2hex((_r,_g,_b))
            ns['vt_line_color'] = colors.rgb2hex((_r,_g,_b))
            _l.img_style = ns
    
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.min_leaf_separation = 2
    ts.draw_guiding_lines = True
    ts.show_scale = False
    ts.scale =  30
    t.show(tree_style = ts)
    
    # color tree based on average of descendents. 
    from ete3 import Tree, TreeStyle, NodeStyle, AttrFace, faces
    from matplotlib import colors
    import numpy as np
    colordc = {'Neu': 'Red', 'Dea': 'grey', 'Str' : 'Magenta' , 'Epi' : 'Green' , 'Mus' : 'MediumBlue' , 'Bla' : 'LightGreen' , 'Gla' : 'Brown' , 'Int' : 'Cyan' , 'Ger' : 'Orange' }
    #read in the tree
    t = Tree(r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\CelegansEmbryonicCellLineageFormat8.nwk",format=8)
    #swap all nodes
    for _l in t.traverse("postorder"):
        if not _l.is_leaf():
            _l.swap_children()
    for _l in t.traverse("postorder"):
        ns = NodeStyle()
        ns['size'] = 0
        if _l.is_leaf():
            _n = _l.Class
            ns['hz_line_width'] = 1
            ns['vt_line_width'] = 1
            ns['hz_line_color'] = colordc[_n]
            ns['vt_line_color'] = colordc[_n]
            _l.img_style = ns
        else: #average color children for internal nodes
            ns['hz_line_width'] = 1
            ns['vt_line_width'] = 1
            ns['hz_line_type'] = 1
            _r,_g,_b = [],[],[]
            for _ll in _l.children:
                _color = colors.ColorConverter.to_rgb(_ll.img_style['hz_line_color'])
                _r.append(_color[0])
                _g.append(_color[1])
                _b.append(_color[2])
            _r = np.mean(_r)
            _g = np.mean(_g)
            _b = np.mean(_b)
#            print(_r,_g,_b,colors.rgb2hex((_r,_g,_b)))
#            break
            ns['hz_line_color'] = colors.rgb2hex((_r,_g,_b))
            ns['vt_line_color'] = colors.rgb2hex((_r,_g,_b))
            _l.img_style = ns
    ts = TreeStyle()
    t.dist = 5
    ## add legend
    ts.legend.clear()
    ts.title.clear()
    _num = 1
    for _n in colordc:
        ts.legend.add_face(faces.TextFace(_n+' ',fgcolor = colordc[_n],ftype = 'Arial'),_num)
        _num += 1
    ts.legend_position = 3
    ts.show_leaf_name = False
    ts.min_leaf_separation = 1
    ts.draw_guiding_lines = True
    ts.show_scale = False
    ts.scale =  30
    ts.mode = 'c'
    ts.draw_guiding_lines = False
    ts.complete_branch_lines_when_necessary = True
    ts.arc_start = 0
    ts.arc_span =  180
    ts.root_opening_factor = 0.5
    t.show(tree_style = ts)
    
def treeComparisonResultChecking20171010():
    '''
    tree comparison
    '''
    # check score matrix
    f_score = r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\Lineages\cost.tsv"
    import pandas as pd
    df_score = pd.DataFrame()
    for _l in open(f_score):
        _r,_c,_s = _l.split()
        df_score.loc[_r,_c] = int(_s)
    ##conclusion: cell types: 'Bla Dea Epi Ger Gla Int Mus Neu Str'
    ##score 2 for identical cell types, and 0 for non-identical pairs
    
    #read embryo Tree
    from ete3 import Tree, TreeStyle, NodeStyle, AttrFace, faces, TreeFace
    from matplotlib import colors
    import numpy as np
    colordc = {'Neu': 'Red', 'Dea': 'grey', 'Str' : 'Magenta' , 'Epi' : 'Green' , 'Mus' : 'MediumBlue' , 'Bla' : 'LightGreen' , 'Gla' : 'Brown' , 'Int' : 'Cyan' , 'Ger' : 'Orange' }
    #read in the tree
    t = Tree(r"C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\CelegansEmbryonicCellLineageFormat8.nwk",format=8)
    #swap all nodes
    ##set up styles for nodes and leaves
    for _l in t.traverse("postorder"):
        if not _l.is_leaf():
            _l.swap_children()
    for _l in t.traverse("postorder"):
        ns = NodeStyle()
        ns['size'] = 0
        if _l.is_leaf():
            _n = _l.Class
            ns['hz_line_width'] = 1
            ns['vt_line_width'] = 1
            ns['hz_line_type'] = 0
            ns['size'] = 0
            ns['hz_line_color'] = colordc[_n]
            ns['vt_line_color'] = colordc[_n]
            _l.img_style = ns
        else: #average color children for internal nodes
            ns['hz_line_width'] = 1
            ns['vt_line_width'] = 1
            ns['hz_line_type'] = 1
            _r,_g,_b = [],[],[]
            for _ll in _l.children:
                _color = colors.ColorConverter.to_rgb(_ll.img_style['hz_line_color'])
                _r.append(_color[0])
                _g.append(_color[1])
                _b.append(_color[2])
            _r = np.mean(_r)
            _g = np.mean(_g)
            _b = np.mean(_b)
#            print(_r,_g,_b,colors.rgb2hex((_r,_g,_b)))
#            break
            ns['hz_line_color'] = colors.rgb2hex((_r,_g,_b))
            ns['vt_line_color'] = colors.rgb2hex((_r,_g,_b))
            _l.img_style = ns
    #create a dictionary with 0101 code as keys
    dcNodes = {}
    for _l in t.iter_descendants():
        dcNodes[_l.Lineage[1:]] = _l
    tt = t.copy()
    #add face of lineage
    for _l in tt.iter_descendants():
        _l.add_face(faces.TextFace(_l.Lineage[1:]), 0 ,position = 'branch-top')
    
    dcNodes2 ={_l.Lineage[1:]:_l for _l in tt.iter_descendants()}
    
    alignAnnotation = '50'
    Score = '43'
    RootS = '00011'
    RootT = '00000'
    PruneS = '000111110 '.split()
    PruneT = '0000010001 '.split()
    MatchS = [RootS] + '000110 000111 0001100 0001101 0001110 0001111 00011000 00011001 00011010 00011011 00011100 00011101 00011110 00011111 000110000 000110001 000110010 000110011 000110100 000110101 000110110 000110111 000111000 000111001 000111010 000111011 000111100 000111101 0001100010 0001100011 0001100110 0001100111 0001101010 0001101011 0001101100 0001101101 0001101110 0001101111 0001110000 0001110001 0001110010 0001110011 0001110100 0001110101 0001110110 0001110111 0001111000 0001111001 0001111010 0001111011 0001111110 0001111111 00011011100 00011011101 00011111100 00011111101 00011111110 00011111111'.split()
    MatchT = [RootT] + '000001 000000 0000010 0000011 0000001 0000000 00000100 00000101 00000110 00000111 00000010 00000011 00000000 00000001 000001000 000001001 000001010 000001011 000001100 000001101 000001110 000001111 000000100 000000101 000000111 000000110 000000000 000000001 0000010010 0000010011 0000010110 0000010111 0000011010 0000011011 0000011100 0000011101 0000011110 0000011111 0000001001 0000001000 0000001011 0000001010 0000001110 0000001111 0000001101 0000001100 0000000000 0000000001 0000000011 0000000010 000000010 000000011 00000111100 00000111101 0000000101 0000000100 0000000110 0000000111'.split()
    
    def plotAlignedTree(fullTreeWithStyle, Score, RootS, RootT, PruneS, PruneT,MatchS, MatchT,alignAnnotation = '', outfile = None, showOrder = False):
        '''
        fullTreeWithStyle is the full tree with the style for each node set up.
        Score is the Score
        RootS, RootT is the root node of the two tree
        PruneS, PruneT are lists of nodes to be pruned.
        alignAnnotation stores the information which the user want to include in the 
        if showOrder = True, add order numbers to the branches. add Lineage numbers
        '''
        tt = fullTreeWithStyle.copy()#copy the input tree to avoid changes to the original tree
        if showOrder:
            for _l in tt.iter_descendants():
                _l.add_face(faces.TextFace(_l.Lineage[1:]), 0 ,position = 'branch-top')
        dcNodes2 = {_l.Lineage[1:]:_l for _l in tt.iter_descendants()}
        Stree = dcNodes2[RootS].copy()
        Ttree = dcNodes2[RootT].copy()
        #remove branches with no match
        ##define a function to decide whether a node to be removed
        def toRemove(_l,matches,prunes):
            if  _l.Lineage[1:] not in matches:
                if _l.Lineage[1:] not in prunes:
                    _lsis = _l.get_sisters()
                    if _lsis:
                        _lsis = _lsis[0]
#                        if _lsis.is_leaf() and _lsis.Lineage[1:] in prunes:
                        if _lsis.Lineage[1:] in prunes:
                            return False
                    return True
                else:
                    _lsis = _l.get_sisters()
                    if _lsis:
                        _lsis = _lsis[0]
                        if _lsis.is_leaf():
                            return True
                        
            return False
        def removeBranches(atree,matches,prunes):#remove unmatched end nodes recursively
            setmatches = set(matches)
            while True:
                _list = [_l for _l in atree.iter_leaves() if toRemove(_l,setmatches,prunes)]
                if not _list:
                    break
                for _l in _list:
                    _l.detach()
            return atree
        Stree = removeBranches(Stree, MatchS, PruneS)
        Ttree = removeBranches(Ttree, MatchT, PruneT)
        
#        def toRemoveInternal(_l,matches,prunes):
#            if _l.is_leaf():
#                return False
#            _l_left, _l_right = _l.children
#            if _l_left.Lineage[1:] not in matches and _l_right.Lineage[1:] not in matches:
#                if _l_left.Lineage[1:] not in prunes and _l_right.Lineage[1:] not in prunes:
#                    for _l_children in _l.iter_descendants:
#                        if _l_children.Lineage[1:] in prunes:
#                            return True
#            return False
#        def removeInternal(atree,matches,prunes):
#            setmatches = set(matches)
#            while True:
#                _list = [_l for _l in atree.traverse() if toRemoveInternal(_l,setmatches,prunes)]
#                if not _list:
#                    break
#                for _l in _list:
#                    _l.dist += _l.children[0].dist
#                    if len(_l.children) <2:
#                        print(_l.Lineage)
#                    _l_left, _l_right = _l.children
#                    _l_left.delete()
#                    _l_right.delete()
#            return atree
#        Stree = removeInternal(Stree,MatchS, PruneS)
#        Ttree = removeInternal(Ttree, MatchT, PruneT)
#        #in some cases, the two or more layers of nodes were trimmed. 
#        #Only the bottom trimmed branch were included. Simplify the tree further by combine nodes with only one children
        
        #find trimmed leaves
        S_Trimedleaves = set()
        T_Trimedleaves = set()
        for _l in Stree.iter_leaves():
            if _l.Lineage[1:] not in MatchS and _l.up.Lineage[1:] in MatchS and not (_l.Lineage[1:] in PruneS or _l.get_sisters()[0].Lineage[1:] in PruneS):
                S_Trimedleaves.add(_l.up)
        for _l in Ttree.iter_leaves():
            if _l.Lineage[1:] not in MatchT and _l.up.Lineage[1:] in MatchT and not (_l.Lineage[1:] in PruneT or _l.get_sisters()[0].Lineage[1:] in PruneT):
                T_Trimedleaves.add(_l.up)
        for _l in S_Trimedleaves:
            PruneS.append(_l.children[0].Lineage[1:])
        for _l in T_Trimedleaves:
            PruneT.append(_l.children[0].Lineage[1:])
        
                
        
        _n = len(MatchS) #add order number to nodes
        for _l in Stree.traverse():
            if _l.Lineage[1:] in MatchS:
                if showOrder:
                    _l.add_face(faces.TextFace(str(MatchS.index(_l.Lineage[1:]))), 0, position = 'branch-bottom')#add face of order
                _l.add_features(orders = MatchS.index(_l.Lineage[1:]))
            else:
                _l.add_features(orders = _n)
                _n += 1
        for _l in Ttree.traverse():
            if _l.Lineage[1:] in MatchT:
                if showOrder:
                    _l.add_face(faces.TextFace(str(MatchT.index(_l.Lineage[1:]))), 0, position = 'branch-bottom')#add face of order
                _l.add_features(orders = MatchT.index(_l.Lineage[1:]))
            else:
                _l.add_features(orders = _n)
                _n += 1
        
        def getNodeStyle(size = 0, vt_line_color='LightGrey', hz_line_color='LightGrey', hz_line_width = 1, vt_line_width = 1, hz_line_type = 2, vt_line_type = 2):
            ns = NodeStyle()
            ns['size'] = size
            ns['vt_line_color'] = vt_line_color
            ns['hz_line_color'] = hz_line_color
            ns['hz_line_width'] = hz_line_width
            ns['vt_line_width'] = vt_line_width
            ns['hz_line_type'] = hz_line_type
            ns['vt_line_type'] = vt_line_type
            return ns
        greyStyle = getNodeStyle()
            
        ##add trimmed nodes back
        dcNodesT = {_l.Lineage[1:]:_l for _l in Ttree.iter_descendants()}
        dcNodesS = {_l.Lineage[1:]:_l for _l in Stree.iter_descendants()}
        PruneS = [e for e in PruneS if e in dcNodesS]
        PruneT = [e for e in PruneT if e in dcNodesT]
        for _p in PruneS:
            _node = dcNodesS[_p]
            _node_sis = _node.get_sisters()[0]
            _node_up = _node.up
            _node_upT = dcNodesT[MatchT[MatchS.index(_node_up.Lineage[1:])]]
            _node2 = _node.copy()
            _node_sis2 = _node_sis.copy()
            for _l in _node2.traverse():#set greystyle for added nodes
                _l.img_style = greyStyle
                _l.del_feature('Class')#remove feature class
            for _l in _node_sis2.traverse():
                _l.img_style = greyStyle
                _l.del_feature('Class')#remove feature class
            for _l in _node.traverse():#set linetype to 2 for non-matched nodes
                _l.img_style['vt_line_type'] = _l.img_style['hz_line_type'] = 2
            _node_sis.img_style['vt_line_type'] = _node_sis.img_style['hz_line_type'] = 2
#            if _node.is_leaf():#add two leaves
#                _node_upT.add_child(_node2)
#                _node_upT.add_child(_node_sis2)
#            else:#add a new layer
            _node_upT_child = [e.copy() for e in _node_upT.children]
            _node_sis2.children = []
            for _e in _node_upT_child:
                _node_sis2.add_child(_e)
            _node_upT.children = []
            _node_upT.add_child(_node_sis2)
            _node_upT.add_child(_node2)
            dcNodesT[_node_sis2.Lineage[1:]] = _node_sis2#update the dictionary
            dcNodesT[_node2.Lineage[1:]] = _node2
            MatchS += [_node2.Lineage[1:], _node_sis2.Lineage[1:]] #update the MatchS/MatchT
            MatchT += [_node2.Lineage[1:], _node_sis2.Lineage[1:]]
        
        for _p in PruneT:
            _node = dcNodesT[_p]
            _node_sis = _node.get_sisters()[0]
            _node_up = _node.up
            _node_upS = dcNodesS[MatchS[MatchT.index(_node_up.Lineage[1:])]]
            _node2 = _node.copy()
            _node_sis2 = _node_sis.copy()
            for _l in _node2.traverse():
                _l.img_style = greyStyle
                _l.del_feature('Class')#remove feature class
            for _l in _node_sis2.traverse():
                _l.img_style = greyStyle
                _l.del_feature('Class')#remove feature class
            for _l in _node.traverse():#set linetype to 2 for non-matched nodes
                _l.img_style['vt_line_type'] = _l.img_style['hz_line_type'] = 2
            _node_sis.img_style['vt_line_type'] = _node_sis.img_style['hz_line_type'] = 2
#            if _node.is_leaf():#add two leaves
#                _node_upS.add_child(_node2)
#                _node_upS.add_child(_node_sis2)
#            else:#add a new layer
            _node_upS_child = [e.copy() for e in _node_upS.children]
            _node_sis2.children = []
            for _e in _node_upS_child:
                _node_sis2.add_child(_e)
            _node_upS.children = []
            _node_upS.add_child(_node_sis2)
            _node_upS.add_child(_node2)
            dcNodesS[_node_sis2.Lineage[1:]] = _node_sis2#update the dictionary
            dcNodesS[_node2.Lineage[1:]] = _node2
            MatchS += [_node2.Lineage[1:], _node_sis2.Lineage[1:]] #update the MatchS/MatchT
            MatchT += [_node2.Lineage[1:], _node_sis2.Lineage[1:]]
    
        
        ##sort node recursivelyfor _l in Stree.traverse():
        Ttree.sort_descendants(attr='orders')
        Stree.sort_descendants(attr='orders')
        
        #change line width for identical leaves 
        for _lT, _lS in zip(Ttree.iter_leaves(),Stree.iter_leaves()):
            if 'Class' in _lT.features and 'Class' in _lS.features:
                if _lT.Class == _lS.Class:
                    if _lT.Lineage[1:] in MatchT and _lS.Lineage[1:] in MatchS:
                        _lT.img_style['hz_line_width'] = 2
                        _lS.img_style['hz_line_width'] = 2
        
        #create TreeFace
        tsT = TreeStyle()
        tsT.show_leaf_name = False
        tsT.min_leaf_separation = 5
    #    tsT.optimal_scale_level = 'full'
        tsT.show_scale = False
        tsT.force_topology = True
        tsT.orientation = 0
        tsS = TreeStyle()
        tsS.show_leaf_name = False
        tsS.min_leaf_separation = 5
    #    tsS.optimal_scale_level = 'full'
        tsS.show_scale = False
        tsS.force_topology = True
        tsS.orientation = 1
        Tbranch = TreeFace(Ttree,tsT)
        Sbranch = TreeFace(Stree,tsS)
        
        midTree = Tree('(S,T);')#to combine two trees together
        ns = NodeStyle()
        ns['size'] = 0
        ns['hz_line_color'] = '#FFFFFF'
        for _l in midTree.traverse():
            _l.img_style = ns
        tsMid = TreeStyle()
        tsMid.mode = 'c'
        tsMid.arc_start = 0
        tsMid.arc_span = 360
        tsMid.show_leaf_name = False
        tsMid.show_scale = False
        tsMid.scale = 20
        tsMid.title.add_face(faces.TextFace(text = alignAnnotation +' Score: {0}  RootS: {1}  RootT: {2}'.format(Score, RootS, RootT)), column = 0)
        S_face,T_face = midTree.children
        S_face.add_face(Sbranch,0)
        S_face.add_face(faces.TextFace('S'),1)
        T_face.add_face(Tbranch,0)
        T_face.add_face(faces.TextFace('T'),1)
        if outfile is None:
            midTree.show(tree_style=tsMid)
        else:
            midTree.render(outfile, tree_style=tsMid)
    
    
    #test
#    plotAlignedTree(t,Score,RootS, RootT, PruneS, PruneT,alignAnnotation, outfile = None)
    
    #plot all alignments from C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\Lineages\20171011CelegansEmbryoLineage.txtl
    filename = r'C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\Lineages\20171011CelegansEmbryoLineage.txtl'
    folder = r"X:\Insects\C_elegans\201710TreeAlignment"+'\\'
    txt = open(filename).read()
    import re
    alignments = re.findall('\d*\{.*?\}',txt,re.DOTALL)
    for alignment in alignments[1:]:
        _ls = alignment.split('\n')
        alignAnnotation = _ls[0].split('{')[0]
        Score = _ls[1].split(':')[1]
        RootS = _ls[2].split(':')[1]
        RootT = _ls[3].split(':')[1]
        PruneS = _ls[4].split(':')[1].split()
        PruneT = _ls[5].split(':')[1].split()
        MatchS = _ls[6].split(':')[1].split()
        MatchT = _ls[7].split(':')[1].split()
#        break
#        if alignAnnotation == '2':
#            break
        plotAlignedTree(t,Score,RootS, RootT, PruneS, PruneT,MatchS, MatchT, alignAnnotation, outfile = folder+alignAnnotation+'.pdf')
#        plotAlignedTree(t,Score,RootS, RootT, PruneS, PruneT,MatchS, MatchT, alignAnnotation, outfile = folder+alignAnnotation+'.pdf')
        
        
    #merge pdf files
    from PyPDF2 import PdfFileReader, PdfFileWriter
    # Creating a routine that appends files to the output file
    def append_pdf(infile,outfile):
        [output.addPage(infile.getPage(page_num)) for page_num in range(infile.numPages)]
    
    # Creating an object where pdf pages are appended to
    output = PdfFileWriter()
    
    # Appending two pdf-pages from two different files
    files = [folder + str(e)+'.pdf' for e in range(2,101)]
    for _f in files:
        append_pdf(PdfFileReader(open(_f,"rb")),output)
    
    # Writing all the collected pages to a file
    output.write(open(folder+"../20171017CombinedPages.pdf","wb"))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    