# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 20:37:55 2017

@author: k
"""

class network(object):
    '''
    a network to show molecules linkage
    knode were stored in dict, with key as the knode label, and content stores knode information
    links were stored in dict, with knode label as key, and a list of linked knode labels
    '''
    def __init__(self):
        from collections import defaultdict
        self.knodes = {}
        self.links = defaultdict(set)
    def add_knodes(self, label, knode_content):
        '''
        add knode to the network
        if label is already used, print warning and do nothing
        '''
        if label in self.knodes:
            print('label already used, cannot add knode')
        else:
            self.knodes[label] = knode_content
    def add_links(self, label1, label2, bidirectional = True):
        if label1 not in self.knodes:
            print('label1 not exist! cannot add link!')
        elif label2 not in self.knodes:
            print('label2 not exist! cannot add link!')
        elif bidirectional:
            self.links[label1].add(label2)
            self.links[label2].add(label1)
        else:
            self.links[label1].add(label2)
    def remove_links(self,label1, label2, bidirectional = True):
        if label1 not in self.knodes:
            print('label1 not exist! cannot add link!')
        elif label2 not in self.knodes:
            print('label2 not exist! cannot add link!')
        elif bidirectional:
            self.links[label1].remove(label2)
            self.links[label2].remove(label1)
        else:
            self.links[label1].remove(label2)
    def get_links(self, bidirectional = True):
        '''
        return self.links in a list of tuple
        '''
        links_in_tuple = []
        for key, values in self.links.items():
            for value in values:
                links_in_tuple.append((key,value))
        if bidirectional:
            links_bi = set()
            for linktuple in links_in_tuple:
                if linktuple not in links_bi:
                    reversetuple = (linktuple[1], linktuple[0])
                    if reversetuple not in links_bi:
                        links_bi.add(linktuple)
            return list(links_bi)
        return links_in_tuple
    def get_independent_map(self):
        '''
        return a list net, each is a list of knodes labels that can be linked together.
        '''
        def getmap(knodes, links):
            knodes = list(knodes)
            if len(knodes) == 0:
                print('empty network!')
                return None
            if len(knodes) == 1:
                return [knodes]
            else:
                independentmaps = getmap(knodes[1:], links)
                knode = knodes[0]
                knodeIn = False
                for independentmap in independentmaps:
                    linked = set()
                    for k in independentmap:
                        linked.update(links[k])
                    if knode in linked:
                        independentmap.append(knode)
                        knodeIn = True
                        break
                if not knodeIn:
                    independentmaps.append([knode])
                return independentmaps
        return getmap(list(self.knodes.keys()),self.links)

def createPGN():
    '''
    create molecule of PGN monomer in network format
    '''
    PGN = network()
    PGN.add_knodes(1,'H')
    PGN.add_knodes(2,'H')
    PGN.add_knodes(3,'N')
    PGN.add_knodes(4,'N')
    PGN.add_knodes(5,'L')
    PGN.add_knodes(6,'A')
    PGN.add_knodes(7,'G')
    PGN.add_knodes(8,'D')
    PGN.add_knodes(9,'A')
    PGN.add_links(1,3)
    PGN.add_links(2,4)
    PGN.add_links(3,4)
    PGN.add_links(4,5)
    PGN.add_links(5,6)
    PGN.add_links(6,7)
    PGN.add_links(7,8)
    PGN.add_links(8,9)
    return PGN

def createPGN4():
    '''
    create molecule of PGN with four PGN monomers
    '''
    PGN4 = network()
    PGN4.add_knodes(1,'H')
    PGN4.add_knodes(2,'H')
    PGN4.add_knodes(3,'N')
    PGN4.add_knodes(4,'N')
    PGN4.add_knodes(5,'L')
    PGN4.add_knodes(6,'A')
    PGN4.add_knodes(7,'G')
    PGN4.add_knodes(8,'D')
    PGN4.add_knodes(9,'A')
    PGN4.add_links(1,3)
    PGN4.add_links(2,4)
    PGN4.add_links(3,4)
    PGN4.add_links(4,5)
    PGN4.add_links(5,6)
    PGN4.add_links(6,7)
    PGN4.add_links(7,8)
    PGN4.add_links(8,9)
    
    PGN4.add_knodes(11,'H')
    PGN4.add_knodes(12,'H')
    PGN4.add_knodes(13,'N')
    PGN4.add_knodes(14,'N')
    PGN4.add_knodes(15,'L')
    PGN4.add_knodes(16,'A')
    PGN4.add_knodes(17,'G')
    PGN4.add_knodes(18,'D')
    PGN4.add_knodes(19,'A')
    PGN4.add_links(11,13)
    PGN4.add_links(12,14)
    PGN4.add_links(13,14)
    PGN4.add_links(14,15)
    PGN4.add_links(15,16)
    PGN4.add_links(16,17)
    PGN4.add_links(17,18)
    PGN4.add_links(18,19)
    PGN4.add_knodes(21,'H')
    PGN4.add_knodes(22,'H')
    PGN4.add_knodes(23,'N')
    PGN4.add_knodes(24,'N')
    PGN4.add_knodes(25,'L')
    PGN4.add_knodes(26,'A')
    PGN4.add_knodes(27,'G')
    PGN4.add_knodes(28,'D')
    PGN4.add_knodes(29,'A')
    PGN4.add_links(21,23)
    PGN4.add_links(22,24)
    PGN4.add_links(23,24)
    PGN4.add_links(24,25)
    PGN4.add_links(25,26)
    PGN4.add_links(26,27)
    PGN4.add_links(27,28)
    PGN4.add_links(28,29)
    
    PGN4.add_knodes(31,'H')
    PGN4.add_knodes(32,'H')
    PGN4.add_knodes(33,'N')
    PGN4.add_knodes(34,'N')
    PGN4.add_knodes(35,'L')
    PGN4.add_knodes(36,'A')
    PGN4.add_knodes(37,'G')
    PGN4.add_knodes(38,'D')
    PGN4.add_knodes(39,'A')
    PGN4.add_links(31,33)
    PGN4.add_links(32,34)
    PGN4.add_links(33,34)
    PGN4.add_links(34,35)
    PGN4.add_links(35,36)
    PGN4.add_links(36,37)
    PGN4.add_links(37,38)
    PGN4.add_links(38,39)
    
    PGN4.add_links(9,18)
    PGN4.add_links(19,28)
    PGN4.add_links(29,38)
    
    PGN4.add_links(4,23)
    PGN4.add_links(13,34)
    return PGN4

from pyteomics import mass
dc = {}
dc['N'] = mass.Composition(formula = 'C8H15NO6')
dc['L'] = mass.Composition(formula = 'C3H6O3')
dc['A'] = mass.Composition(formula = 'C3H7NO2')
dc['G'] = mass.Composition(formula = 'C5H9NO4')
dc['D'] = mass.Composition(formula = 'C7H14N2O4')
dc['H'] = mass.Composition()
dc_DAPPG=dc
dc = {}
dc['N'] = mass.Composition(formula = 'C8H15NO6')
dc['L'] = mass.Composition(formula = 'C3H6O3')
dc['A'] = mass.Composition(formula = 'C3H7NO2')
dc['G'] = mass.Composition(formula = 'C5H9NO4')
dc['D'] = mass.Composition(formula = 'C6H14N2O2')
dc['H'] = mass.Composition()
dc_LysPG = dc

def getCompositionAndMW(s,waterCount,dc=dc):
    '''
    given a string of single letter of dc keys, and water to remove, return composition
    '''
    m = mass.Composition()
    for k in s:
        m = m+dc[k]
    m = m - waterCount * mass.Composition(formula = 'H2O')
    f = ''
    for a,n in m.items():
        f = f + a +str(n)
    mw = mass.calculate_mass(formula = f)
    return f,mw

def getPGNproduct(PGN, length, lsHelper = None):
    '''
    get all possible product from molecule PGN with defined length
    lsHelper is result of length-1. 
    '''
    if length == 1:
        return set([(product,) for product in PGN.knodes.keys()])
    else:
        if lsHelper is None:
            products = getPGNproduct(PGN, length -1)
        else:
            products = lsHelper
        newproducts = set()
        for product in products:
            linked = set()#store possible linked knode_labels
            for knode in product:
                linked.update(PGN.links[knode])
            for knode in product:#remove labels already used
                if knode in linked:
                    linked.remove(knode)
            for newknode in linked:#add one more knode
                newproduct = list(product)
                newproduct.append(newknode)
                newproduct = tuple(sorted(newproduct))
                newproducts.add(newproduct)
        return newproducts

def getPGNallProducts(PGN):
    '''
    return a list of tuple, with string of composition and water molecules
    '''
    productlist = []
    productlist.append(getPGNproduct(PGN,1))
    maxlength = len(PGN.knodes)
    for i in range(2,maxlength+1):
        productlist.append(getPGNproduct(PGN,i,productlist[-1]))
    products = []
    for _e in productlist:
        products += _e
    productsComp = set()
    two_circleProduct = 'AAAAAAADDDDGGGGLLLLNNNNNN'
    one_circleProduct = 'AAAADDDGGLLNNN'
    from collections import Counter
    counter_1circle = Counter(one_circleProduct)
    counter_2circle = Counter(two_circleProduct)
    for product in products:
        s = ''
        for knode_key in product:
            s += PGN.knodes[knode_key] #get molecule in str
        s = ''.join(sorted(s))#sort string for reduce redundancy
        watercount = len(s)-1
        productsComp.add((s,watercount))
        if len(s)>= 25:
            counter_s = Counter(s)
            if False not in [counter_s[k] >= counter_2circle[k] for k in counter_2circle]:
                productsComp.add((s,watercount-1))
                productsComp.add((s,watercount-2))
        elif len(s)>= 14:
            counter_s = Counter(s)
            if False not in [counter_s[k] >= counter_1circle[k] for k in counter_1circle]:
                productsComp.add((s,watercount-1))
    return productsComp





def getCompositionFromMWsimple(mw, ppm=10, water_error=True):
    '''
    dc is a dictionary with sub molecule names and their composition
    return possible compositions of with molecules from dc, within allowed error range
    water removed is the number of subunits -1. if water_error, can add or remove water equal or less than the number of GlcN
    '''
    from pyteomics import mass
#    mw = 161.0688078
#    ppm = 10
#    water_error = True
    dc_PGcomplex = {}
    dc_PGcomplex['GlcN'] = mass.Composition(formula = 'C6H13NO5') #sugar
    dc_PGcomplex['Ac'] = mass.Composition(formula = 'C2H4O2') #acetic acid
    dc_PGcomplex['Lac'] = mass.Composition(formula = 'C3H6O3') #lactate
    dc_PGcomplex['Ala'] = mass.Composition(formula = 'C3H7NO2')
    dc_PGcomplex['Glu'] = mass.Composition(formula = 'C5H9NO4')
    dc_PGcomplex['Gln'] = mass.Composition(formula = 'C5H10N2O3')
    dc_PGcomplex['Lys'] = mass.Composition(formula = 'C6H14N2O2')
    dc_PGcomplex['DAP'] = mass.Composition(formula = 'C7H14N2O4')
    dc_PGcomplex['Gly'] = mass.Composition(formula = 'C2H5NO2')
    dc = dc_PGcomplex
    if water_error:
        dc['GlcN-H2O'] = mass.Composition(formula = 'C6H11NO4')
        dc['GlcN+H2O'] = mass.Composition(formula = 'C6H15NO6')
        dc['GlcNH'] = mass.Composition(formula = 'C6H15NO5')
        dc['GlcNHH'] = mass.Composition(formula = 'C6H17NO5')
        dc['Asp'] = mass.Composition(formula = 'C4H7NO4')
        dc['Leu'] = mass.Composition(formula = 'C6H13NO2')
        dc['Met'] = mass.Composition(formula = 'C5H11NO2S')
        dc['Phe'] = mass.Composition(formula = 'C9H11NO2')
        dc['Thr'] = mass.Composition(formula = 'C4H9NO3')
#        dc['Trp'] = mass.Composition(formula = 'C11H12N2O2')
        dc['Val'] = mass.Composition(formula = 'C5H11NO2')
        dc['Arg'] = mass.Composition(formula = 'C6H14N4O2')
        dc['His'] = mass.Composition(formula = 'C6H9N3O2')
        dc['Asn'] = mass.Composition(formula = 'C4H8N2O3')
#        dc['Pro'] = mass.Composition(formula = 'C5H9NO2')
        dc['Ser'] = mass.Composition(formula = 'C3H7NO3')
        dc['Tyr'] = mass.Composition(formula = 'C9H11NO3')
        
    
    dcMW = {k:mass.calculate_mass(v) for k,v in dc.items()}
    waterMW = mass.calculate_mass(formula = 'H2O')
    mwl = mw*(1-ppm/1000000)
    mwh = mw*(1+ppm/1000000)
    
    def mw_in(m,mw=mw,ppm=ppm):
        return m>=mwl and m<=mwh
    
    def calculateMW(comp):
        M = 0
        comp = list(comp)
        water_RM = max(sum(comp)-1,0)
        _mws = list(dcMW.values())
        M = sum([a*b for a,b in zip(_mws,comp)])
        return M-waterMW*water_RM
    
    keys = list(dcMW.keys())
    
    def getComp(mw=mw,comps=None,good=[]):
        if comps is None:
            comps =[[0 for e in keys]]
        for comp in comps:
            mw_comp = calculateMW(comp)
#            print(mw_comp)
            if mw_in(mw_comp):
                _good = comp.copy()
                good.append(_good)
#                print(comp)
                
        comps = [comp for comp in comps if calculateMW(comp)<mwl]
#        print(len(comps))
        if len(comps) == 0:
            return good
        else:
            newcomps = []
            for comp in comps:
                for i in range(len(comp)):
                    comp2 = comp.copy()
                    comp2[i] += 1
                    newcomps.append(comp2)
            newcomps = set([tuple(e) for e in newcomps])
            newcomps = [list(e) for e in newcomps]
#            print(len(newcomps))
            return getComp(mw,newcomps,good)
    
    resultgood = getComp(mw)
    products = []
    for comp in resultgood:
        description = ''
        product = mass.Composition()
        for i in range(len(comp)):
            product += dc[keys[i]]*comp[i]
            if comp[i]>0:
                description = description + keys[i]+str(comp[i])
        product = product - mass.Composition(formula = 'H2O') * (sum(comp)-1)
        f = ''
        for a,n in product.items():
            f = f + a +str(n)
        productms = mass.calculate_mass(product)
        productppm = abs(productms-mw)/mw*1000000
        
        #for products, only keep the most possible ones
        if water_error:
            if comp[-1]>0 and comp[-2]>0:
                continue
        products.append([description, f,productms,productppm])
#    print('construct, formula, MW, ppm')
    return products

from pyteomics import mass
test1 = getCompositionFromMWsimple(319.2355679)
test2 = getCompositionFromMWsimple(419.2514375)
test3 = getCompositionFromMWsimple(602.4138936)
test4 = getCompositionFromMWsimple(881.5606639)
test5 = getCompositionFromMWsimple(1164.738799)
test6 = getCompositionFromMWsimple(mass.calculate_mass('DDEE'))
