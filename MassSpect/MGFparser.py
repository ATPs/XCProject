# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:40:52 2015

@author: k
"""
import numpy as np

class MGF(object):
    """
    store mgf file elements.
    mgf file looks like:
    "
    BEGIN IONS
    TITLE=1380_Control_m100min090215.49776.49776.2 File:"1380_Control_m100min090215.raw", NativeID:"controllerType=0 controllerNumber=1 scan=49776"
    RTINSECONDS=6054.542088
    PEPMASS=825.2567 50774.640625
    CHARGE=2+
    257.8935242 17.8418960571
    284.0291443 12.8084545135
    319.2088013 14.620434761
    ...
    1950.239746 34.9507827759
    END IONS
    "
    Store the data properly, with title, rt, pepmass, charge of line 2 to line 5. for msms, store in mzs and mzis
    pepmass contains two element, pepmz and pepmzi
    """
    def __init__(self, title = None, rt = None,pepmass = None, charge = None, mzs = [], mzis = []):
        self.title = title
        self.rt = rt
        self.pepmass = pepmass
        self.charge = charge
        self.mzs = mzs
        self.mzis = mzis

    def towrite(self):
        """
        return a string of the original mgf format
        """
        s = "BEGIN IONS"
        s = s + "TITLE=" + self.title+'\n'
        s = s + "RTINSECONDS=" +str(self.rt) + '\n'
        s = s + "PEPMASS=" + str(self.pepmass[0]) + ' ' + str(self.pepmass[1]) + '\n'
        s = s + "CHARGE=" + str(self.charge) + '+\n'
        msms = [str(self.mzs[n]) +str(self.mzis[n]) for n in range(len(self.mzs))]
        msms = '\n'.join(msms)
        s = s + msms + '\n'
        s = s + 'END IONS\n'
        return s
    
def mgfiter(filename):
    """
    given a filename, return a list of class MGF
    """
    fo = open(filename,"r")
    onemgf = None
    line = fo.readline()
    onemgf = None
    while line:
        if line == 'BEGIN IONS\n':
            if onemgf != None:
                onemgf.mzs = np.array(onemgf.mzs)
                onemgf.mzi = np.array(onemgf.mzis)
                yield onemgf
            onemgf = MGF()
            line = fo.readline()
            onemgf.title = line.split('=')[1][:-1]
            line = fo.readline()
            onemgf.rt = float(line.split('=')[1][:-1])
            line = fo.readline()
            pep = line.split('=')[1][:-1].split()
            if len(pep) == 1:
                print(line)
                pep = [pep[0], 0]
            pepmz, pepmzi = pep
            pepmz = float(pepmz)
            pepmzi = float(pepmzi)
            onemgf.pepmass = (pepmz, pepmzi)
            line = fo.readline()
            onemgf.charge = int(line.split('=')[1][:-2])
            line = fo.readline()
            while line:
                if line == 'END IONS\n' or line == 'END IONS':
                    line = fo.readline()
                    break
                else:
                    mz, mzi = line.split()
                    mz = float(mz)
                    mzi = float(mzi)
                    onemgf.mzs.append(mz)
                    onemgf.mzis.append(mzi)
                    line = fo.readline()
    yield onemgf

def mgf2npls(filename):
    '''
    given a filename, return a np array of MGF elements
    '''
    fo = open(filename)
    s = fo.read()
    length = s.count('BEGIN IONS')
    print(s, ' total MGFs')
    del s
    mgfnp = np.array([MGF() for n in range(length)])
    for i, mgf in enumerate(mgfiter(filename)):
        mgfnp[i] = mgf
    return mgfnp

def mgfinfo(filename):
    """
    given a mgf filename (the mgf file read as str), 
    return a str of the whole file,
    lists of mgf start/end/mz/charge/rt
    """
    import re
    s = open(filename).read()
    locs = [(m.start(0), m.end(0)) for m in re.finditer('BEGIN IONS\n.*?END IONS', s,re.DOTALL)]
    start = np.array([ele[0] for ele in locs])
    end = np.array([ele[1] for ele in locs])
    mzs = []
    charges = []
    rts = []
    for num in range(len(locs)):
        loc = locs[num]
        mgf = s[loc[0]:loc[1]]
        mgflines = mgf.split('\n')
        mz = float(mgflines[3].split('=')[1].split()[0])
        charge = int(mgflines[4].split('=')[1][:-1])
        rt = float(mgflines[2].split('=')[1])
        mzs.append(mz)
        charges.append(charge)
        rts.append(rt)
    mzs = np.array(mzs)
    charges = np.array(charges)
    rts = np.array(rts)
    return s, start, end, mzs, charges, rt
    
    

def mzfilter(filename, targetmzcharges,accuracy = 0.01, outfile = None):
    """
    given a mgf filename, a list of tuple with targetmz and charge, output a file of targeted mgfs.
    if outfile not provided, mgffilename + .Target.
    if target charge == 0, then only match targetmz.
    """
    if outfile is None:
        outfile = filename + '.Target'
    fout = open(outfile, 'w')
    mgf_s, mgf_start, mgf_end, mgf_mzs, mgf_charges, mgf_rt = mgfinfo(filename)
    good_index = []
    for num in range(len(mgf_start)):
        for ele in targetmzcharges:
            targetmz, targetcharge = ele
            if targetcharge == 0:
                if abs(mgf_mzs[num] - targetmz) < accuracy:
                    good_index.append(num)
                    break
            else:
                if mgf_charges[num] == targetcharge:
                    if abs(mgf_mzs[num] - targetmz) < accuracy:
                        good_index.append(num)
                        break
    print(len(good_index), ' targets found')
    for ele in good_index:
        fout.write(mgf_s[mgf_start[ele]: mgf_end[ele]] +'\n')
    fout.close()
    
    
        
        
    
