# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 21:22:43 2017

@author: k
"""

def getMs2List20170713():
    '''
    get the charge, mzs from mgf files
    '''
    import os
    folder = 'D:\\Insects\\ManducaSexta\\2017Ms_fatbody_hemocyte_membrane\\mgf\\'
    files = os.listdir(folder)
    from pyteomics import mgf
    mzs = []
    for file in files:
        MGFs = mgf.read(folder+file)
        for _m in MGFs:
            mzs.append([file, _m['params']['charge'][0],_m['params']['pepmass'][0]])
    fout = open(folder+'2017mz_summary.txt','w')
    for e in mzs:
        fout.write('{0}\t{1}\t{2}\n'.format(*e))
    fout.close()
    