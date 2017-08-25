# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 20:20:13 2016

@author: k
"""

def cleanMS1toPeaks20161128():
    import ms1Class
    import glob
    folder = 'F:\\Insects\\ManducaSexta\\20161128Mansi\\ms1\\'
    filenames = glob.glob(folder + '*.ms1')
    for filename in filenames:
        ms1Class.ms1ToPeaks(filename,None,10000)

def find_DAP_PG_cuttingEvidence20161202():
    """
    try to find cutting evidence of DAP_PG or LYS_PG
    """
    folder = 'F:\\Insects\\ManducaSexta\\20161128Mansi\\ms1peaks\\'
    dcFile = {} # store file names of ms1.Peak
    dcFile['DAP'] = folder + '1475_MPP-DAP.ms1.Peak'
    dcFile['Lys'] = folder + '1475_MPP-Lys.ms1.Peak'
    dcFile['PGRP3'] = folder + '1475_PGRP3_r.ms1.Peak'
    dcFile['PGRP3+DAP'] = folder + '1475_PGRP3plusMPP_DAP.ms1.Peak'
    dcFile['PGRP3+DAP+Zn'] = folder + '1475_PGRP3plusMPP_DAPplusZnSO4.ms1.Peak'
    dcFile['PGRP3+Lys'] = folder + '1475_PGRP3plusMPP_Lys.ms1.Peak'
    dcFile['PGRP3+Lys+Zn'] = folder + '1475_PGRP3plusMPP_LysplusZnSO4.ms1.Peak'
    dcFile['blank'] = folder + '1475_S1dil.ms1.Peak'
    lsFile = ['PGRP3','DAP','PGRP3+DAP','PGRP3+DAP+Zn','Lys','PGRP3+Lys','PGRP3+Lys+Zn']
    
    import ms1Class
    dcMzs = {}#for keys in dcFile, store Mzs
    dcMzi = {}#for keys in dcFile, store Mzi
    dcRT = {}#for keys in dcFile, store scan_num and rt
    for key in dcFile:
        dcMzs[key], dcMzi[key], dcRT[key] =  ms1Class.readms1Peaks2mzsmzi(dcFile[key])
    
    dcTarget = {} #store target mzs
    dcTarget['DAP, +1'] = 806.3890
    dcTarget['DAP, -1AA, +1'] = 736.3359
    dcTarget['DAP, -2AA, +1'] = 665.2988
    dcTarget['DAP, -3AA, +1'] = 493.2140
    dcTarget['DAP, +3AA, +1'] = 332.1928
    dcTarget['DAP, -4AA, +1'] = 365.1554
    dcTarget['DAP, +4AA, +1'] = 460.2514
    dcTarget['DAP, -5AA, +1'] = 294.1183
    dcTarget['DAP, +5AA, +1'] = 531.2885
    dcTarget['Lys, +1'] = 762.3992
    dcTarget['Lys, -1AA, +1'] = 692.3461
    dcTarget['Lys, -2AA, +1'] = 621.3090
    dcTarget['Lys, -3AA, +1'] = 493.2140
    dcTarget['Lys, +3AA, +1'] = 288.2030
    dcTarget['Lys, -4AA, +1'] = 365.1554
    dcTarget['Lys, +4AA, +1'] = 416.2616
    dcTarget['Lys, -5AA, +1'] = 294.1183
    dcTarget['Lys, +5AA, +1'] = 487.2987
    lsTarget = ['DAP, +1',
                'DAP, -1AA, +1',
                'DAP, -2AA, +1',
                'DAP, -3AA, +1',
                'DAP, +3AA, +1',
                'DAP, -4AA, +1',
                'DAP, +4AA, +1',
                'DAP, -5AA, +1',
                'DAP, +5AA, +1',
                'Lys, +1',
                'Lys, -1AA, +1',
                'Lys, -2AA, +1',
                'Lys, -3AA, +1',
                'Lys, +3AA, +1',
                'Lys, -4AA, +1',
                'Lys, +4AA, +1',
                'Lys, -5AA, +1',
                'Lys, +5AA, +1']
    
    dcTgMzi = {} #store mzis of target mzs in different files
    for target in dcTarget:
        dcTgMzi[target] = {}
        for file in dcFile:
            dcTgMzi[target][file] = ms1Class.targetmzIntensityFromNpmzsmzi(dcTarget[target], dcMzs[file], dcMzi[file])
    
    #plot for each target
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    outfile = '20161202MansiPGRP3CuttingPG_MSMS_HPLC.pdf'
    pp = PdfPages(outfile)
    plt.rcParams["figure.figsize"] = [20.0,10]
    cm = plt.get_cmap('nipy_spectral')
    n_color = 9.0
    for target in lsTarget:
        n = 0
        for file in lsFile:
            n += 1
            x = [ele[1] for ele in dcRT[file]]
            y = dcTgMzi[target][file]
            colorCode = n/n_color
            label =file + ', peak mz: %.4g'%(max([y[_n] for _n in range(len(x)) if x[_n] < 40]))
            plt.subplot(7,1,n)
            plt.plot(x,y,label = label,color = cm(colorCode))
            plt.xlim(0,40)
            plt.legend(loc = 0)
            if n == 1:
                plt.title(target + ": "+str(dcTarget[target]))
        plt.xlabel('retension time')
        plt.ylabel('intensity')
        
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()
    
    #top 20 peaks in each file
    import itertools
    dctopmzs = {}
    import numpy as np

    topn = 20
    for file in lsFile:
        topmzs = []
        _lsmzs = np.array(list(itertools.chain.from_iterable(dcMzs[file])))
        _lsmzi = np.array(list(itertools.chain.from_iterable(dcMzi[file])))
        _lsmzs = _lsmzs[_lsmzi.argsort()[::-1]]
        topmzs.append(_lsmzs[0])
        for ele in _lsmzs:
            if min(abs(np.array(topmzs) - ele)) > 0.003:
                if len(topmzs) < topn:
                    topmzs.append(ele)
                else:
                    break
        dctopmzs[file] = topmzs
    topmzs = []
    for file in dctopmzs:
        topmzs += dctopmzs[file]
    topmzs = [round(float(ele),3) for ele in topmzs]
    topmzs = list(set(topmzs))
    topmzs.sort()
    
    dcTopTgmzi = {}
    for ele in topmzs:
        dcTopTgmzi[ele] = {}
        for file in dcFile:
            dcTopTgmzi[ele][file] = ms1Class.targetmzIntensityFromNpmzsmzi(ele, dcMzs[file], dcMzi[file])
    #plot for each target
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    outfile = '20161202MansiPGRP3CuttingPG_MSMS_HPLC_topMz.pdf'
    pp = PdfPages(outfile)
    plt.rcParams["figure.figsize"] = [20.0,10]
    cm = plt.get_cmap('nipy_spectral')
    n_color = 9.0
    for target in topmzs:
        n = 0
        for file in lsFile:
            n += 1
            x = [ele[1] for ele in dcRT[file]]
            y = dcTopTgmzi[target][file]
            colorCode = n/n_color
            label =file + ', peak mz: %.4g'%(max(y))
            plt.subplot(7,1,n)
            plt.plot(x,y,label = label,color = cm(colorCode))
#            plt.xlim(0,40)
            plt.legend(loc = 0)
            if n == 1:
                plt.title(target)
        plt.xlabel('retension time')
        plt.ylabel('intensity')
        
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()
            
        

def mansi_PGN_PGRP3_dealwith20170419():
    '''
    '''
    from pyteomics import mass
    molecules = ['C3H7N1O2', 'C8H14N2O5','C28H48N8O13','C15H26N4O8','C21H36N6O10','C10H19N3O5','C33H55N9O16','C18H31N5O9','C36H60N10O17']
    m_weight = {}
    for _m in molecules:
        m_weight[_m] = mass.calculate_mass(formula = _m)
    
    import pandas as pd
    #read in three files
    df_pgn = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1485_PGNonly_2.PeakPickerHiRes.FeatureFinderIsotopeWavelet.featureXML.mzTab",sep = '\t',skiprows = 9)
    df_p3 = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1475_PGRP3_r.PeakPickerHiRes.FeatureFinderIsotopeWavelet.featureXML.mzTab",sep = '\t',skiprows = 9)
    df_mix = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1485_PGN_PGRP3_Zn_2.PeakPickerHiRes.FeatureFinderIsotopeWavelet.featureXML.mzTab",sep = '\t',skiprows = 9)
    #keep only several columns
    columns_keep = ['retention_time','retention_time_window', 'charge', 'mass_to_charge', 'peptide_abundance_study_variable[1]', 'opt_global_spectrum_index','opt_global_spectrum_native_id']
    df_pgn = df_pgn.loc[:,columns_keep]
    df_p3 = df_p3.loc[:,columns_keep]
    df_mix = df_mix.loc[:,columns_keep]
    
    #add a column of molecular weight
    proton = 1.00727646677
    def add_mw(df):
        for i in df.index:
            df.loc[i,'molecular weight'] = (df.loc[i,'mass_to_charge'] - proton)* df.loc[i,'charge']
    add_mw(df_p3)
    add_mw(df_pgn)
    add_mw(df_mix)
    
    #combine all dataframe to one
    _df_pg3 = df_p3.copy()
    _df_pg3.loc[:,'file'] = 'PGRP3'
    _df_pgn = df_pgn.copy()
    _df_pgn.loc[:,'file'] = 'PGN'
    _df_mix = df_mix.copy()
    _df_mix.loc[:,'file'] = 'Mix'
    df_all = pd.concat([_df_pg3,_df_pgn, _df_mix])
    df_all.index = range(len(df_all.index))
    
    #match molecular weight to df_all
    ppn = 100
    for i in df_all.index:
        mw = df_all.loc[i,'molecular weight']
        for m,w in m_weight.items():
            if w*(1-ppn/1000000) <= mw and mw <= w*(1+ppn/1000000):
                df_all.loc[i,'molecule'] = m
                break
    import numpy as np
    df_matched = df_all.dropna()
    
    #find mix unique molecular weights
    ppn = 50
    df_mix_sort = df_mix.sort(columns = ['peptide_abundance_study_variable[1]','molecular weight'],ascending=False)
    df_premix = df_all[df_all.loc[:,'file'] != 'Mix']
    def within(data,up,low):
        return data <=up and data >= low

    for i in df_mix_sort.index:
        mw = df_mix_sort.loc[i,'molecular weight']
        df_mix_sort.loc[i,'found'] = ''
#            for j in df_premix.index:
#                w = df_premix.loc[j,'molecular weight']
#                if w*(1-ppn/1000000) <= mw and mw <= w*(1+ppn/1000000):
#                    df_mix_sort.loc[i,'found'] += df_premix.loc[j,'file']
        ws = df_premix.loc[:,'molecular weight']
        found = ws.apply(within, args = (mw*(1+ppn/1000000), mw*(1-ppn/1000000)))
        df_mix_sort.loc[i,'found'] = ' '.join(df_premix[found].loc[:,'file'])
    
    
    df_mix_sort.to_csv('result.csv')
    
    #
def findMoleculebyMW20170420():
    '''
    '''
    target1 = 2411.490297
    
    m1 = 'C8H15O6N1'
    m2 = 'C11H19O8N1'
    m3 = 'C3H7O2N'
    m4 = 'C5H9O4N1'
    m5 = 'C7H14O4N2'
    
    from pyteomics import mass
    w1 = mass.calculate_mass(formula = m1)
    w2 = mass.calculate_mass(formula = m2)
    w3 = mass.calculate_mass(formula = m3)
    w4 = mass.calculate_mass(formula = m4)
    w5 = mass.calculate_mass(formula = m5)
    water = mass.calculate_mass(formula = 'H2O1')
    
    for a1 in range(int((target1-water)/(w1-water)) + 1):
        for a2 in range(int((target1-water -a1*(w1-water))/(w2-water)) + 1):
            for a3 in range(int((target1-water -a1*(w1-water) - a2*(w2-water))/(w3-water)) + 1):
                for a4 in range(int((target1-water -a1*(w1-water) - a2*(w2-water)- a3*(w3-water))/(w4-water)) + 1):
                    for a5 in range(int((target1-water -a1*(w1-water) - a2*(w2-water)- a3*(w3-water) - a4*(w4-water))/(w5-water)) + 1):
                        if abs(a1 * w1 + a2 * w2 + a3*w3 +a4*w4 + a5*w5 + -target1 -(a1+a2+a3+a4+a5-1) * water) <1:
                            print(a1, a2, a3, a4, a5)
    
def mzMLfileTest20170423():
    '''
    OpenMS not accurate enough 
    '''
#    f_peaks = r"D:\Insects\ManducaSexta\2017Mansi\20170302MansiPGRP_EcPG\1485_PGN_PGRP3_Zn_2\p0\1485_PGN_PGRP3_Zn_2.peaks"
#    fo_peaks = open(f_peaks,'rb')
#    a = fo_peaks.read(8)
#    b = fo_peaks.read(40)
#    c = a+b+fo_peaks.read(1000)
    
#    f_iso = r"D:\Insects\ManducaSexta\2017Mansi\20170302MansiPGRP_EcPG\1485_PGN_PGRP3_Zn_1\p0\1485_PGN_PGRP3_Zn_1.iso"
#    fo_iso = open(f_iso,'rb')
#    a = fo_iso.read()
#    import struct
#    count = struct.unpack('i',a[:4])[0]
#    members = struct.unpack('i'*count,a[4:4+4*count])
#    f_mz = r"D:\Insects\ManducaSexta\2017Mansi\20170302MansiPGRP_EcPG\1485_PGN_PGRP3_Zn_1\p0\1485_PGN_PGRP3_Zn_1.peaks_mzrtinfo"
#    fo_mz = open(f_mz,'rb')
#    b = fo_mz.read()
    
    import pymzml
    files ={}
    files['mix'] = r"D:\Insects\ManducaSexta\2017Mansi\mzML\1485PGNPGRP3Zn2.PeakPickerHiRes.mzML"
    files['P3'] = r"D:\Insects\ManducaSexta\2017Mansi\mzML\1475PGRP3r.PeakPickerHiRes.mzML"
    files['PG'] = r"D:\Insects\ManducaSexta\2017Mansi\mzML\1485PGN2.PeakPickerHiRes.mzML"
    run = {}
    for k in files:
        run[k] = pymzml.run.Reader(files[k])
    def extractPeaks(filename):
        run = pymzml.run.Reader(filename,obo_version='3.71.0')
        run2 = pymzml.run.Writer(filename = filename+'ms1peaks', run= run , overwrite = True)
        for spec in run:
            if spec['ms level'] == 1:
                spec.peaks = spec.centroidedPeaks
                run2.addSpec(spec)
        run2.save()
    for k in files:
        extractPeaks(files[k])
        break

def mansiCheck20170426():
    '''
    re-run functions in mansi_PGN_PGRP3_dealwith20170419
    '''
    import pandas as pd
    dfs = {}
    dfs['PG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1485PGN2.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['P3'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1475PGRP3r.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['mix'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1485PGNPGRP3Zn2.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    columns_keep = ['retention_time','retention_time_window', 'charge', 'mass_to_charge', 'peptide_abundance_study_variable[1]', 'opt_global_spectrum_index']
    for k in dfs:
        dfs[k] = dfs[k].loc[:,columns_keep]
        dfs[k].columns = ['rt','rt_window','charge','mz','intensity','MS1id']
    
    #calculate molecular weight and add mz ppm
    ppm = 10
    for k in dfs:
        dfs[k].loc[:,'mw'] = dfs[k].apply(lambda row:(row['mz'] - 1.00727646677) * row['charge'], axis=1)
        dfs[k].loc[:,'mwl'] = dfs[k].loc[:,'mw'] *(1-dfs[k].loc[:,'charge']*ppm/1000000)
        dfs[k].loc[:,'mwh'] = dfs[k].loc[:,'mw'] *(1+dfs[k].loc[:,'charge']*ppm/1000000)
    #add column showing source file and combine
    for k in dfs:
        dfs[k].loc[:,'file'] = k
    dfall = pd.concat(dfs.values())
    #sort by mw
    dfall = dfall.sort(columns = 'mw')
    dfall.index = range(len(dfall.index))
    
    from pyteomics import mass
    molecules = ['C3H7N1O2', 'C8H14N2O5','C28H48N8O13','C15H26N4O8','C21H36N6O10','C10H19N3O5','C33H55N9O16','C18H31N5O9','C36H60N10O17','C23H41N7O10']
    m_weight = {}
    for _m in molecules:
        m_weight[_m] = mass.calculate_mass(formula = _m)
    mws = list(m_weight.items())
    mws.sort(key=lambda x:x[1])
    for mw in mws:
        print(mw[0],mw[1],mw[1]*(1-10/1000000),mw[1]*(1+10/1000000))
    for mw in mws:
        print(mw[0],mw[1],mw[1]*(1-20/1000000),mw[1]*(1+20/1000000))
    for mw in mws:
        print(mw[0],mw[1],mw[1]*(1-30/1000000),mw[1]*(1+30/1000000))
    for mw in mws:
        print(mw[0],mw[1],mw[1]*(1-40/1000000),mw[1]*(1+40/1000000))
    
    dc = {}
    dc['N'] = mass.Composition(formula = 'C8H15NO6')
    dc['L'] = mass.Composition(formula = 'C3H6O3')
    dc['A'] = mass.Composition(formula = 'C3H7NO2')
    dc['G'] = mass.Composition(formula = 'C5H9NO4')
    dc['D'] = mass.Composition(formula = 'C7H14N2O4')
    
    def getCompositon(s, returnstr = True):
        '''
        given a string of single letter of dc keys, return composition
        remove water molecules
        if returnstr, return formula, otherwise return composition
        '''
        m = mass.Composition()
        for k in s:
            m = m+dc[k]
        m = m - (len(s)-1) * mass.Composition(formula = 'H2O')
        if returnstr:
            f = ''
            for a,n in m.items():
                f = f + a +str(n)
            m = f
        return m
    
    target1 = []
    s1 = 'NNLAGDA'
    for i in range(len(s1)):
        for j in range(i+1,len(s1)+1):
            target1.append(s1[i:j])
    s2 = 'ADGALNN'
    target2 = []
    for i in range(len(s2)):
        for j in range(i+1,len(s2)+1):
            target2.append(s2[i:j])
    import itertools
    targets = list(itertools.product([t for t in target1 if 'D' in t],['A']+[t for t in target2 if 'AD' in t]))
    targets = [''.join(t) for t in targets] +target1
    #keep unique composition
    targets = list(set([''.join(sorted(t)) for t in targets]))
    
    dcFormula = {}
    dcMW = {}
    for t in targets:
        formula = getCompositon(t)
        dcFormula[t] = formula
        dcMW[t] = mass.calculate_mass(formula = formula)
    
    #add column of molecule to dfall
    def getMolecule(row):
        molecules = []
        mws = []
        formulas = []
        for m, mw in dcMW.items():
            if row['mwl'] <= mw <= row['mwh']:
                molecules.append(m)
                mws.append(str(mw))
                formulas.append(dcFormula[m])
        return ' '.join(molecules)
    def getMoleculeMW(row):
        molecules = []
        mws = []
        formulas = []
        for m, mw in dcMW.items():
            if row['mwl'] <= mw <= row['mwh']:
                molecules.append(m)
                mws.append(str(mw))
                formulas.append(dcFormula[m])
        return ' '.join(mws)
    def getMoleculeFormulas(row):
        molecules = []
        mws = []
        formulas = []
        for m, mw in dcMW.items():
            if row['mwl'] <= mw <= row['mwh']:
                molecules.append(m)
                mws.append(str(mw))
                formulas.append(dcFormula[m])
        return ' '.join(formulas)
    dfall.loc[:,'match_m'] = dfall.apply(getMolecule,axis = 1)
    dfall.loc[:,'match_f'] = dfall.apply(getMoleculeFormulas,axis = 1)
    dfall.loc[:,'match_MW'] = dfall.apply(getMoleculeMW,axis = 1)
    
    dfmatched = dfall[dfall.loc[:,'match_m'] != '']
    dfmatched.to_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\2017matched.csv")
    
    
    #find high abundant new molecules in mix
    dfmix = dfs['mix'].copy()
    dfmix = dfmix.sort(columns = 'intensity',ascending = False)
    dfmix = dfmix.iloc[:1000,:]
    import numpy as np
    def getMatchInOtherRun(row, df_other, time_error =600):
        '''
        get values in df_other, the highest intensity within the allowed time error
        '''
        df_other = df_other.sort(columns = 'rt')
        rts = np.array(df_other.loc[:,'rt'])
        locmin = max(0,np.argmin(abs(rts - row['rt'] +time_error)))
        locmax = min(len(rts),np.argmin(abs(rts - row['rt'] -time_error)))
        df_otherTime = df_other.iloc[locmin:locmax,:]
        df_otherTime = df_otherTime.sort(columns = 'mw')
        mws = np.array(df_otherTime.loc[:,'mw'])
        locmin = max(0,np.argmin(abs(mws - row['mw'] +0.5)))
        locmax = min(len(rts),np.argmin(abs(mws - row['mw'] -0.5)))
        df_otherMW = df_otherTime.iloc[locmin:locmax,:]
        df_otherMW = df_otherMW[df_otherTime.apply(lambda r:r['mw'] >= row['mwl'] and r['mw']<=row['mwh'], axis = 1)]
        df_otherMW = df_otherMW.sort(columns = 'intensity',ascending = False)
        if df_otherMW.shape[0] >=1:
            return df_otherMW.iloc[0,:]
    
    df_InPG = dfmix.apply(getMatchInOtherRun,axis = 1, args = (dfs['PG'],600))
    df_InP3 = dfmix.apply(getMatchInOtherRun,axis = 1, args = (dfs['P3'],600))
    
    dfmixFound = pd.concat([dfmix, df_InPG, df_InP3], axis = 1)
    dfmixFound.to_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\2017mixHigh.csv")
    
    
    def predictStructure(mw, charge = 1, ppm = 10):
        '''
        given a mw, predict composition
        '''
        from pyteomics import mass
        dc = {}
        dc['N'] = mass.Composition(formula = 'C8H15NO6')
        dc['L'] = mass.Composition(formula = 'C3H6O3')
        dc['A'] = mass.Composition(formula = 'C3H7NO2')
        dc['G'] = mass.Composition(formula = 'C5H9NO4')
        dc['D'] = mass.Composition(formula = 'C7H14N2O4')
        H = mass.Composition(formula = 'H2O')
        Hw = mass.calculate_mass(composition = H)
        mwl = mw *(1-charge * ppm/1000000)
        mwh = mw *(1+charge * ppm/1000000)
        dcMW = {}
        for k in dc:
            dcMW[k] = mass.calculate_mass(composition = dc[k])
        
            
        items = [(k,mass.calculate_mass(composition = dc[k]),int(mw/(mass.calculate_mass(composition = dc[k]) - Hw)) +1) for k in dc]
        combinations = itertools.product(*[list(range(n[2])) for n in items], list(range(sum(m[2] for m in items))))
        mwl = mw *(1-charge * ppm/1000000)
        mwh = mw *(1+charge * ppm/1000000)
        items.append(('H',-Hw,0))
        for c in combinations:
            if c[-1] <= sum(c[:-1]):
                predictedMW = sum(items[i][1] * c[i] for i in range(6))
                if mwl < predictedMW and predictedMW < mwh:
                    comp = ''
                    for i in range(6):
                        comp += items[i][0]+str(c[i])
                    print(mw,charge, comp)
        
        target = 1333.891745
        import itertools
        combinations = itertools.product(list(range(42,62)), list(range(55,115)), list(range(5,16)), list(range(18,38)), list(range(3)), list(range(4)))
        from pyteomics import mass
        for c in combinations:
            formula = 'C{0}H{1}N{2}O{3}P{4}S{5}'.format(*c)
            m = mass.calculate_mass(formula = formula)
            if target * (1 - 10/1000000) < m and m < target * (1 + 10/1000000):
                print(formula)

def getAllPossibleTargetSequences20170427():
    '''
    monomer of E. coli PG: NNLAGDA, H can be linked to N, to show that the sugar can loss one more water molecule.
    '''
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
    
    from pyteomics import mass
    dc = {}
    dc['N'] = mass.Composition(formula = 'C8H15NO6')
    dc['L'] = mass.Composition(formula = 'C3H6O3')
    dc['A'] = mass.Composition(formula = 'C3H7NO2')
    dc['G'] = mass.Composition(formula = 'C5H9NO4')
    dc['D'] = mass.Composition(formula = 'C7H14N2O4')
    dc['H'] = mass.Composition()
    
    def getCompositionAndMW(s,waterCount):
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
    
#    def getsubCompositonFromPGNnetwork(PGN, unordered = True):
#        '''
#        given a string of single letter of dc keys, return composition
#        remove water molecules
#        if returnstr, return formula, otherwise return composition
#        return a list of tuple, with ADGLNH compostion and water to remove
#        '''
#        pgn_subs = PGN.get_independent_map()
#        links = PGN.get_links()
#        ms = []
#        for pgn_sub in pgn_subs:
#            s = ''
#            watercount = 0
#            for link in links:
#                if link[0] in pgn_sub or link[1] in pgn_sub:
#                    watercount += 1
#            for knode_key in pgn_sub:
#                s += PGN.knodes[knode_key] #get molecule in str
#            if unordered:
#                s = ''.join(sorted(s))
#            ms.append((s,watercount))
#        return ms
#    
#    def getPGNallsubProducts(PGN, unordered = True, unique = True):
#        '''
#        use getsubCompositonFromPGNnetwork, get all possible subproducts
#        '''
#        products = set()
#        import copy
#        import itertools
#        links = PGN.get_links()
#        for i in range(len(links)):
#            for to_removes in itertools.combinations(links,i):
#                newpgn = copy.deepcopy(PGN)
#                for to_remove in to_removes:
#                    newpgn.remove_links(to_remove[0],to_remove[1])
#                    newproducts = getsubCompositonFromPGNnetwork(newpgn)
#                    products.update(newproducts)
#        return products
    
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
    
    import pandas as pd
    df_mixHigh =pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\2017mixHigh.csv",index_col  = 0)
    products = getPGNallProducts(PGN4)
    productsFormulaMW = []#list of tuple, with string description, formula, and molecular weight
    for product in products:
        _s, _waterCount = product
        _f,_mw = getCompositionAndMW(_s,_waterCount)
        productsFormulaMW.append((_s+'-%dH2O'%_waterCount,_f,_mw))
    productsFormulaMW.sort(key = lambda x:x[2])
    import numpy as np
    productsMW = np.array([_e[2] for _e in productsFormulaMW])
    
    def getMolecule(row):
        molecules = []
        mws = []
        formulas = []
        loc_min = max(0,np.argmin(abs(productsMW -row['mwl']))-1)
        loc_max = min(np.argmin(abs(productsMW-row['mwh'])) +2, len(productsMW))
        for m, f, mw in productsFormulaMW[loc_min:loc_max]:
            if row['mwl'] <= mw <= row['mwh']:
                molecules.append(m)
                mws.append(str(mw))
                formulas.append(f)
        return ' '.join(molecules), ' '.join(formulas), ' '.join(mws)
    
    dfmatch = df_mixHigh.apply(getMolecule, axis = 1)
    for _n in dfmatch.index:
        _m, _f, _mw = dfmatch[_n]
        df_mixHigh.loc[_n,'match_m'] = _m
        df_mixHigh.loc[_n,'match_f'] = _f
        df_mixHigh.loc[_n,'match_mw'] = _mw
    
    df_all = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\2017All3_molecularWeights_intensity.csv", index_col = 0)
    df_all.loc[:,'mwl'] = df_all.apply(lambda row:row['mw'] * (1-row['charge']*10/1000000), axis = 1)
    df_all.loc[:,'mwh'] = df_all.apply(lambda row:row['mw'] * (1+row['charge']*10/1000000), axis = 1)
    dfmatchall = df_all.apply(getMolecule, axis = 1)
    df_all.loc[:,'match_m'] = [e[0] for e in dfmatchall]
    df_all.loc[:,'match_f'] = [e[1] for e in dfmatchall]
    df_all.loc[:,'match_mw'] = [e[2] for e in dfmatchall]