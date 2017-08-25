# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 20:43:24 2017

@author: k
"""

def DAP_LysPG_analysis20170723():
    '''
    PGRP1 PGRP3 cut DAP/Lys-PG in the presence of lysozyme
    '''
    import pandas as pd
    dfs = {}#stores the .PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab files in dataframe
    dfs['PGRP1'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1-1402_PGRP1_40hr.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['PGRP3'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1475_PGRP3.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['DAPPG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1525_1_DAP_Lys.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['LysPG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1525_2_Lys_PG.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['PGRP3-LysPG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1525_3_PGRP3_Lys_Zn.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['PGRP1-LysPG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1525_4_PGRP1_Lys_Zn.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['PGRP3-DAPPG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1525_5_PGRP3_DAP_Zn.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)
    dfs['PGRP1-DAPPG'] = pd.read_csv(r"D:\Insects\ManducaSexta\2017Mansi\mzML\1525_6_PGRP1_DAP_Zn.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab",sep = '\t',skiprows = 11)

    #keep columns with v
    columns_keep = ['retention_time','retention_time_window', 'charge', 'mass_to_charge', 'peptide_abundance_study_variable[1]', 'opt_global_spectrum_index']
    for k in dfs:
        dfs[k] = dfs[k].loc[:,columns_keep]
        dfs[k].columns = ['rt','rt_window','charge','mz','intensity','MS1id']
    
    #add a column to indicate the origination of the file
    for k in dfs:
        dfs[k]['file'] = k
    
    #normalize rt of PGRP1, as PGRP1 retention time is  different from others
    dfs['PGRP1']['rt'] = dfs['PGRP1']['rt'] / 65.1 * 114
    
    #use mixture as the base, connect data from PG and protein to the dataframe
    def connect_reactions(df_mix, df_PGRP, df_PG, ppm = 10, rt_error = 600):
        '''
        return a new dataframe based on df_mix
        connect ids in df_PGRP and df_PG to df_mix
        for mz's in the allowed ppm and rt_error range, keep the most abundant one
        also match the charge
        '''
        df_mix = df_mix.sort(columns = 'rt')
        df_PGRP = df_PGRP.sort(columns = 'rt')
        df_PG = df_PG.sort(columns = 'rt')
        index_PGRP = []
        index_PG = []
        for index, row in df_mix.iterrows():
            df_keep = df_PGRP.copy()
            df_keep['keep'] = df_keep.apply(lambda x: abs(x['rt'] - row['rt']) < rt_error, axis = 1)
            df_keep = df_keep[df_keep['keep']]
            if df_keep.shape[0]>=1:
                df_keep['keep'] = df_keep.apply(lambda x: x['charge'] == row['charge'], axis = 1)
                df_keep = df_keep[df_keep['keep']]
                if df_keep.shape[0]>=1:
                    df_keep['keep'] = df_keep.apply(lambda x: abs(x['mz']-row['mz'])/row['mz'] < ppm/1000000, axis = 1)
                    df_keep = df_keep[df_keep['keep']]
            if df_keep.shape[0]>=1:
                max_intensity = max(df_keep.loc[:,'intensity'])
                df_keep['score'] = df_keep.apply(lambda x: (1-abs(x['rt'] - row['rt']) / rt_error) + (1 - abs(x['mz']-row['mz'])/row['mz'] *1000000 / ppm) + (x['intensity']/max_intensity), axis=1)
                df_keep = df_keep.sort(columns = 'score',ascending=False)
                index_PGRP.append(df_keep.index[0])
            else:
                index_PGRP.append(-1)
            
            df_keep = df_PG.copy()
            df_keep['keep'] = df_keep.apply(lambda x: abs(x['rt'] - row['rt']) < rt_error, axis = 1)
            df_keep = df_keep[df_keep['keep']]
            if df_keep.shape[0]>=1:
                df_keep['keep'] = df_keep.apply(lambda x: x['charge'] == row['charge'], axis = 1)
                df_keep = df_keep[df_keep['keep']]
                if df_keep.shape[0]>=1:
                    df_keep['keep'] = df_keep.apply(lambda x: abs(x['mz']-row['mz'])/row['mz'] < ppm/1000000, axis = 1)
                    df_keep = df_keep[df_keep['keep']]
            if df_keep.shape[0]>=1:
                max_intensity = max(df_keep.loc[:,'intensity'])
                df_keep['score'] = df_keep.apply(lambda x: (1-abs(x['rt'] - row['rt']) / rt_error) + (abs(x['mz']-row['mz'])/row['mz'] *1000000 / ppm) + (x['intensity']/max_intensity), axis=1)
                df_keep = df_keep.sort(columns = 'score',ascending=False)
                index_PG.append(df_keep.index[0])
            else:
                index_PG.append(-1)
#            break
        df_PGRP.loc[-1,:] = ''
        df_PG.loc[-1,:] = ''
        df_PGRP_keep = df_PGRP.loc[index_PGRP,:]
        df_PG_keep = df_PG.loc[index_PG,:]
        df_PGRP_keep.columns = ['PGRP_'+e for e in df_PGRP_keep.columns]
        df_PGRP_keep.index = df_mix.index
        df_PG_keep.columns = ['PG_'+e for e in df_PG_keep.columns]
        df_PG_keep.index = df_mix.index
        df_all = pd.concat([df_mix, df_PGRP_keep,df_PG_keep], axis = 1, join_axes = [df_mix.index])
        return df_all
    
    #test
#    df_mix = dfs['PGRP1-DAPPG'].iloc[range(1000),:]
#    df_PGRP = dfs['PGRP1']
#    df_PG = dfs['DAPPG']
#    ppm = 10
#    rt_error = 600
#    print(len([index_PGRP[i] for i in range(len(index_PGRP)) if index_PGRP[i] != df_PGRP.index[i]]))
    
    df_PGRP1DAP = connect_reactions(dfs['PGRP1-DAPPG'],dfs['PGRP1'],dfs['DAPPG'])
    df_PGRP1DAP.to_csv('20170724PGRP1DAP.csv')
    
    df_PGRP1Lys = connect_reactions(dfs['PGRP1-LysPG'],dfs['PGRP1'],dfs['LysPG'])
    df_PGRP1Lys.to_csv('20170724PGRP1Lys.csv')
    
    df_PGRP3DAP = connect_reactions(dfs['PGRP3-DAPPG'],dfs['PGRP3'],dfs['DAPPG'])
    df_PGRP3DAP.to_csv('20170724PGRP3DAP.csv')
    
    df_PGRP3Lys = connect_reactions(dfs['PGRP3-LysPG'],dfs['PGRP3'],dfs['LysPG'])
    df_PGRP3Lys.to_csv('20170724PGRP3Lys.csv')
    
    #calculate possible mz's for DAPPG and LysPG
    
    
def DAP_LysPG_ReanalysisAlicaData20170726():
    '''
    re-analysis of Alicia's data with the new methods
    '''
    folder = 'D:\\Insects\\ManducaSexta\\20160212AliciaMSMS\\20170726mzML\\'
    import os
    import pandas as pd
    files = os.listdir(folder)
    files = [e for e in files if '.mzTab' in e]
    files.sort(key=lambda x:int(x.split('-')[0]))
    dfs ={}
    dfs['PGRP1'] = files[0]
    dfs['DAPPG'] = files[1]
    dfs['LysPG'] = files[2]
    dfs['PGRP1-DAPPG'] = files[3]
    dfs['PGRP1-LysPG'] = files[4]
    dfs['PGRP2'] = files[7]
    dfs['PGRP2-DAPPG'] = files[8]
    dfs['PGRP2-LysPG'] = files[9]
    
    columns_keep = ['retention_time','retention_time_window', 'charge', 'mass_to_charge', 'peptide_abundance_study_variable[1]', 'opt_global_spectrum_index']
    for _k in dfs:
        _df = pd.read_csv(folder+dfs[_k],sep = '\t',skiprows = 11)
        _df = _df.loc[:,columns_keep]
        _df.columns = ['rt','rt_window','charge','mz','intensity','MS1id']
        dfs[_k] =_df
        print(_k,dfs[_k].shape)
    
    #run connect_reactions
    def connect_reactions(df_mix, df_PGRP, df_PG, ppm = 10, rt_error = 600):
        return 'use funtion in previous code'
    
    dfPGRP1_DAPPG = connect_reactions(dfs['PGRP1-DAPPG'], dfs['PGRP1'], dfs['DAPPG'], 20)
    dfPGRP1_LysPG = connect_reactions(dfs['PGRP1-LysPG'], dfs['PGRP1'], dfs['LysPG'], 20)
    dfPGRP2_DAPPG = connect_reactions(dfs['PGRP2-DAPPG'], dfs['PGRP2'], dfs['DAPPG'], 20)
    dfPGRP2_LysPG = connect_reactions(dfs['PGRP2-LysPG'], dfs['PGRP2'], dfs['LysPG'], 20)
    
    dfPGRP1_DAPPG.to_csv(folder+'../20170727PGRP1_DAPPG.csv')
    dfPGRP1_LysPG.to_csv(folder+'../20170727PGRP1_LysPG.csv')
    dfPGRP2_DAPPG.to_csv(folder+'../20170727PGRP2_DAPPG.csv')
    dfPGRP2_LysPG.to_csv(folder+'../20170727PGRP2_LysPG.csv')
    
    
def MansiMSMSPGRP13cuttingReAnalysis20170803():
    '''
    re-analysis of all Mansi's MS data
    '''
    dcFiles = {}
    _fileinfo = r'''empty	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_S1dil.raw
DAPPGmono	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_MPP-DAP.raw
LysPGmono	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_MPP-Lys.raw
PGRP3	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475PGRP3r.raw
PGRP3+DAPPGmono	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_PGRP3plusMPP_DAP.raw
PGRP3+DAPPGmono+Zn	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_PGRP3plusMPP_DAPplusZnSO4.raw
PGRP3+LysPGmono	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_PGRP3plusMPP_Lys.raw
PGRP3+LysPGmono+Zn	D:\Insects\ManducaSexta\2017Mansi\20161128Mansi\raw\1475_PGRP3plusMPP_LysplusZnSO4.raw
PGRP3+DAPPG+Zn	D:\Insects\ManducaSexta\2017Mansi\20170302MansiPGRP_EcPG\1485PGNPGRP3Zn2.raw
DAPPG	D:\Insects\ManducaSexta\2017Mansi\20170302MansiPGRP_EcPG\1485PGN2.raw
DAPPG+Lysozyme	D:\Insects\ManducaSexta\2017Mansi\20170626Mansi\raw\1525_1_DAP_Lys.raw
LysPG+Lysozyme	D:\Insects\ManducaSexta\2017Mansi\20170626Mansi\raw\1525_2_Lys_PG.raw
PGRP3+LysPG+Lysozyme	D:\Insects\ManducaSexta\2017Mansi\20170626Mansi\raw\1525_3_PGRP3_Lys_Zn.raw
PGRP1+LysPG+Lysozyme	D:\Insects\ManducaSexta\2017Mansi\20170626Mansi\raw\1525_4_PGRP1_Lys_Zn.raw
PGRP3+DAPPG+Lysozyme	D:\Insects\ManducaSexta\2017Mansi\20170626Mansi\raw\1525_5_PGRP3_DAP_Zn.raw
PGRP1+DAPPG+Lysozyme	D:\Insects\ManducaSexta\2017Mansi\20170626Mansi\raw\1525_6_PGRP1_DAP_Zn.raw'''
    for _e in _fileinfo.split('\n'):
        _e1,_e2 = _e.split('\t')
        dcFiles[_e1] = _e2
    
    #write scripts to covert raw files to mzML format
    def save2scripts(txt):
        fout = open(r"D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170803scriptsReAnalyze16files.txt",'a')
        fout.write(txt+'\n')
        fout.close()
    save2scripts('#convert 2 mzML format')
    for _k in dcFiles:
        save2scripts(r'"C:\P\ProteoWizard 3.0.10730\msconvert.exe" -o "D:\Insects\ManducaSexta\2017Mansi\20170803mzML" "{0}"'.format(dcFiles[_k]))
    
    #extract MS1 only from mzML files, keep only data from 15min to 75min
    save2scripts('#extract MS1 only from mzML files, keep only data from 15min to 75min')
    for _k in dcFiles:
        _f = dcFiles[_k].split('\\')[-1].split('.')[0]
        save2scripts(r'"C:\P\OpenMS-2.1.0\bin\FileFilter.exe" -in "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.mzML" -out "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.MS1.mzML" -peak_options:level 1 -rt 900:4500 '.format(_f))
        
    #Extract peaks
    save2scripts('\n#extract peaks')
    for _k in dcFiles:
        _f = dcFiles[_k].split('\\')[-1].split('.')[0]
        save2scripts(r'"C:\P\OpenMS-2.1.0\bin\PeakPickerHiRes.exe"  -in "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.MS1.mzML" -out "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.PeakPickerHiRes.mzML" -threads 4  -algorithm:ms_levels 1  -algorithm:SignalToNoise:bin_count 3 -algorithm:SignalToNoise:min_required_elements 7 -algorithm:signal_to_noise 0 '.format(_f))
    
    
    #find features
    save2scripts('\n#find features')
    for _k in dcFiles:
        _f = dcFiles[_k].split('\\')[-1].split('.')[0]
        save2scripts(r'"C:\P\OpenMS-2.1.0\bin\FeatureFinderCentroided.exe"   -in "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.PeakPickerHiRes.mzML" -out "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.PeakPickerHiRes.FeatureFinderCentroided.featureXML" -threads 4 -algorithm:isotopic_pattern:charge_high 4 -algorithm:fit:max_iterations 1000 -algorithm:mass_trace:mz_tolerance 0.004  -algorithm:mass_trace:min_spectra 10 -algorithm:mass_trace:max_missing 0 -algorithm:feature:min_rt_span 0.1 -algorithm:feature:max_rt_span 10000000 -algorithm:feature:rt_shape symetric -algorithm:isotopic_pattern:mz_tolerance 0.02 -algorithm:isotopic_pattern:optional_fit_improvement 100 '.format(_f))
    
    #5.	Extract data use MzTabExporter
    save2scripts('\n#5.	Extract data use MzTabExporter')
    for _k in dcFiles:
        _f = dcFiles[_k].split('\\')[-1].split('.')[0]
        save2scripts(r'"C:\P\OpenMS-2.1.0\bin\MzTabExporter.exe" -in "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.PeakPickerHiRes.FeatureFinderCentroided.featureXML" -out "D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab"'.format(_f))
    
    # Analyize monomers of DAP or Lys PG
    ## 
    import pandas as pd
    columns_keep = ['retention_time', 'charge', 'mass_to_charge', 'peptide_abundance_study_variable[1]']
    dcDF = {}
    for _k in dcFiles:
        _f = dcFiles[_k].split('\\')[-1].split('.')[0]
        _f = r'D:\Insects\ManducaSexta\2017Mansi\20170803mzML'+'\\'+_f+'.PeakPickerHiRes.FeatureFinderCentroided.featureXML.mzTab'
        _df = pd.read_csv(_f,sep = '\t',skiprows = 11)
        _df = _df.loc[:,columns_keep]
        _df.columns = ['rt','charge','mz','intensity']
        dcDF[_k] =_df
        print(_k,dcDF[_k].shape)
    
    def connect_reactions(df1, df2, df2_label = 'df2', ppm = 10, rt_error = 600):
        '''
        return a new dataframe based on df1
        connect ids in df2 to df1
        both df1 and df2 have rt, charge, mz, intensity. in the returned df, change id of df2 to 'df2rt','df2charge'...
        for mz's in the allowed ppm and rt_error range, keep the most abundant one
        also match the charge. change
        '''
        df1 = df1.copy()
        df2 = df2.copy()
        index_keep = []
        def f1(x1,x2,x3):
            return abs(x1-x2)<x3
        df2_filter = df2.copy()
        for index, row in df1.iterrows():
            df2_filter1 = df2_filter[df2_filter.loc[:,'rt'].between(row['rt']-rt_error, row['rt']+rt_error)]
            if df2_filter1.shape[0]>=1:
                df2_filter1 = df2_filter1[df2_filter1.loc[:,'charge'] == row['charge']]
                if df2_filter1.shape[0]>=1:
                    df2_filter1 = df2_filter1[df2_filter1.loc[:,'mz'].between((1-ppm/1000000)*row['mz'], (1+ppm/1000000)*row['mz'])]
            if df2_filter1.shape[0]>=1:
                max_intensity = max(df2_filter1.loc[:,'intensity'])
                df2_filter1.loc[:,'score'] = df2_filter1.apply(lambda x: (1-abs(x['rt'] - row['rt']) / rt_error) + (1 - abs(x['mz']-row['mz'])/row['mz'] *1000000 / ppm) + (x['intensity']/max_intensity), axis=1)
                df2_filter1 = df2_filter1.sort(columns = 'score',ascending=False)
                index_keep.append(df2_filter1.index[0])
            else:
                index_keep.append(-1)
            
        df2.loc[-1,:] = ''
        df2_keep = df2.loc[index_keep,:]
        df2_keep.columns = [df2_label+'_' + e for e in df2_keep.columns]
        df2_keep.index = df1.index
        df_all = pd.concat([df1, df2_keep], axis = 1, join_axes = [df1.index])
        return df_all
    
    #keys: 'empty', 'DAPPGmono', 'LysPGmono', 'PGRP3', 'PGRP3+DAPPGmono', 'PGRP3+DAPPGmono+Zn', 'PGRP3+LysPGmono', 'PGRP3+LysPGmono+Zn', 'PGRP3+DAPPG+Zn', 'DAPPG', 'DAPPG+Lysozyme', 'LysPG+Lysozyme', 'PGRP3+LysPG+Lysozyme', 'PGRP1+LysPG+Lysozyme', 'PGRP3+DAPPG+Lysozyme', 'PGRP1+DAPPG+Lysozyme'
    
    #connect 'PGRP3+DAPPGmono+Zn' 'PGRP3+DAPPGmono' 'PGRP3' 'DAPPGmono' 'empty'
    df_PGRP3_DAPPGmono = dcDF['PGRP3+DAPPGmono+Zn'].copy()
    df_PGRP3_DAPPGmono = connect_reactions(df_PGRP3_DAPPGmono,dcDF['PGRP3+DAPPGmono'],df2_label ='PGRP3+DAPPGmono')
    df_PGRP3_DAPPGmono = connect_reactions(df_PGRP3_DAPPGmono,dcDF['PGRP3'],df2_label ='PGRP3')
    df_PGRP3_DAPPGmono = connect_reactions(df_PGRP3_DAPPGmono,dcDF['DAPPGmono'],df2_label ='DAPPGmono')
    df_PGRP3_DAPPGmono = connect_reactions(df_PGRP3_DAPPGmono,dcDF['empty'],df2_label ='empty')
    df_PGRP3_DAPPGmono.to_csv(r'D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170811MansiPGRP3_DAPPGmono.csv')
    
    df_PGRP3_LysPGmono = dcDF['PGRP3+LysPGmono+Zn'].copy()
    df_PGRP3_LysPGmono = connect_reactions(df_PGRP3_LysPGmono,dcDF['PGRP3+LysPGmono'],df2_label ='PGRP3+LysPGmono')
    df_PGRP3_LysPGmono = connect_reactions(df_PGRP3_LysPGmono,dcDF['PGRP3'],df2_label ='PGRP3')
    df_PGRP3_LysPGmono = connect_reactions(df_PGRP3_LysPGmono,dcDF['LysPGmono'],df2_label ='LysPGmono')
    df_PGRP3_LysPGmono = connect_reactions(df_PGRP3_LysPGmono,dcDF['empty'],df2_label ='empty')
    df_PGRP3_LysPGmono.to_csv(r'D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170811MansiPGRP3_LysPGmono.csv')
    
    # 20170814 the openMS module is not good enough. test new methods.
    #use pyteomics mzml to read files
    from pyteomics import mzml
    import numpy as np
#    lsmzml = mzml.read(r"D:\Insects\ManducaSexta\2017Mansi\20170803mzML\1475_MPP-DAP.PeakPickerHiRes.mzML")
#    lsmzml = np.array(list(lsmzml))
#    scan = lsmzml[0]
    def getTargetmzIntensity(targetmz, scan, ppm=10):
        '''
        scan is a element of mzml from pyteomics
        return 0 if targetmz is not in e_mzml
        '''
        mzs = scan['m/z array']
        mzi = scan['intensity array']
        idx = np.abs(mzs-targetmz).argmin()
        mz = mzs[idx]
        if abs(mz - targetmz)/targetmz < ppm/1000000:
            return mzi[idx]
        return 0
        
    def getTargetmzFromlist(targetmz, lsmzml, ppm=10,timewindow = None):
        '''
        return two numpy array of rt and intensity
        lsmzml is a list of mzml element from pyteomics
        '''
        rts = []
        mzi = []
        if timewindow is None:
            rt_min = 0
            rt_max = float('inf')
        else:
            rt_min,rt_max = timewindow
        for scan in lsmzml:
            rt = scan['scanList']['scan'][0]['scan start time']
            if rt >rt_min and rt < rt_max:
                rts.append(rt)
                mzi.append(getTargetmzIntensity(targetmz,scan,ppm))
        rts = np.array(rts)
        mzi = np.array(mzi)
        return rts, mzi
    
    ## read in all mzml files to dcMzml
    dcMzml = {}
    for _k in dcFiles:
        _f = dcFiles[_k].split('\\')[-1].split('.')[0]
        _mzmlFile = r"D:\Insects\ManducaSexta\2017Mansi\20170803mzML\{0}.PeakPickerHiRes.mzML".format(_f)
        dcMzml[_k] = np.array(list(mzml.read(_mzmlFile)))
    
    ## get mzs for possible products of DAPPGmono
    ls_DAPPGmono = 'C32H55N9O15 C29H49N7O15 C3H8N2O C26H44N6O14 C6H13N3O2 C19H32N4O11 C13H25N5O5 C14H24N2O9 C18H33N7O7 C11H19N1O8 C21H38N8O8 C8H15N1O6 C24H42N8O10'.split()
    ls_LysPGmono = 'C31H55N9O13 C28H49N7O13 C3H8N2O C25H44N6O12 C6H13N3O2 C19H32N4O11 C12H25N5O3 C14H24N2O9 C17H33N7O5 C11H19N1O8 C20H38N8O6 C8H15N1O6 C23H42N8O8'.split()
    from pyteomics import mass
    dc_DAPPGmono = {_k: mass.calculate_mass(formula = _k, charge = 1) for _k in ls_DAPPGmono}
    dc_LysPGmono = {_k: mass.calculate_mass(formula = _k, charge = 1) for _k in ls_LysPGmono}
    ## remove targetmz which is not in the reange of 200 to 1400
    dc_DAPPGmono = {_k:_v for _k,_v in dc_DAPPGmono.items() if _v>200 and _v<1400}
    dc_LysPGmono = {_k:_v for _k,_v in dc_LysPGmono.items() if _v>200 and _v<1400}
    
    f_DAPPGmono = ['empty', 'DAPPGmono', 'PGRP3', 'PGRP3+DAPPGmono', 'PGRP3+DAPPGmono+Zn']
    f_LysPGmono = ['empty', 'LysPGmono', 'PGRP3', 'PGRP3+LysPGmono', 'PGRP3+LysPGmono+Zn']
    dc_PGmono = {} #key is tuple, of target molecule and library, like ('C32H55N9O15','empty'), value is rts and mzi
    for _k in dc_DAPPGmono:
        mz = dc_DAPPGmono[_k]
        for _f in f_DAPPGmono:
            dc_PGmono[(_k,_f)] = getTargetmzFromlist(mz,dcMzml[_f],ppm=10,timewindow=[1000,1400])
    for _k in dc_LysPGmono:
        mz = dc_LysPGmono[_k]
        for _f in f_LysPGmono:
            dc_PGmono[(_k,_f)] = getTargetmzFromlist(mz,dcMzml[_f],ppm=10,timewindow=[1000,1400])
    
    #plot figures
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('20170815MansiDAPPGmono.pdf')
    plt.rcParams["figure.figsize"] = [10,5]
    for _k in dc_DAPPGmono:
        for _f in f_DAPPGmono:
            rts, mzi = dc_PGmono[(_k,_f)]
            plt.plot(rts, mzi, label = str(_f),linewidth = 1)
        plt.legend(loc = 0)
        plt.title('PGRP3 DAPPGmono '+_k+' mz(+1): '+str(dc_DAPPGmono[_k]))
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()
    
    #plot figures

    pp = PdfPages('20170815MansiLysPGmono.pdf')
    plt.rcParams["figure.figsize"] = [10.0,5]
    for _k in dc_LysPGmono:
        for _f in f_LysPGmono:
            rts, mzi = dc_PGmono[(_k,_f)]
            plt.plot(rts, mzi, label = str(_f),linewidth = 1)
        plt.legend(loc = 0)
        plt.title('PGRP3 LysPGmono '+_k+' mz(+1): '+str(dc_LysPGmono[_k]))
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()
    
    
    # 20170815 reaction with polymer PGN
    from pyteomics import mass
    dc = {}
    dc['s'] = mass.Composition(formula = 'C6H13NO5') #sugar
    dc['l'] = mass.Composition(formula = 'C3H6O3') #lactate
    dc['e'] = mass.Composition(formula = 'C2H4O2') #ethanoic acid
    dc['A'] = mass.Composition(formula = 'C3H7NO2') #alanine
    dc['D'] = mass.Composition(formula = 'C5H9NO4') #glutamate
    ## generate a table of possible combinations of elements in dc.keys(), with maximum MW of 3034
    def getStructureFromMW(mw,ppm=10,dc=dc):
        '''
        given a mw, return the structure made from molecules in dc with allowed ppm
        '''
        dcmw = {k:mass.calculate_mass(v) for k,v in dc.items()}
        water_mw = mass.calculate_mass(formula = 'H2O')
        dcmax = {k:int(mw/(dcmw[k]-water_mw))+1 for k in dcmw}
    
    #connect files
    ## connect PGRP3 LysPG
    df_PGRP3_LysPG = dcDF['PGRP3+LysPG+Lysozyme'].copy()
    df_PGRP3_LysPG = connect_reactions(df_PGRP3_LysPG,dcDF['LysPG+Lysozyme'],df2_label ='LysPG+Lysozyme')
    df_PGRP3_LysPG = connect_reactions(df_PGRP3_LysPG,dcDF['PGRP3'],df2_label ='PGRP3')
    df_PGRP3_LysPG = connect_reactions(df_PGRP3_LysPG,dcDF['empty'],df2_label ='empty')
    df_PGRP3_LysPG.to_csv(r'D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170816MansiPGRP3_LysPG.csv')
    ## connect PGRP3 DAPPG
    df_PGRP3_DAPPG = dcDF['PGRP3+DAPPG+Lysozyme'].copy()
    df_PGRP3_DAPPG = connect_reactions(df_PGRP3_DAPPG,dcDF['DAPPG+Lysozyme'],df2_label ='DAPPG+Lysozyme')
    df_PGRP3_DAPPG = connect_reactions(df_PGRP3_DAPPG,dcDF['PGRP3'],df2_label ='PGRP3')
    df_PGRP3_DAPPG = connect_reactions(df_PGRP3_DAPPG,dcDF['PGRP3+DAPPG+Zn'],df2_label ='PGRP3+DAPPG+Zn')
    df_PGRP3_DAPPG = connect_reactions(df_PGRP3_DAPPG,dcDF['DAPPG'],df2_label ='DAPPG')
    df_PGRP3_DAPPG = connect_reactions(df_PGRP3_DAPPG,dcDF['empty'],df2_label ='empty')
    df_PGRP3_DAPPG.to_csv(r'D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170816MansiPGRP3_DAPPG.csv')
    ##connect PGRP1 LysPG
    df_PGRP1_LysPG = dcDF['PGRP1+LysPG+Lysozyme'].copy()
    df_PGRP1_LysPG = connect_reactions(df_PGRP1_LysPG,dcDF['LysPG+Lysozyme'],df2_label ='LysPG+Lysozyme')
    df_PGRP1_LysPG = connect_reactions(df_PGRP1_LysPG,dcDF['empty'],df2_label ='empty')
    df_PGRP1_LysPG.to_csv(r'D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170816MansiPGRP1_LysPG.csv')
    ##connect PGRP1 DAPPG
    df_PGRP1_DAPPG = dcDF['PGRP1+DAPPG+Lysozyme'].copy()
    df_PGRP1_DAPPG = connect_reactions(df_PGRP1_DAPPG,dcDF['DAPPG+Lysozyme'],df2_label ='DAPPG+Lysozyme')
    df_PGRP1_DAPPG = connect_reactions(df_PGRP1_DAPPG,dcDF['DAPPG'],df2_label ='DAPPG')
    df_PGRP1_DAPPG = connect_reactions(df_PGRP1_DAPPG,dcDF['empty'],df2_label ='empty')
    df_PGRP1_DAPPG.to_csv(r'D:\Insects\ManducaSexta\2017Mansi\20170724Result\20170816MansiPGRP1_DAPPG.csv')
    
    #plot target highest enriched peaks of PGRP in reaction with PG+lysozyme
    targets = '229.9566002 252.9725252 270.5177977 270.6849651 292.1920056 302.2142228 313.143101 317.7143497 320.2428439 329.2145042 329.9755298 334.478932 334.7290839 336.7196867 336.8867265 345.2387211 355.7246694 359.4824185 362.9867031 364.9105428 365.077728 383.9160362 384.0831597 384.5816593 389.2535023 391.4943277 405.2731593 412.2737855 420.2587135 427.2638369 433.7810068 437.691256 441.787608 445.1197367 445.635785 460.4973251 462.2875333 468.2952801 476.0665702 483.6461871 504.5745839 546.8610337'.split()
    targets = [float(e) for e in targets]
    
    
    f_PGpoly = ['empty',  'PGRP3',  'PGRP3+DAPPG+Zn', 'DAPPG', 'DAPPG+Lysozyme', 'LysPG+Lysozyme', 'PGRP3+LysPG+Lysozyme', 'PGRP1+LysPG+Lysozyme', 'PGRP3+DAPPG+Lysozyme', 'PGRP1+DAPPG+Lysozyme']
    dc_PGpoly = {} #key is tuple, of target molecule and library, like ('C32H55N9O15','empty'), value is rts and mzi
    for mz in targets:
        for _f in f_PGpoly:
            dc_PGpoly[(str(mz),_f)] = getTargetmzFromlist(mz,dcMzml[_f],ppm=10,timewindow=None)
    
    #plot figures
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('20170816MansiDAPPGpoly.pdf')
    plt.rcParams["figure.figsize"] = [20,10]
    for mz in targets:
        for _f in f_PGpoly:
            rts, mzi = dc_PGpoly[(str(mz),_f)]
            plt.plot(rts, mzi, label = str(_f),linewidth = 1)
        plt.legend(loc = 0)
        plt.title('PGpoly '+' mz: '+str(mz))
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()
    
    dcTargetRange = {}#after initial checking of targets peaks, choose those "real" good peaks to plot. Store proper rt range
    dcTargetRange[229.9566002] = [900,1800]
    dcTargetRange[252.9725252] = [900,1800]
    dcTargetRange[270.5177977] = [2400,3100]
    dcTargetRange[292.1920056] = [2700,3500]
    dcTargetRange[302.2142228] = [1300,1900]
    dcTargetRange[317.7143497] = [2500,3200]
    dcTargetRange[320.2428439] = [900,1200]
    dcTargetRange[329.2145042] = [3100,3800]
    dcTargetRange[329.9755298] = [2400,3100]
    dcTargetRange[334.478932] = [2100,2900]
    dcTargetRange[336.7196867] = [3000,3800]
    dcTargetRange[345.2387211] = [1000,3000]
    dcTargetRange[355.7246694] = [3500,4200]
    dcTargetRange[359.4824185] = [2500,3200]
    dcTargetRange[362.9867031] = [2800,3500]
    dcTargetRange[364.9105428] = [2600,3400]
    dcTargetRange[365.077728] = [2600,3400]
    dcTargetRange[383.9160362] = [3100,3800]
    dcTargetRange[384.0831597] = [3100,3800]
    dcTargetRange[389.2535023] = [2700,3500]
    dcTargetRange[391.4943277] = [3500,4100]
    dcTargetRange[405.2731593] = [2400,3100]
    dcTargetRange[412.2737855] = [2800,3500]
    dcTargetRange[420.2587135] = [1100,1700]
    dcTargetRange[427.2638369] = [3500,4200]
    dcTargetRange[433.7810068] = [2900,3600]
    dcTargetRange[437.691256] = [2600,3400]
    dcTargetRange[441.787608] = [2400,3300]
    dcTargetRange[445.635785] = [2100,2900]
    dcTargetRange[460.4973251] = [3100,3800]
    dcTargetRange[462.2875333] = [3500,4100]
    dcTargetRange[468.2952801] = [3000,4000]
    dcTargetRange[476.0665702] = [2500,3300]
    dcTargetRange[483.6461871] = [2800,3500]
    dcTargetRange[504.5745839] = [3000,3800]
    dcTargetRange[546.8610337] = [2600,3400]
    
    f_PGpoly = ['PGRP3',  'PGRP3+DAPPG+Zn', 'DAPPG', 'DAPPG+Lysozyme', 'LysPG+Lysozyme', 'PGRP3+LysPG+Lysozyme', 'PGRP1+LysPG+Lysozyme', 'PGRP3+DAPPG+Lysozyme', 'PGRP1+DAPPG+Lysozyme']#remove 'empty', since there is always nothing there
    f_PGpolyControl = ['PGRP3', 'DAPPG', 'DAPPG+Lysozyme', 'LysPG+Lysozyme']
    dc_PGpoly = {} #key is tuple, of target molecule and library, like ('C32H55N9O15','empty'), value is rts and mzi
    for mz in dcTargetRange:
        for _f in f_PGpoly:
            dc_PGpoly[(str(mz),_f)] = getTargetmzFromlist(mz,dcMzml[_f],ppm=10,timewindow=dcTargetRange[mz])
    
    #plot figures
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('20170816MansiDAPPGpolyRangeAjusted2.pdf')
    plt.rcParams["figure.figsize"] = [10,5]
    for mz in dcTargetRange:
        for _f in f_PGpoly:
            rts, mzi = dc_PGpoly[(str(mz),_f)]
            if _f in f_PGpolyControl:
                plt.plot(rts, mzi, label = str(_f),linewidth = 2)
            else:
                plt.plot(rts, mzi, label = str(_f),linewidth = 1)
        plt.legend(loc = 0)
        plt.title('PGpoly '+' mz: '+str(mz))
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()
    
    
#plot some target highest enriched peaks of PGRP in reaction with PG+lysozyme
#20170822
    dcTargetRange = {}#after initial checking of targets peaks, choose those "real" good peaks to plot. Store proper rt range
#    dcTargetRange[299.6949029] = [0,5000]
#    dcTargetRange[316.2114] = [1200,2600]
    dcTargetRange[975.3776] = [0,1500]
    
    
    
    f_PGpoly = ['PGRP3',  'PGRP3+DAPPG+Zn', 'DAPPG', 'DAPPG+Lysozyme', 'LysPG+Lysozyme', 'PGRP3+LysPG+Lysozyme', 'PGRP1+LysPG+Lysozyme', 'PGRP3+DAPPG+Lysozyme', 'PGRP1+DAPPG+Lysozyme']#remove 'empty', since there is always nothing there
    f_PGpolyControl = ['PGRP3', 'DAPPG', 'DAPPG+Lysozyme', 'LysPG+Lysozyme']
    dc_PGpoly = {} #key is tuple, of target molecule and library, like ('C32H55N9O15','empty'), value is rts and mzi
    for mz in dcTargetRange:
        for _f in f_PGpoly:
            dc_PGpoly[(str(mz),_f)] = getTargetmzFromlist(mz,dcMzml[_f],ppm=10,timewindow=dcTargetRange[mz])
    
    #plot figures
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('20170822Test.pdf')
    plt.rcParams["figure.figsize"] = [10,5]
    for mz in dcTargetRange:
        for _f in f_PGpoly:
            rts, mzi = dc_PGpoly[(str(mz),_f)]
            if _f in f_PGpolyControl:
                plt.plot(rts, mzi, label = str(_f),linewidth = 2)
            else:
                plt.plot(rts, mzi, label = str(_f),linewidth = 1)
        plt.legend(loc = 0)
        plt.title('PGpoly '+' mz: '+str(mz))
        plt.tight_layout()
        pp.savefig()
        plt.close()
    pp.close()

    
