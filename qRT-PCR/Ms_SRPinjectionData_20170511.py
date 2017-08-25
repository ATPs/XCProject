# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:05:15 2017

@author: k
"""

def get_amplication_efficiency_20170511():
    '''
    qRT-PCR, quantification amplification result, contains RFU (relative fluorescence units) of each well in each cycle
    use two points, round threshold, like 200 RFU, to estimate amplication efficiency for each well.
    use the median or average value as the amplication efficiency
    '''
    #test
    import pandas as pd
    file_test = r"D:\Lab\works\20170511qPCR\2017-05-02 rps3 CIFGH heat-treatment 6h 24h with replicates -  Quantification Amplification Results_SYBR.csv"
    df_test = pd.read_csv(file_test)
    
    threshold = 150
    def get_efficiency(column, threshold = 150):
        '''
        return efficiency, RFU greater than threshold /RFU smaller than threshold
        '''
        for n in range(len(column)-1):
            if column[n+1]>=threshold and column[n] <= threshold:
                return column[n+1]/column[n]
#        print('threshold too high!')
        return None
    
    efficiency_test = df_test.loc[:,df_test.columns[2:]].apply(get_efficiency,axis = 0,args = ([threshold]))
    print('mean:',efficiency_test.mean(), '\nmedian:',efficiency_test.median())
    
    
    #calculate efficiency for each sample
    folder = 'D:\\Lab\\works\\20170511qPCR\\'
    import os
    files = os.listdir(folder)
    file_sets = []
    for _f in files:
        _fs = _f.split(' - ')[0]
        if _fs not in file_sets and '.csv' not in _fs:
            file_sets.append(_fs)
    
    targets = set()
    sample = set()
    for _f in file_sets:
        _fe = folder + _f+' -  End Point Results_SYBR.csv'
        _fo = open(_fe).readlines()
        for _l in _fo[1:]:
            _t = _l.split(',')[3]
            _s = _l.split(',')[5]
            targets.add(_t)
            sample.add(_s)
    
    targets.remove('')
    sample.remove('')
    #readin all amplication data
    dcdf_amplication = {}
    dcdf_endpoint = {}
    for _f in file_sets:
        dcdf_amplication[_f] = pd.read_csv(folder + _f+' -  Quantification Amplification Results_SYBR.csv')
        dcdf_endpoint[_f] = pd.read_csv(folder + _f+' -  End Point Results_SYBR.csv',index_col = 1)
    #calculate efficiency for each sample in dataframe
    def endpointindexchange(e):
        '''
        change A01 to A1
        '''
        if e[1] == '0' and len(e) == 3:
            return e[0]+e[-1]
        return e
    dcdf_efficiency = {}
    for _f in file_sets:
        dcdf_efficiency[_f] = dcdf_amplication[_f].loc[:,df_test.columns[2:]].apply(get_efficiency,axis = 0,args = ([threshold]))
        dcdf_efficiency[_f] = dcdf_efficiency[_f].drop([e for e in dcdf_efficiency[_f].index if e not in dcdf_endpoint[_f].index])
        dcdf_endpoint[_f].index = [endpointindexchange(e) for e in dcdf_endpoint[_f].index]
        dcdf_efficiency[_f].index = [dcdf_endpoint[_f].loc[e,'Target'] for e in dcdf_efficiency[_f].index]
    
    efficiency_all = pd.concat(dcdf_efficiency.values())
    efficiency_all_noNA = efficiency_all.dropna()
    efficiency_all_filter1 = efficiency_all_noNA[efficiency_all_noNA>1.5]
    efficiency_all_filter2 = efficiency_all_filter1[efficiency_all_filter1<2.5]
    
    dc_efficiency = {}
    for _target in set(efficiency_all_filter2.index):
        dc_efficiency[_target] = efficiency_all_filter2[_target]
        print(_target, dc_efficiency[_target].mean()-1, dc_efficiency[_target].median()-1)
    
    #mean middle 50% values
    for _target in dc_efficiency:
        _s = dc_efficiency[_target]
        _low = _s.quantile(q=0.25)
        _high = _s.quantile(q=0.75)
        _s = _s[_s>_low]
        _s = _s[_s<_high]
        print(_target, _s.mean()-1)
    

def calculation_of_expression20170512():
    '''
    sum up Ct values, caluclate if needed
    '''
    import pandas as pd
    samples = 'CF IF CH IH Cuticle Midgut MT Muscle Nerve SG Trachea FB-5th-D2 FB-5th-D5 FB-wandering FB-prepupa-D2 FB-prepupa-D5 FB-pupa-D1 FB-heat-6h-1 FB-heat-6h-2 FB-heat-6h-3 FB-heat-24h FB-PBS-6h-1 FB-PBS-6h-2 FB-PBS-6h-3 FB-PBS-24h-1 FB-PBS-24h-2 FB-PBS-24h-3 FB-mix-6h-1 FB-mix-6h-2 FB-mix-6h-3 FB-mix-24h-1 FB-mix-24h-2 FB-mix-24h-3 FB-srp1-6h-1 FB-srp1-6h-2 FB-srp1-6h-3 FB-srp1-24h FB-srp2-6h-1 FB-srp2-6h-2 FB-srp2-6h-3 FB-srp2-24h FB-srp3-6h-1 FB-srp3-6h-2 FB-srp3-6h-3 FB-srp3-24h FB-srp4-6h-1 FB-srp4-6h-2 FB-srp4-6h-3 FB-srp4-24h FB-srp6-6h-1 FB-srp6-6h-2 FB-srp6-6h-3 FB-srp6-24h FB-uENF1-6h-1 FB-uENF1-6h-2 FB-uENF1-6h-3 FB-uENF1-24h HC-5th-D2 HC-5th-D5 HC-wandering HC-prepupa-D2 HC-prepupa-D5 HC-pupa-D1 HC-heat-6h-1 HC-heat-6h-2 HC-heat-6h-3 HC-heat-24h HC-PBS-6h-1 HC-PBS-6h-2 HC-PBS-6h-3 HC-PBS-24h-1 HC-PBS-24h-2 HC-PBS-24h-3 HC-mix-6h-1 HC-mix-6h-2 HC-mix-6h-3 HC-mix-24h-1 HC-mix-24h-2 HC-mix-24h-3 HC-srp1-6h-1 HC-srp1-6h-2 HC-srp1-6h-3 HC-srp1-24h HC-srp2-6h-1 HC-srp2-6h-2 HC-srp2-6h-3 HC-srp2-24h HC-srp3-6h-1 HC-srp3-6h-2 HC-srp3-6h-3 HC-srp3-24h HC-srp4-6h-1 HC-srp4-6h-2 HC-srp4-6h-3 HC-srp4-24h HC-srp6-6h-1 HC-srp6-6h-2 HC-srp6-6h-3 HC-srp6-24h HC-uENF1-6h-1 HC-uENF1-6h-2 HC-uENF1-6h-3 HC-uENF1-24h MG-5th-D2 MG-5th-D5 MG-wandering MG-prepupa-D2 MG-prepupa-D5 MG-pupa-D1 MG-heat-6h-1 MG-heat-6h-2 MG-heat-6h-3 MG-heat-24h MG-PBS-6h-1 MG-PBS-6h-2 MG-PBS-6h-3 MG-PBS-24h-1 MG-PBS-24h-2 MG-PBS-24h-3 MG-mix-6h-1 MG-mix-6h-2 MG-mix-6h-3 MG-mix-24h-1 MG-mix-24h-2 MG-mix-24h-3 MG-srp1-6h-1 MG-srp1-6h-2 MG-srp1-6h-3 MG-srp1-24h MG-srp2-6h-1 MG-srp2-6h-2 MG-srp2-6h-3 MG-srp2-24h MG-srp3-6h-1 MG-srp3-6h-2 MG-srp3-6h-3 MG-srp3-24h MG-srp4-6h-1 MG-srp4-6h-2 MG-srp4-6h-3 MG-srp4-24h MG-srp6-6h-1 MG-srp6-6h-2 MG-srp6-6h-3 MG-srp6-24h MG-uENF1-6h-1 MG-uENF1-6h-2 MG-uENF1-6h-3 MG-uENF1-24h'''.split()
    targets = 'rpS3 SRP1 SRP2 SRP3 SRP4 SRP6 uENF1 AFP4 attacin1 PP cecropin6 gloverin lebocinD lysozyme MBP moricine PGRP5 serpin6'.split()
    folder = 'D:\\Lab\\works\\20170511qPCR\\'
    import os
    files = os.listdir(folder)
    file_sets = []
    for _f in files:
        _fs = _f.split(' - ')[0]
        if _fs not in file_sets and '.csv' not in _fs:
            file_sets.append(_fs)
    
    targets2 = []
    for _e in targets:
        targets2.append(_e)
        targets2.append(_e)
    
    import numpy as np
    from collections import defaultdict
    df = pd.DataFrame(data = np.NaN,columns = targets2, index = samples)
    for _f in file_sets:
        _df = pd.read_csv(folder + _f+' -  Quantification Cq Results_0.csv')
        _dc = defaultdict(list)
        for _r, _row in _df.iterrows():
            _rr = _row['Sample']
            _rc = _row['Target']
            _cq = _row['Cq']
            _dc[(_rr,_rc)].append(_cq)
        for _k in _dc:
                if _k[0] in samples and _k[1] in targets:
                    df.loc[_k[0], _k[1]] = _dc[_k]
    
    