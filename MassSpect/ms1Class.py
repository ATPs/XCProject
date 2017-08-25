# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 13:57:20 2016

@author: k


store in a class MS1
"""
import numpy as np

class MS1fromFMSConvertGUI(object):
    """
    MS1fromFMSConvertGUI stores scan_num, rt, bpi, bpm, tic, mz list and intensity list
    ms1 data like:
    H	CreationDate Wed Feb 24 15:36:27 2016
    H	Extractor	ProteoWizard
    H	Extractor version	Xcalibur
    H	Source file	1-1402_PGRP1_40hr.RAW
    S	1072	1072
    I	RTime	15.00011
    I	BPI	223564.5
    I	BPM	149.0231
    I	TIC	6592032
    130.0008 0
    130.0011 0
    130.0013 0
    130.0016 0
    130.0842 0
    130.0845 0
    130.0848 0
    130.0851 0
    130.0853 2270.266
    130.0856 5149.655
    130.0859 6727.027
    130.0862 5799.265
    130.0865 2918.76
    130.0867 108.4825
    130.087 0
    """
    def __init__(self,scan_num = None, rt=None, bpi=None, bpm=None, tic=None):
        self.scan_num = scan_num
        self.rt = rt
        self.bpi = bpi
        self.bpm = bpm
        self.tic = tic
        mzs =[]
        mzi =[]
        self.mzs =mzs
        self.mzi = mzi
    def peaks(self, min_intensity = 10000,peak_error = 0.2):
        """
        return a list of peaks by finding local maximun. for example, 
        130.0851 0
        130.0853 2270.266
        130.0856 2000
        130.0860 3000
        130.0866 0
        
        peak is [130.0853, 130.0860] if min_intensity<2000, 
        for targetmzs, if two difference between two targetmz is less than error, keep the larger one (the one with higher intensity)
        """
        if len(self.mzi) == 0:
            return []
        xval = np.asarray(self.mzs)
        yval = np.asarray(self.mzi)
        gradient = np.diff(yval)
        maxima = np.diff((gradient >0).view(np.int8))
        maxindex = np.concatenate((([0],) if gradient[0] < 0 else ())+(np.where(maxima == -1)[0] + 1,)+(([len(yval)-1],) if gradient[-1] > 0 else ()))
        ypeak = yval[maxindex]
        maxindex_above_min = maxindex[np.where(ypeak>min_intensity)]
        targetx = xval[maxindex_above_min]
        if len(targetx) == 0:
            return None
        if peak_error ==0:
            return targetx
        else:
            targety = yval[maxindex_above_min]
            targetx2 =[]
            targety2 =[]
            targetx2.append(targetx[0])
            targety2.append(targety[0])
            for num in range(1,len(targetx)):
                if targetx[num] - targetx2[-1] < peak_error:
                    if targety[num] > targety2[-1]:
                        targetx2[-1] = targetx[num]
                        targety2[-1] = targety[num]
                else:
                    targetx2.append(targetx[num])
                    targety2.append(targety[num])
            return targetx2
    
    def simplePeaks(self, min_intensity = 1000000):
        """
        function similar to peaks. return peaks and its intensity directly, the only filter is min_intensity.
        """
        if len(self.mzi) == 0:
            return [],[]
        xval = np.asarray(self.mzs)
        yval = np.asarray(self.mzi)
        gradient = np.diff(yval)
        maxima = np.diff((gradient >0).view(np.int8))
        maxindex = np.concatenate((([0],) if gradient[0] < 0 else ())+(np.where(maxima == -1)[0] + 1,)+(([len(yval)-1],) if gradient[-1] > 0 else ()))
        ypeak = yval[maxindex]
        maxindex_above_min = maxindex[np.where(ypeak>min_intensity)]
        targetx = xval[maxindex_above_min]
        targety = yval[maxindex_above_min]
        if targetx is None:
            targetx =[]
            targety =[]
        return targetx, targety
    
    def peak_intensity(self, min_intensity = 10000, error = 0.01,peak_error = 0.2):
        """
        the same as self.peaks. return mz list and mz intensity list.
        intensity is calculated if |mzs - targetmz| < error
        
        """
        targetmzs = self.peaks(min_intensity,peak_error)
        if targetmzs == None:
            return [],[]
        num = 0
        targetmzi = []
        for targetmz in targetmzs:
            targetmzi.append(0)
            while num < len(self.mzs):
                if self.mzs[num] - targetmz < -error:
                    num +=1
                elif self.mzs[num] - targetmz < error:
                    targetmzi[-1] += self.mzi[num]
                    num += 1
                else:
                    num = max(num -20, 0)
                    break
        return targetmzs, targetmzi
                

def readMS1File2list(filename):
    """
    given a filename, return a list of element if class MS1fromFMSConvertGUI
    """
    fo = open(filename,"r")
    mylist = fo.readlines()
    fo.close()
    length = len(mylist)
    for num in range(length):
        if mylist[num][0] != "H":
            break
    ms1lis =[]
    scan = None
    mzi =[]
    mzs =[]
    for num2 in range(num,length):
        line = mylist[num2]
        if line[0] == "S":
            if scan != None:
                scan.mzi = mzi
                scan.mzs = mzs
                ms1lis.append(scan)
            scan = MS1fromFMSConvertGUI()
            mzs =[]
            mzi =[]
            scan.scan_num = int(line.split()[2])
        elif "RTime" in line:
            scan.rt = float(line.split()[2])
        elif "BPM" in line:
            scan.bpm = float(line.split()[2])
        elif "BPI" in line:
            scan.bpi = float(line.split()[2])
        elif "TIC" in line:
            scan.tic = float(line.split()[2])
        else:
            mzs.append(float(line.split()[0]))
            mzi.append(float(line.split()[1]))
    if scan != None:
        if scan.scan_num !=None:
            scan.mzi = mzi
            scan.mzs = mzs
            ms1lis.append(scan)
    return ms1lis


def ms12mgflike(filename, outfilename = None,min_intensity = 100000, error = 0.01,peak_error = 0.2,simple = False):
    """
    given a filename, convert it to like mgf file, with peaks and peak intensity
    if outfilename is not provided, save filename+mgflike
    file is ms1 from MSConvertGUI
    output looks like:
    H	CreationDate Thu Nov 12 14:16:49 2015
    H	Extractor	ProteoWizard
    H	Extractor version	Xcalibur
    H	Source file	1380_Control_m100min090215.raw
    S	1	1
    I	RTime	18.0058
    I	BPI	1.1871e+007
    I	BPM	646.2391
    I	TIC	1.097377e+008
    410.1888	1020061
    482.1901	3030660
    487.1455	3765159
    503.1192	1222837
    S	2	2
    ....
    """
    fo = open(filename)
    if outfilename == None:
        outfilename = filename + "mgflike"
    fout = open(outfilename,"w")
    line = fo.readline()
    while line != "":
        if line[0] == "H":
            fout.write(line)
            line = fo.readline()
        else:
            break
    while line != "":
        if line[0] == "S":
            scan = MS1fromFMSConvertGUI()
            scan.scan_num = int(line.split()[2])
            to_write = line
            line = fo.readline()
            scan.rt = float(line.split()[2])
            to_write = to_write + line
            line = fo.readline()
            scan.bpi = float(line.split()[2])
            to_write = to_write + line
            line = fo.readline()
            scan.bpm = float(line.split()[2])
            to_write = to_write + line
            line = fo.readline()
            scan.tic = float(line.split()[2])
            to_write = to_write + line
            mzs = []
            mzi = []
            line = fo.readline()
            if line != "":
                while line[0] != "S":
                    mzs.append(float(line.split()[0]))
                    mzi.append(float(line.split()[1]))
                    line = fo.readline()
                    if line == "":
                        break
            scan.mzs = mzs
            scan.mzi =mzi
            if simple:
                targetmzs, targetmzi = scan.simplePeaks(min_intensity)
            else:
                targetmzs, targetmzi = scan.peak_intensity(min_intensity,error,peak_error)
            if len(targetmzs) !=0:
                fout.write(to_write)
                for num in range(len(targetmzs)):
                    fout.write("%.4f\t%.0f\n"%(targetmzs[num],targetmzi[num]))
#            print(len(targetmzi))
#        break
    fout.close()



def txtThermo2ms1MGFlike(filename, outfilename = None,min_intensity = 100000, error = 0.01,peak_error = 0.2,simple = False):
    """
    Thermo Xcalibur can convert raw file to txt file, with is very big.
    Extract ms1 information from it. the output is similar like ms12mgflike
    """
    fo = open(filename)
    if outfilename == None:
        outfilename = filename + "ms1mgflike"
    fout = open(outfilename,"w")
    line = fo.readline()
    while line != "":
        if line[:10] != "ScanHeader":
            line = fo.readline()
        else:
            break
    while line != "":
        if line[:12] =="ScanHeader #":
            scan = MS1fromFMSConvertGUI()
            scan.scan_num = int(line.split()[-1])
            fo.readline()
            line = fo.readline()
            scan.rt = float(line.split()[2][:-1])
            for _num in range(6):#skip 6 lines
                fo.readline()
            line = fo.readline()
            if "MS2 Scan" not in line:
                for _num in range(5): # skip 5 lines
                    fo.readline()
                mzs =[]
                mzi =[]
                line = fo.readline()
                if line != "":
                    while "ScanHeader #" not in line:
                        if "Packet" in line:
                            mzs.append(float(line.split()[8]))
                            mzi.append(float(line.split()[5][:-1]))
                        line = fo.readline()
                        if line == "":
                            break
                scan.mzi = mzi
                scan.mzs = mzs
                if simple:
                    targetmzs, targetmzi = scan.simplePeaks(min_intensity)
                else:
                    targetmzs, targetmzi = scan.peak_intensity(min_intensity,error,peak_error)
                if len(targetmzs) !=0:
                    fout.write("Scan\t"+str(scan.scan_num)+"\n")
                    fout.write("RTime\t"+str(scan.rt)+"\n")
                    for num in range(len(targetmzs)):
                        fout.write("%.4f\t%.0f\n"%(targetmzs[num],targetmzi[num]))
            else:
                while "ScanHeader #" not in line and line != "":
                    line = fo.readline()
#        break
    fout.close()

                
def ms1ToIndividualScans(filename, outfolder = None):
    """
    given a .ms1 file, break it to individual scan files.
    filename: scan000001.txt.
    each file contsins mz and intensity only.
    also a file with scan number and rt time.
    """
    fo = open(filename)
    if outfolder is None:
        outfolder = ''
    outfilename = 'scan2rt.txt'
    fout = open(outfilename,"w")
    line = fo.readline()
    while line != "":
        if line[0] == "H":
            line = fo.readline()
        else:
            break
    fout2 = None
    to_write = ''
    while line:
        if 'S' in line:
            scan_num = int(line.split()[1])
            line = fo.readline()
            scan_rt = float(line.split()[2])
            fout.write('%d\t%f\n'%(scan_num, scan_rt))
            if fout2 is not None:
                fout2.write(to_write)
                fout2.close()
            fout2 = open(outfolder+'scan%06d.txt'%(scan_num),'w')
            to_write = ''
            fo.readline()
            fo.readline()
            fo.readline()
            line = fo.readline()
            while line:
                if 'S' not in line:
                    to_write = to_write + line
                    line = fo.readline()
                else:
                    break
    fout2.write(to_write)
    fout2.close()
    fout.close()
            
def ms1CleanUp(filename, mzi_min = 10000, outname = None):
    """
    for ms1 file from MSConvertGUI, usually the values with mzi <10000 are noises, 
    delete these lines to make the file a little smaller.
    also, based on the ms1 file from Orbitrap or Fusion machine,
    mz's resolution is about 0.001.
    outname will be filename + '.clear' in default setting.
    """
    fo = open(filename,'r')
    if outname is None:
        outname = filename + '.clear'
    fout = open(outname,'w')
    numbers = '0123456789'
    for ele in fo:
        if ele[0] in numbers:
            mz, mzi = ele.split()
            if float(mzi) > mzi_min:
                fout.write(ele)
        else:
            fout.write(ele)
    fo.close()
    fout.close()

def readms1_iter(filename,mzi_min = 10000):
    """
    given a ms1 file or ms1.clear, return generator, yield scan_num, rt, lsmzs, lsmzi
    mzi_min is the minimum value for mzi to include in the result
    """
    fo = open(filename)
    line = fo.readline()
    while line != "":
        if line[0] == "H":
            line = fo.readline()
        else:
            break
    lsmzs = {}
    lsmzi = {}
    scan_num = None
    rt = None
    scan_info = []
    while line:
        if 'S' in line:
            if scan_num is not None:
                lsmzs = np.array(lsmzs)
                lsmzi = np.array(lsmzi)
                yield scan_num, rt, lsmzs, lsmzi
            scan_num = int(line.split()[1])
            line = fo.readline()
            rt = float(line.split()[-1])
            scan_info.append((scan_num, rt))
            fo.readline()
            fo.readline()
            fo.readline()
            lsmzs = []
            lsmzi = []
            line = fo.readline()
            while line:
                if 'S' not in line:
                    if len(line.split()) == 2:
                        mzv, mzi = line.split()
                    else:
                        print(line)
                        print('something wrong with this line')
                        line = ''
                        break
                    mzv = float(mzv)
                    mzi = float(mzi)
                    if mzi >= mzi_min:
                        lsmzs.append(mzv)
                        lsmzi.append(mzi)
                    line = fo.readline()
                else:
                    break
    fo.close()
    lsmzs = np.array(lsmzs)
    lsmzi = np.array(lsmzi)
    yield scan_num, rt, lsmzs, lsmzi


def readms1ToMzsMzi(filename,mzi_min = 10000, outputScanNum = False):
    """
    given a ms1 file or ms1.clear, return 2 dictionary, 1 for mzs, 1 for mzi
    mzs and mzi are stored in np.array format to save memory
    if outputScanNum is True, return a list, with scan_num and rt
    """
    ls = list(readms1_iter(filename, mzi_min))
    scan_info = [ele[:2] for ele in ls]
    npmzs = np.array([ele[2] for ele in ls])
    npmzi = np.array([ele[3] for ele in ls])
    if outputScanNum:
        return npmzs, npmzi, scan_info
    else:
        return npmzs, npmzi


def ms1ToPeaks(filename, outfilename = None,mzi_min = 1000, charge_max = 10, peakfoldmin = 5 ,peakonly = False):
    """
    given a ms1, convert it mgf like file, with peaks and peak intensity
    charge_max = 10, which means peaks whose m/z difference is less than 1/11, only the major peak will be kept.
    only peaks with mzi greater than mzi_min will be kept.
    if outfilename is None, outfilename = filename + '.Peak'
    output file looks like:
    '
    Scan: 1
    RTime: 14
    mzs: 360 362 399
    mzis: 1000 500 1000
    '
    if peakonly is True, output peaks directly with no quality control with charge_max and peakfoldmin.
    charge_max and peakfoldmin can help to remove small peaks around a big peak.
    """
    f = readms1_iter(filename, 0)
    if outfilename is None:
        outfilename = filename + '.Peak'
    fout = open(outfilename, 'w')
    diffmin = 1/(charge_max+1)        
    for scan_num, rt, lsmzs, lsmzi in f:
#        print(lsmzs)
        fout.write('Scan: %d\n'%scan_num)
        fout.write('RTime: %f\n'%rt)
        if len(lsmzi) < 3:
            fout.write('mzs: \nmzis: \n')
        else:
            peaks = []
            for num in range(1, len(lsmzi) -1):
                if lsmzi[num] > 0:
                    if lsmzi[num] > mzi_min:
                        if lsmzi[num] > lsmzi[num -1] and lsmzi[num] > lsmzi[num +1]:
                            if lsmzi[num -1] >0 and lsmzi[num +1] >0:
                                peaks.append((lsmzs[num], lsmzi[num]))
#            print('peaks top 10 ',peaks[:10])
            if peakonly:
                mzs_str = ' '.join([str(ele[0]) for ele in peaks])
                mzi_str = ' '.join([str(ele[1]) for ele in peaks])
                fout.write("mzs: " + mzs_str+'\n')
                fout.write("mzi: " + mzi_str+'\n')
            elif len(peaks) == 0:
                fout.write('mzs: \nmzis: \n')
            elif (len(peaks)) == 1:
                fout.write('mzs: %f\nmzis: %f\n'%(peaks[0][0],peaks[0][1]))
            else:
                peaks_keep = []
                mzs = np.array([ele[0] for ele in peaks])
                mzis = np.array([ele[1] for ele in peaks])
                while len(mzs) > 0:
                    to_removeID = []
                    max_id = mzis.argmax()
                    to_removeID.append(max_id)
                    peaks_keep.append((mzs[max_id], mzis[max_id]))
                    if len(mzs) > 0:
                        for num in range(max_id +1, len(mzs)):
                            if abs(mzs[num] - peaks_keep[-1][0]) < diffmin:
                                if peaks_keep[-1][1] / mzis[num] > peakfoldmin:
                                    to_removeID.append(num)
                            else:
                                break
                        for num in range(max_id -1, -1,-1):
                            if abs(mzs[num] - peaks_keep[-1][0]) < diffmin:
                                if peaks_keep[-1][1] / mzis[num] > peakfoldmin:
                                    to_removeID.append(num)
                            else:
                                break
                        mzs = np.delete(mzs,to_removeID)
                        mzis = np.delete(mzis,to_removeID)
                peaks_keep.sort(key = lambda x:x[0])
                mzs_str = ' '.join([str(ele[0]) for ele in peaks_keep])
                mzi_str = ' '.join([str(ele[1]) for ele in peaks_keep])
                fout.write("mzs: " + mzs_str+'\n')
                fout.write("mzi: " + mzi_str+'\n')
#        break
    fout.close()
                    
                        
def readms1Peaks2mzsmzi(filename,mzi_min = 0):
    '''
    given a file from ms1ToPeaks, return a numpy array of mzs, mzi and scan_info.
    scan_info have two field. first is scan num, second RTtime.
    '''
    fo = open(filename,'r')
    ls = fo.readlines()
    scan_info = []
    npmzs = []
    npmzi = []
    for num in range(0, len(ls), 4):
        scan_num = int(ls[num].split()[1])
        rt = float(ls[num+1].split()[1])
        mzs = [float(i) for i in ls[num+2].split()[1:]]
        mzi = [float(i) for i in ls[num+3].split()[1:]]
        scan_info.append((scan_num,rt))
        npmzs.append(np.array(mzs))
        npmzi.append(np.array(mzi))
    scan_info = np.array(scan_info)
    npmzs = np.array(npmzs)
    npmzi = np.array(npmzi)
    return npmzs, npmzi, scan_info
    

        
    

def clean_up_targetmzis(targetmzis, clean_up = 3):
    """
    for a list of targetmzis, if adjacent 3 are more values are large than 0, keep them. else, set them to 0
    """
    if clean_up == 0 or clean_up == 1:
        return targetmzis

    l = len(targetmzis)
    if l < clean_up:
        return [0 for i in range(l)]
    
    lsls =[]
    sublis = []
    import itertools
    sublis.append(targetmzis[0])
    for num in range(1,len(targetmzis)):
        mzi = targetmzis[num]
        if sublis[-1] > 0 and mzi >0 or sublis[-1] <=0 and mzi <=0:
            sublis.append(mzi)
        else:
            lsls.append(sublis)
            sublis =[]
            sublis.append(mzi)
    lsls.append(sublis)
    for num in range(len(lsls)):
        if len(lsls[num]) < clean_up:
            lsls[num] = [0 for ele in lsls[num]]
    return list(itertools.chain.from_iterable(lsls))

def targetmzIntensityFromNpmzsmzi(targetmz, scans_mzs, scans_mzi, error = 0.003,clean_up = 3):
    """
    given a target_mz, return a list of it's intensity in each scan.
    scans_mzs and scan_mzi are from the function readms1ToMzsMzi, npmzs and npmzi.
    if no mz in scan_mzs is very close to targetmz (with error of 0.003 in default condition), intensity is 0
    clean_up = 3 means at least 3 adjacent values should >0 to report, else the value will be set to 0. 
    clean_up 3 can help to reduce noises.
    """
    import numpy as np
    targetmzis = []
    for num in range(len(scans_mzs)):
        scan_mz = scans_mzs[num]
        scan_mi = scans_mzi[num]
        abs_array = np.abs(np.array(scan_mz) - targetmz)
        if len(abs_array) > 0:
            index_bestmz = abs_array.argmin()
            bestmz = scan_mz[index_bestmz]
            if abs(bestmz - targetmz) > error:
                targetmzis.append(0)
            else:
                targetmzis.append(scan_mi[index_bestmz])
        else:
            targetmzis.append(0)
    return clean_up_targetmzis(targetmzis, clean_up)
            
def targetPepmzIntensityFromNpmzsmzi(peptide, scans_mzs, scans_mzi, charge = 2, min_extra_neutron = 0, error = 0.003, clean_up = 3):
    """
    given a peptide sequence, and its charge, return a list of it's intensity in each scan. 
    Function pretty similar as targetmzIntensityFromNpmzsmzi
    """
    import chemicalAllPossibleMsms
    targetmzs = []
    ls_targetmzis = []
    for num in range(min_extra_neutron+1):
        targetmz = chemicalAllPossibleMsms.getPepmassWithModi(peptide, charge, num)
        targetmzs.append(targetmz)
#        print(targetmz)
        targetmzis = targetmzIntensityFromNpmzsmzi(targetmz, scans_mzs, scans_mzi, error, 0)
        ls_targetmzis.append(targetmzis)
    mzis = []
    n = len(ls_targetmzis)
    for num in range(len(targetmzis)):
        mzs = [ls_targetmzis[i][num] for i in range(n)]
        if min(mzs) != 0:
            mzis.append(ls_targetmzis[0][num])
        else:
            mzis.append(0)
    return clean_up_targetmzis(mzis, clean_up)

    

    