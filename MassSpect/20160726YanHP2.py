# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:19:38 2016

@author: k
"""
def extractMZms1FromFMSConvertGUI(filename,targetmz,error = 0.01):
    """
    given a filename of .ms1 converted by MSConvertGUI
    return a list of rt and intensity
    intensity is the sum up of intensities of m/z's with differences compared to targetmz smaller than error.
    for example, below is the content in the file.
    extractMZms1FromFMSConvertGUI(filename, 130.391) 
    142.0917+299.7396+375.319+337.8319+206.0124 = 
    should return lcRT =[64.99898] lcIntensity = [1360.9946]
    I	RTime	64.99898
    I	BPI	212604.2
    I	BPM	149.0234
    I	TIC	2118927
    130.3899 0
    130.3902 0
    130.3904 142.0917
    130.3907 299.7396
    130.391 375.319
    130.3913 337.8319
    130.3916 206.0124
    130.3919 0
    """
    fo = open(filename,"r")
    line = fo.readline()
    intensity = 0
    rt = 0
    lcRT=[]
    lcIntensity =[]
    while line:
        if "RTime" in line:
            if intensity > 0:
                lcRT.append(rt)
                lcIntensity.append(intensity)
            rt = float(line.split("\t")[2][:-1])
            intensity = 0
        if "\t" not in line:
            mz, mzi = line.split()
            mz = float(mz)
            mzi = float(mzi)
            if abs(mz - targetmz) < 0.01:
                intensity += mzi
        line = fo.readline()
    if intensity > 0:
        lcRT.append(rt)
        lcIntensity.append(intensity)
    return lcRT, lcIntensity


def extractMZms1FromFMSConvertGUISave(filename,targetmz,error = 0.01, outname='out.txt'):
    """
    save result to outname
    """
    lcRT, lcIntensity = extractMZms1FromFMSConvertGUI(filename, targetmz, error)
    fo = open(outname,'w')
    for num in range(len(lcRT)):
        fo.write("%f\t%f\n"%(lcRT[num],lcIntensity[num]))
    fo.close()
        



folder = 'E:\\Lab\\works\\20160614YangMSMS\\new\\'
f_control = folder + '1436_control_norm.ms1'
f_fulllength = folder + '1436_proHP2_norm.ms1'
f_cleaved = folder + '1436_cleaved_HP2_norm.ms1'


target = 436.75106
error = 0.005
extractMZms1FromFMSConvertGUISave(f_control, target, error, folder + 'control.txt')
extractMZms1FromFMSConvertGUISave(f_fulllength,target, error, folder + 'full_length.txt')
extractMZms1FromFMSConvertGUISave(f_cleaved,target, error, folder + 'cleaved.txt')
