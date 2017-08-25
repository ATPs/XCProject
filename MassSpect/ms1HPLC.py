# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:12:33 2016

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


folder = "E:\\Lab\\works\\20160212AliciaMSMS\\RAW\\"
files = ["1-1402_PGRP1_40hr.ms1","2-1402_DAP_PG.ms1","3-1402_LYS_PG.ms1","4-1402_PGRP1_DAP_PG_40hr.ms1",\
"5-1402_PGRP1_LYS_PG_40hr.ms1","6-1402_PGRP1_DAP_PG_70hr.ms1","7-1402_PGRP1_LYS_PG_70hr.ms1"]
filenamedic ={"1":files[0],"D":files[1],"L":files[2],"1D4":files[3],"1L4":files[4],"1D7":files[5],"1L7":files[6]}


dap_products=["C32H55N9O15","C29H49N7O15","C3H8N2O","C26H44N6O14","C6H13N3O2","C19H32N4O11",\
"C13H25N5O5","C14H24N2O9","C18H33N7O7","C11H19N1O8","C21H38N8O8"]
lys_products=["C31H55N9O13","C28H49N7O13","C3H8N2O","C25H44N6O12","C6H13N3O2","C19H32N4O11",\
"C12H25N5O3","C14H24N2O9","C17H33N7O5","C11H19N1O8","C20H38N8O6"]
# a list of all possible molecules from MPP-DAP, the first one is MPP-DAP itself
def listmoleculemz(moleculelist,ioncharge = 1):
    """
    given a list of molecules, return list of m/z's.
    """
    from pyteomics import mass
    mzlist =[]
    for ele in moleculelist:
        mz = mass.calculate_mass(formula = ele, charge = ioncharge)
        mzlist.append(mz)
    return mzlist
dap_mzs = listmoleculemz(dap_products)
lys_mzs = listmoleculemz(lys_products)


#plot dappg whole molecule
dappg_lcdic={}
for key in ["1","D","1D4","1D7"]:
    dappg_lcdic[key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key],dap_mzs[0])
##plot in 4 subplots
import matplotlib.pyplot as plt
fig = plt.figure()
plt.rcParams["figure.figsize"] = [12.0,10]
lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
daplis = ["1","D","1D4","1D7"]
for num in range(4):
    plt.subplot(4,1,num+1)
    plt.plot(dappg_lcdic[daplis[num]][0],dappg_lcdic[daplis[num]][1],zorder = num+3,color =lisColor[num],label = daplis[num])
    plt.xlim(xmin=10,xmax=50)
    plt.legend()
fig.tight_layout()
plt.savefig("20160224DAP-PG.pdf")
plt.savefig("20160224DAP-PG.png",dpi=300)
plt.close()
##plot in one figure
import matplotlib.pyplot as plt
fig = plt.figure()
plt.rcParams["figure.figsize"] = [12.0,10]
lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
daplis = ["1","D","1D4","1D7"]
for num in range(4):
    plt.plot(dappg_lcdic[daplis[num]][0],dappg_lcdic[daplis[num]][1],"--",dashes =[4,2],zorder = num+3,color =lisColor[num],label = daplis[num])
    plt.xlim(xmin=5,xmax=45)
fig.tight_layout()
plt.savefig("20160224DAP-PG_all.pdf")
plt.savefig("20160224DAP-PG_all.png",dpi=300)
plt.close()

#plot all dappg
dappg_lclis=[{} for num in range(11)]
for num in range(11):
    for key in ["1","D","1D4","1D7"]:
        dappg_lclis[num][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key],dap_mzs[num])

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for num1 in range(11):
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    daplis = ["1","D","1D4","1D7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(dappg_lclis[num1][daplis[num]][0],dappg_lclis[num1][daplis[num]][1],zorder = num+3,color =lisColor[num],label = daplis[num])
        plt.xlim(xmin=10,xmax=50)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160224DAP-PG"+dap_products[num1]+".pdf")
    plt.savefig("20160224DAP-PG"+dap_products[num1]+".png",dpi=300)
    plt.close()


#plot all lyspg
lyspg_lclis=[{} for num in range(11)]
for num in range(11):
    for key in ["1","L","1L4","1L7"]:
        lyspg_lclis[num][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key],lys_mzs[num])

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for num1 in range(11):
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    lyslis = ["1","L","1L4","1L7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(lyspg_lclis[num1][lyslis[num]][0],lyspg_lclis[num1][lyslis[num]][1],zorder = num+3,color =lisColor[num],label = lyslis[num])
        plt.xlim(xmin=10,xmax=50)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160224LYS-PG"+lys_products[num1]+".pdf")
    plt.savefig("20160224LYS-PG"+lys_products[num1]+".png",dpi=300)
    plt.close()


for num1 in [6]:
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    lyslis = ["1","L","1L4","1L7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(lyspg_lclis[num1][lyslis[num]][0],lyspg_lclis[num1][lyslis[num]][1],zorder = num+3,color =lisColor[num],label = lyslis[num])
        plt.xlim(xmin=26,xmax=27)
        plt.ylim(ymax = 1e6)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160224LYS-PG"+lys_products[num1]+"2530.pdf")
    plt.savefig("20160224LYS-PG"+lys_products[num1]+"2530.png",dpi=300)
    plt.close()

#clean the ms1 files. RT from 10min to 30min. also remove those pairs with intensity ==0
#20160226
folder = "E:\\Lab\\works\\20160212AliciaMSMS\\RAW\\"
files = ["1-1402_PGRP1_40hr.ms1","2-1402_DAP_PG.ms1","3-1402_LYS_PG.ms1","4-1402_PGRP1_DAP_PG_40hr.ms1",\
"5-1402_PGRP1_LYS_PG_40hr.ms1","6-1402_PGRP1_DAP_PG_70hr.ms1","7-1402_PGRP1_LYS_PG_70hr.ms1"]
filenamedic ={"1":files[0],"D":files[1],"L":files[2],"1D4":files[3],"1L4":files[4],"1D7":files[5],"1L7":files[6]}

def ms1FileClean(filename,timemin = 600,timemax=1800, minIntensity = 0,outfile = None):
    """
    given a filename of ms1 data, output a ms1 file, based on those filters
    if outfile name is not provided, outfile = filename+ "filtered"
    timemin and timemax is in seconds of retension time
    """
    fo = open(filename,"r")
    if outfile == None:
        fout = open(filename+".filtered","w")
    else:
        fout = open(outfile,"w")
    outlist=[]
    scan =[]
    for ele in fo:
        if ele[0] == "H":
            fout.write(ele)
        elif ele[0] == "S":
            if scan != []:
                rt =float(scan[1].split()[2])
                if rt*60 > timemin and rt*60<timemax:
                    outlist.append(scan)
                if rt*60 > timemax:
                    break
            scan=[]
            scan.append(ele)
        elif ele[0] == "I":
            scan.append(ele)
        else:
            if float(ele.split()[1])>minIntensity:
                scan.append(ele)
    #the last scan is not processed in the for loop
    if scan != []:
        rt =float(scan[1].split()[2])
        rtlast = float(outlist[-1][1].split()[2])
        if rt != rtlast:
            if rt*60 > timemin and rt*60<timemax:
                outlist.append(scan)
    for scan in outlist:
        for line in scan:
            fout.write(line)
    print("scans in output is ",len(outlist))
    fout.close()

#clean all 7 files
for file in files:
    ms1FileClean(folder+file,900,1500,-1)

#clean 3 PGRP2 files
files = ["8-1402_PGRP2_zinc_40hr.ms1","9-1402_PRGP2_MPP_DAP_zinc_40hrs.ms1","10-1402_PGRP2_MPP_LYS_zinc_40hr.ms1"]
for file in files:
    ms1FileClean(folder+file,900,1500,-1)



#20160229
#products including dehydrate, methylated and +2 ions
dap_products=["C32H55N9O15","C29H49N7O15","C3H8N2O","C26H44N6O14","C6H13N3O2","C19H32N4O11",\
"C13H25N5O5","C14H24N2O9","C18H33N7O7","C11H19N1O8","C21H38N8O8"]
# a list of all possible molecules from MPP-DAP, the first one is MPP-DAP itself
dap_mzs = listmoleculemz(dap_products)
lys_mzs = listmoleculemz(lys_products)

from pyteomics import mass
#plot all dappg
#plot all dappg -H2o
dappg_lclis=[{} for num in range(11)]
for num in range(11):
    for key in ["1","D","1D4","1D7"]:
        dappg_lclis[num][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key]+".filtered",\
        dap_mzs[num] -mass.calculate_mass(formula ="H2O",charge = 0))

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for num1 in range(11):
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    daplis = ["1","D","1D4","1D7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(dappg_lclis[num1][daplis[num]][0],dappg_lclis[num1][daplis[num]][1],zorder = num+3,\
        color =lisColor[num],label = daplis[num]+" "+dap_products[num1]+"_-H2O")
        plt.xlim(xmin=14,xmax=26)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_-H2O.pdf")
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_-H2O.png",dpi=300)
    plt.close()


#plot all dappg +CH2
dappg_lclis=[{} for num in range(11)]
for num in range(11):
    for key in ["1","D","1D4","1D7"]:
        dappg_lclis[num][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key]+".filtered",\
        dap_mzs[num] + mass.calculate_mass(formula ="CH2",charge = 0))

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for num1 in range(11):
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    daplis = ["1","D","1D4","1D7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(dappg_lclis[num1][daplis[num]][0],dappg_lclis[num1][daplis[num]][1],zorder = num+3,\
        color =lisColor[num],label = daplis[num]+" "+dap_products[num1]+"_CH2")
        plt.xlim(xmin=14,xmax=26)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_CH2.pdf")
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_CH2.png",dpi=300)
    plt.close()

#plot all dappg -H2o +2
dappg_lclis=[{} for num in range(11)]
for num in range(11):
    for key in ["1","D","1D4","1D7"]:
        dappg_lclis[num][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key]+".filtered",\
        mass.calculate_mass(dap_products[num],charge = 2) -mass.calculate_mass(formula ="H2O",charge = 0)/2)

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for num1 in range(11):
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    daplis = ["1","D","1D4","1D7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(dappg_lclis[num1][daplis[num]][0],dappg_lclis[num1][daplis[num]][1],zorder = num+3,\
        color =lisColor[num],label = daplis[num]+" "+dap_products[num1]+"_-H2Ocharge2")
        plt.xlim(xmin=14,xmax=26)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_-H2Ocharge2.pdf")
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_-H2Ocharge2.png",dpi=300)
    plt.close()


#plot all dappg +CH2 +2
dappg_lclis=[{} for num in range(11)]
for num in range(11):
    for key in ["1","D","1D4","1D7"]:
        dappg_lclis[num][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key]+".filtered",\
        mass.calculate_mass(dap_products[num],charge = 2) +mass.calculate_mass(formula ="CH2",charge = 0)/2)

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for num1 in range(11):
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    daplis = ["1","D","1D4","1D7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(dappg_lclis[num1][daplis[num]][0],dappg_lclis[num1][daplis[num]][1],zorder = num+3,\
        color =lisColor[num],label = daplis[num]+" "+dap_products[num1]+"_CH2charge2")
        plt.xlim(xmin=14,xmax=26)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_CH2charge2.pdf")
    plt.savefig("20160229DAP-PG"+dap_products[num1]+"_CH2charge2.png",dpi=300)
    plt.close()




#plot all lyspg
lys_products=["C31H55N9O13","C28H49N7O13","C3H8N2O","C25H44N6O12","C6H13N3O2","C19H32N4O11",\
"C12H25N5O3","C14H24N2O9","C17H33N7O5","C11H19N1O8","C20H38N8O6"]
lysmzdic ={}
from pyteomics import mass
for ele in lys_products:
    lysmzdic[ele +"-H2O"] = mass.calculate_mass(formula = ele, charge = 1) -mass.calculate_mass(formula ="H2O",charge = 0)
    lysmzdic[ele +"+CH2"] = mass.calculate_mass(formula = ele, charge = 1) + mass.calculate_mass(formula ="CH2",charge = 0)
    lysmzdic[ele +"-H2O charge2"] = mass.calculate_mass(formula = ele, charge = 2) -mass.calculate_mass(formula ="H2O",charge = 0)/2
    lysmzdic[ele +"+CH2 charge2"] = mass.calculate_mass(formula = ele, charge = 2) + mass.calculate_mass(formula ="CH2",charge = 0)/2


lyspg_lcdic={}
for key2 in lysmzdic:
    lyspg_lcdic[key2] = {}
    for key in ["1","L","1L4","1L7"]:
        lyspg_lcdic[key2][key] = extractMZms1FromFMSConvertGUI(folder+filenamedic[key]+".filtered",lysmzdic[key2])

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,10]

for key in lyspg_lcdic:
    fig = plt.figure()
    lisColor = [(0.2,0.2,0.2,0.5),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7)]
    lyslis = ["1","L","1L4","1L7"]
    for num in range(4):
        plt.subplot(4,1,num+1)
        plt.plot(lyspg_lcdic[key][lyslis[num]][0],lyspg_lcdic[key][lyslis[num]][1],zorder = num+6,color =lisColor[num],label = lyslis[num] + key)
        plt.xlim(xmin=14,xmax=26)
        plt.legend()
    fig.tight_layout()
    plt.savefig("20160229LYS-PG"+key+".pdf")
    plt.savefig("20160229LYS-PG"+key+".png",dpi=300)
    plt.close()

#plot all PG for PGRP2
lys_products=["C31H55N9O13","C28H49N7O13","C3H8N2O","C25H44N6O12","C6H13N3O2","C19H32N4O11",\
"C12H25N5O3","C14H24N2O9","C17H33N7O5","C11H19N1O8","C20H38N8O6"]
dap_products=["C32H55N9O15","C29H49N7O15","C3H8N2O","C26H44N6O14","C6H13N3O2","C19H32N4O11",\
"C13H25N5O5","C14H24N2O9","C18H33N7O7","C11H19N1O8","C21H38N8O8"]
pgmzdic = {}
from pyteomics import mass
for ele in lys_products:
    pgmzdic["LYSPG " + ele +" charge2"] = mass.calculate_mass(formula = ele, charge = 2)
    pgmzdic["LYSPG " + ele] = mass.calculate_mass(formula = ele, charge = 1)
    pgmzdic["LYSPG " + ele +"-H2O"] = mass.calculate_mass(formula = ele, charge = 1) -mass.calculate_mass(formula ="H2O",charge = 0)
    pgmzdic["LYSPG " + ele +"+CH2"] = mass.calculate_mass(formula = ele, charge = 1) + mass.calculate_mass(formula ="CH2",charge = 0)
    pgmzdic["LYSPG " + ele +"-H2O charge2"] = mass.calculate_mass(formula = ele, charge = 2) -mass.calculate_mass(formula ="H2O",charge = 0)/2
    pgmzdic["LYSPG " + ele +"+CH2 charge2"] = mass.calculate_mass(formula = ele, charge = 2) + mass.calculate_mass(formula ="CH2",charge = 0)/2
for ele in dap_products:
    pgmzdic["DAPPG " + ele +" charge2"] = mass.calculate_mass(formula = ele, charge = 2)
    pgmzdic["DAPPG " + ele] = mass.calculate_mass(formula = ele, charge = 1)
    pgmzdic["DAPPG " + ele +"-H2O"] = mass.calculate_mass(formula = ele, charge = 1) -mass.calculate_mass(formula ="H2O",charge = 0)
    pgmzdic["DAPPG " + ele +"+CH2"] = mass.calculate_mass(formula = ele, charge = 1) + mass.calculate_mass(formula ="CH2",charge = 0)
    pgmzdic["DAPPG " + ele +"-H2O charge2"] = mass.calculate_mass(formula = ele, charge = 2) -mass.calculate_mass(formula ="H2O",charge = 0)/2
    pgmzdic["DAPPG " + ele +"+CH2 charge2"] = mass.calculate_mass(formula = ele, charge = 2) + mass.calculate_mass(formula ="CH2",charge = 0)/2

folder = "E:\\Lab\\works\\20160212AliciaMSMS\\RAW\\"
files = ["1-1402_PGRP1_40hr.ms1","2-1402_DAP_PG.ms1","3-1402_LYS_PG.ms1","4-1402_PGRP1_DAP_PG_40hr.ms1",\
"5-1402_PGRP1_LYS_PG_40hr.ms1","6-1402_PGRP1_DAP_PG_70hr.ms1","7-1402_PGRP1_LYS_PG_70hr.ms1",\
"8-1402_PGRP2_zinc_40hr.ms1","9-1402_PRGP2_MPP_DAP_zinc_40hrs.ms1","10-1402_PGRP2_MPP_LYS_zinc_40hr.ms1"]
filenamedic ={"1":files[0],"D":files[1],"L":files[2],"1D4":files[3],"1L4":files[4],"1D7":files[5],"1L7":files[6],"2":files[7],"2D":files[8],"2L":files[9]}

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [12.0,12]

for key in pgmzdic:
    key2s = ["2","D","2D","L","2L"]
    fig = plt.figure()
    lisColor = ["r","g","b","c","y","m"]
    for num in range(len(key2s)):
        key2 = key2s[num]
        mzs,mzi = extractMZms1FromFMSConvertGUI(folder+filenamedic[key2]+".filtered",pgmzdic[key])
        plt.subplot(len(key2s),1,num+1)
        plt.plot(mzs,mzi,color = lisColor[num], zorder = num+6 )
        plt.xlabel(key2s[num]+" "+key+" mz " + str(pgmzdic[key]))
        plt.xlim(xmin=14,xmax=26)
    fig.tight_layout()
    plt.savefig("20160301"+key+".pdf")
    plt.savefig("20160301"+key+".png",dpi=300)
    plt.close()
#    break
        