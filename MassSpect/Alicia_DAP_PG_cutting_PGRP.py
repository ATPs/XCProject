from ms1Class import *

folder = "E:\\Lab\\works\\20160212AliciaMSMS\\RAW\\"
files = ["1-1402_PGRP1_40hr.ms1","2-1402_DAP_PG.ms1","3-1402_LYS_PG.ms1","4-1402_PGRP1_DAP_PG_40hr.ms1",\
"5-1402_PGRP1_LYS_PG_40hr.ms1","6-1402_PGRP1_DAP_PG_70hr.ms1","7-1402_PGRP1_LYS_PG_70hr.ms1"]
filenamedic ={"1":files[0],"D":files[1],"L":files[2],"1D4":files[3],"1L4":files[4],"1D7":files[5],"1L7":files[6]}
import numpy
import peakutils
from peakutils.plot import plot as pplot
from matplotlib import pyplot as plt


testms1 = readMS1File2list(folder + files[0]+".filtered")

mzs,mzi = testms1[0].peak_intensity(1e4,0.01)
print(len(mzs))
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1,aspect='equal')
plt.axis([-10,max(mzs),-10,10])
import numpy as np
for num in range(len(mzs)):
    ax.add_patch(plt.Circle(xy=(mzs[num],0),radius=np.sqrt(mzi[num])/100,color = "g",alpha = 0.3,fill = False,lw=0.05))
#plt.scatter(mzs,mzi,marker="o", s=300,zorder=10)
plt.savefig("test.pdf")
plt.close()

"""
plot testms1
"""
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1,aspect='equal')
plt.axis([100,900,800,7000])
plt.rcParams =plt.rcParamsDefault
plt.rcParams["figure.figsize"] = [20,60]
import numpy as np
for num2 in range(len(testms1)):
    mzs,mzi = testms1[num2].peak_intensity(1e5,0.01)
    if mzs !=None:
        for num in range(len(mzs)):
            ax.add_patch(plt.Circle(xy=(mzs[num],900+10*(testms1[num2].rt*60-900)),radius=np.sqrt(mzi[num])/400,color = "g",alpha = 0.3,fill = False,lw=0.2))
    #plt.scatter(mzs,mzi,marker="o", s=300,zorder=10)
plt.savefig("test.pdf",bbox_inches='tight')
plt.close()


"""
plot DAP-PG
"""
ms1dic ={}
for key in ["1","D","1D4","1D7"]:
     ms1dic[key] = readMS1File2list(folder + filenamedic[key]+".filtered")
colordic ={"1":"b","D":"b","1D4":"r","1D7":"r"}

fig = plt.figure(1)
ax = fig.add_subplot(1,1,1,aspect='equal')
plt.axis([100,900,9000,15000])
plt.rcParams =plt.rcParamsDefault
plt.rcParams["figure.figsize"] = [20,60]
import numpy as np
for key in colordic:
    for num2 in range(len(ms1dic[key])):
        mzs,mzi = ms1dic[key][num2].peak_intensity(1e5,0.01)
        if mzs !=None:
            for num in range(len(mzs)):
                ax.add_patch(plt.Circle(xy=(mzs[num],10*(ms1dic["1"][num2].rt*60)),radius=np.sqrt(mzi[num])/1000,color = colordic[key],alpha = 0.4,fill = False,lw=0.2))
    #plt.scatter(mzs,mzi,marker="o", s=300,zorder=10)
plt.savefig("20160229DAP-PG3.pdf",bbox_inches='tight')
plt.close()


allmzs =[]
for key in colordic:
    for num2 in range(len(ms1dic[key])):
        mzs,mzi = ms1dic[key][num2].peak_intensity(1e5,0.01)
        if mzs !=None:
            allmzs = allmzs + mzs
print(len(allmzs))
allmzs = set(allmzs)
print(len(allmzs))
allmzs2 =set()
for ele in allmzs:
    allmzs2.add(round(ele,2))
print(len(allmzs2))

mzdic ={}
for ele in allmzs2:
    ele = round(ele,2)
    mzdic[ele] = [0,0,0,0]
lis = ["1","D","1D4","1D7"]
for key in range(4):
    for ele in ms1dic[lis[key]]:
        mzs,mzi = ele.peak_intensity(1e5,0.005)
        for num3 in range(len(mzs)):
            mz = round(mzs[num3],2)
            if mz in mzdic:
                mzdic[mz][key] += mzi[num3]
allmzs = list(allmzs2)

selected ={}
for ele in mzdic:
    mzis = mzdic[ele]
    if mzis[2] / max((mzis[1]+1),(mzis[0]+1)) >4 or mzis[3] / max((mzis[1]+1),(mzis[0]+1)) >4:
        selected[ele] = mzdic[ele]
print(len(selected))

fo = open("list.txt","w")
for key in selected:
    fo.write(str(round(key,2)) +"\t"+str(selected[key])+"\n")
fo.close()