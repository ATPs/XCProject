#20151214#

"""
for signal-3L result processing
"""

my3L = open("D:\\mine\\OneDrive\\Lab\\Jiang Lab Pubilic\\ProteomeAg\\Signal3L result.txt").read().split("Hong-Bin")

fo = open("D:\\mine\\OneDrive\\Lab\\Jiang Lab Pubilic\\ProteomeAg\\ProteinSeq_all.txt","r")

from Bio import SeqIO
myFa = list(SeqIO.parse(fo,"fasta"))

my3Lk =[]
for ele in my3L:
    if "According to Signal-3L engine for your selected species, your input sequence does not include a signal peptide." not in ele:
       my3Lk.append(ele)

myK = []
for ele2 in myFa:
    for ele in my3Lk:
        if str(ele2.seq[:30]) in ele:
            myK.append((ele2.id,ele))

fo = open("signal3L-out_with.txt","w")
for ele in myK:
    fo.write(ele[0]+"\t"+ele[1])
fo.close()





"""
plot peptides on protein
"""
import matplotlib.pyplot as plt
import numpy as np

NUM_COLORS = 12

cm = plt.get_cmap('Accent')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
for i in range(NUM_COLORS):
    ax.plot(np.arange(10)*(i+1))

fig.savefig('moreColors.pdf')
plt.show()






import matplotlib.pyplot as plt
import numpy as np

folder = "E:\\store_for_D\\XuesongMassRawData\\txt\\band\\"
templist1 = open(folder+"test.txt").readlines()
mylist =[]
for ele in templist1[1:]:
    mylist.append(ele.split("\t"))
mydic={}
for ele in mylist:
    abundance = int(ele[3])
    if abundance not in mydic:
        mydic[abundance] =[]
    mydic[abundance].append(ele)



from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('AgProteinSlides5o.pdf')
NUM_COLORS = 12
cm = plt.get_cmap('nipy_spectral')
TOP = 300
totalcontrol = 1
totalinduced = 1
minamount = 1e7
for proteins in mydic:
    if proteins <= TOP:
        plt.title(str(proteins)+" "+ mydic[proteins][0][2]+" "+mydic[proteins][0][31],fontsize =6)
        for peptides in mydic[proteins]:
            plt.plot([0,int(peptides[4])],[0,0],'k',linewidth=0.1)
            baseamount = 0
            baseamount2 = 0
            totalcontrol = 0
            totalinduced = 0
            for slides in range(12):
                totalcontrol += int(peptides[slides+7])
                totalinduced += int(peptides[slides+7+12])
            aglengend1 =[]
            aglengend2 =[]
            if totalcontrol+totalinduced >= minamount:
                totalcontrol = 1
                totalinduced = 1
                for slides in range(12):
                    if totalcontrol >0:
                        p1 = plt.bar(left = int(peptides[5]), width = int(peptides[6]),\
                        height = int(peptides[slides+7])/totalcontrol, bottom = baseamount, \
                        linewidth=0.05, color = cm(slides/NUM_COLORS))
                        baseamount = baseamount + int(peptides[slides+7])/totalcontrol
                    if totalinduced >0:
                        plt.bar(left = int(peptides[5]), width = int(peptides[6]),\
                        height = -int(peptides[slides+7+12])/totalinduced, bottom = baseamount2, \
                        linewidth=0.05, color = cm(slides/NUM_COLORS))
                        baseamount2 = baseamount2 - int(peptides[slides+7+12])/totalinduced
            import matplotlib.patches as mpatches
            myhandles = []
            for slides in range(NUM_COLORS):
                myhandles.append(mpatches.Patch(color = cm(slides/NUM_COLORS), label = str(slides+1)))
            plt.legend(handles=myhandles,loc=2,borderaxespad=0.,bbox_to_anchor=(1, 1),fontsize=4)
#            plt.legend(handles=myhandles,markerscale=0.3,fontsize=4)
        pp.savefig()
    plt.close()
pp.close()
exit()





    
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
cm = plt.get_cmap('Accent')
NUM_COLORS = 12

myhandles = []
for slides in range(NUM_COLORS):
    myhandles.append(mpatches.Patch(color = cm(slides/NUM_COLORS), label = str(slides+1)))
red_patch = mpatches.Patch(color='red', label='The red data')
green_patch = mpatches.Patch(color='green', label='The green data')

plt.legend(handles=myhandles,loc=2,borderaxespad=0.)

plt.show()    