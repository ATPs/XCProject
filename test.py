
txt = open('list.txt').read()
import re
ls = re.split('[\n:\t]',txt)
lsnt = [e for e in ls if re.findall('.gz$',e) and re.findall('^nr',e)]
fo = open('list.txt','w')
for e in lsnt:
    fo.write('wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/'+e+'\n')
fo.close()

from ete3 import Tree
shape = []
for i in range(10000):
    t1 = Tree()
    t1.populate(672)
    _a = t1.get_closest_leaf()[1]
    _b = t1.get_farthest_leaf()[1]
    shape.append((_a,_b))
import matplotlib.pyplot as plt
import numpy as np
x = np.array([e[0] for e in shape])
y = np.array([e[1] for e in shape])
plt.plot(x, y,'o')

import subprocess
subprocess.run(r'''"C:\Users\ATPs\OneDrive\Lab\YangLab\201709BackgroundInfo\v0.3\windows\HSA.exe"  "C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\Lineages\fun.alm" "C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\Lineages\fun.alm" "C:\Users\ATPs\OneDrive\Lab\YangLab\YangLabSharedProject_Xiaolong\Celegans\Lineages\cost.tsv" l 100''')