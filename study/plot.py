# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:15:03 2015

@author: k
"""

import numpy as np
import matplotlib.pyplot as plt

# evenly sampled time at 200ms intervals
t = np.arange(0., 5., 0.2)

# red dashes, blue squares and green triangles
plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

mu, sigma = 100, 15
x = mu + sigma * np.random.randn(10000)

# the histogram of the data
n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)


plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title('Histogram of IQ')
plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.axis([40, 160, 0, 0.03])
plt.grid(True)
plt.show()


import numpy as np
import matplotlib.pyplot as plt

p1= plt.bar(range(10),[5 for i in range(10)],width=0.5,color=['0.8',"r","r","y"],linewidth=0)
p2= plt.bar([1,2,3,4,5],[3,4,3,3,3],width=[1,0.25,1.25,0.5,0.5],color='g',bottom = [5,5,5,5,5])
p3= plt.bar([1,2,3,4,5],[-3,-4,-3,-3,-3],width=1,color='g')
p3= plt.bar(0,5,width=7,color='c',bottom = 6)
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('multipage.pdf')
pp.savefig()
pp.savefig()
p4= plt.bar(0,5,width=7,color='m',bottom = 12)
pp.savefig()
pp.close()


mylist = []
templist1 = open("test.txt").readlines()
for ele in templist1:
    mylist.append(ele.split())



import matplotlib.pyplot as plt
import numpy as np
x = list(range(10))
y = np.random.randint(1,100,10)
fig = plt.figure()
plt.vlines(x,0,y)
