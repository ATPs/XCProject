# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:09:43 2016

@author: k
"""

import numpy as np
check = np.zeros((9,9))
check[::2,1::2] = 1
check[1::2,::2] = 1
import matplotlib.pylab as plt
plt.imshow(check,cmap='gray',interpolation = 'nearest')

import skimage
from skimage import data
camera = data.camera()
from skimage import restoration
filtered_camera = restoration.denoise_bilateral(camera)

import os
from skimage import io
filename = os.path.join(skimage.data_dir,'camera.png')
camera = io.imread(filename)

camera = data.camera()
from skimage import img_as_float
camera_float = img_as_float(camera)


import scipy
face = scipy.misc.face()
face.shape

dimen = 1
a = np.zeros((dimen,dimen,3))
a[::2,::2,2] = 100
plt.imshow(a)

import random
a = [[[i/255,i/255,i/255] for i in range(0,255,1)] for i in range(200)]
plt.imshow(a,interpolation = "nearest")


import os
from skimage import io
filename = "E:\\store_for_D\\Transcriptome\\KEGG_Annotation2016\\test\\KEGG PATHWAY_ Glycolysis _ Gluconeogenesis - Reference pathway + ame api aga dme tca bmor_files\\map00010.png"
fig = io.imread(filename)
scipy.misc.imsave('outfile.png', fig)

fig[398:415,460:502] = [255,255,0]
scipy.misc.imsave('outfile.png', fig)

testp = fig[169:186,460:502]
