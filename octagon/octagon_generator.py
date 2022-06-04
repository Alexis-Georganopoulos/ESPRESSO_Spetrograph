# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 09:38:30 2021

@author: User
"""


from os import path
import numpy as np
from PIL import Image
#C:\Users\User\Documents\EPFL\Robotics-2019_2020-Year 3\Masters Thesis\
image_source = path.dirname(__file__)

image_file = Image.open(image_source+'\octagon-deluxe-spa-cover-50.png') # open colour image
thresh = 200
fn = lambda x : 255 if x > thresh else 0
image_file = image_file.convert('L').point(fn, mode='1')


image_file.show()

mtrx = np.array(image_file)
W,L= np.shape(mtrx)[0],np.shape(mtrx)[1]

top, bottom = 0 ,0
left, right = 0 ,0

current_pixel = True
for i in range(L):
    if mtrx[i][int(W/2)] != current_pixel:
        if i < int(L/2):
            top = i
            current_pixel = False
        else:
            bottom = i
            break

current_pixel = True
for i in range(W):
    if mtrx[int(L/2)][i] != current_pixel:
        if i < int(W/2):
            left = i
            current_pixel = False
        else:
            right = i
            break
        
image_file = image_file.crop((left,top,right,bottom))
image_file.show()
mtrx = np.array(image_file)

newsize = (101,101)
image_file = image_file.resize(newsize)

docs_dir=os.path.expanduser('~\Documents\EPFL\Robotics-2019_2020-Year 3\Masters Thesis')
image_file.save(os.path.join(docs_dir,'OCTAGON.png'))

image_file.show()
mtrx = np.array(image_file)
mtrx = (~mtrx)*1


np.savetxt(os.path.join(docs_dir,'OCTAGON.txt'),mtrx, fmt='%d',delimiter='')
del i
del fn
