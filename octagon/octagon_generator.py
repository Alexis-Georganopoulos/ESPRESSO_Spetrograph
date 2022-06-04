# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 09:38:30 2021

@author: User
"""



import numpy as np
from PIL import Image
import os

image_file = Image.open(r"C:\Users\User\Documents\EPFL\Robotics-2019_2020-Year 3\Masters Thesis\octagon-deluxe-spa-cover-50.png") # open colour image
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
# tol = 1e-6
# theta_octagon = np.pi/4

# #counterclockwise rotation, but makes no difference in our deffinition
# R = np.array([[np.cos(theta_octagon), -np.sin(theta_octagon)], [np.sin(theta_octagon), np.cos(theta_octagon)]])

# #♣keep this at 1 for simplicity
# square_side_length = 1

# divisions = 100

# #the step in our coordinate system, also coupoled to the resolution of the pixels of the grid
# #1/coord_step_per_pixel = number of grid pixels per unit length
# coord_step_per_pixel = square_side_length/divisions

# x_square_coord = np.arange(-square_side_length/2, square_side_length/2+coord_step_per_pixel, coord_step_per_pixel)

# #clean up python rounding errors
# for i in range(len(x_square_coord)):
#     x_square_coord[i] = round(x_square_coord[i],int(np.ceil(np.log10(divisions))))
# x_square_coord[int(1/coord_step_per_pixel/2)] = 0
# y_square_coord = np.flip(x_square_coord)



# square_coord = np.ndarray(shape = (len(x_square_coord),len(y_square_coord), 2))
# square_coord[::]=0

# for i in range(len(x_square_coord)):
#     for j in range(len(y_square_coord)):
#         square_coord[j][i] = [x_square_coord[i],y_square_coord[j]] 
        
# rhombus_coord = np.ndarray(shape = (len(x_square_coord),len(y_square_coord), 2))
# rhombus_coord[::]=0

# for i in range(len(x_square_coord)):
#     for j in range(len(y_square_coord)):
#         rhombus_coord[j][i] = R@square_coord[j][i]
#         # if abs(rhombus_coord[j][i][0]) < tol:
#         #     rhombus_coord[j][i][0] = 0
#         # if abs(rhombus_coord[j][i][1]) < tol:
#         #     rhombus_coord[j][i][1] = 0
        
# #def coord_to_grid(coord_matrix, pixels_per_coord_step)

# del square_side_length
# del theta_octagon
# del divisions
# del i
# del j