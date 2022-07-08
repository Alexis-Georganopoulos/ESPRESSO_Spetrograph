# -*- coding: utf-8 -*-
"""
Created on Mon May 31 08:45:25 2021

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt


do_i_plot = False
#PSF plots raw
n_w = 9
n_c = 12
rootpath = 'C:/Users/User/Documents/EPFL/Robotics-2019_2020-Year 3/Masters Thesis/zemax_psf/'
w_range = [i+1 for i in range(n_w)]
lst = []

for j in range(n_c):
    for i in range(n_w):
        lst.append(rootpath+'C'+str(j+1)+'_W'+str(i+1)+'.txt')
    
del i
del w_range
del rootpath


psf_all = []
if do_i_plot:
    plt.figure(figsize = (18,6))
for i,f in enumerate(lst):
    data = np.loadtxt(f,skiprows=17)
    psf_all.append(data)
    if do_i_plot:
        plt.subplot(n_c,n_w,1+i)
        plt.imshow(data,extent = [-150,150,-150,150])
    
    
#del n_w
#del n_c
del i
del j
del f 
del data
del lst
#plt.close('all')
#%%Dispersion plots
rootpath = 'C:/Users/User/Documents/EPFL/Robotics-2019_2020-Year 3/Masters Thesis/zemax_psf/RayTrace_Iterate.txt'
data = np.loadtxt(rootpath, skiprows = 9)
orders = [((i)*n_w) for i in range(n_c)]
#data format:
#ORDER---WAVELENGTH(um)---XPOS(mm)---YPOS(mm)

order = orders[0]
y = data[order:order+9, 1]
x = data[order:order+9, 2]

#8th degree polynomial creates the required accuracy
if do_i_plot:
    plt.figure()
    plt.grid('on')
    plt.plot(x,y, '-*')



n_p = []
for i in range(n_c):
    order = orders[i]
    x = data[order:order+n_w, 2]
    y = data[order:order+n_w, 1]
    p = np.polyfit(x, y, 8)
    n_p.append(p)
    if do_i_plot:
        plt.figure()
        plt.plot(x,y, '*', label  = 'Ray-Traced data')
    yp = np.polyval(n_p[i], x)
#    print(min(x), max(x))
    if do_i_plot:
        #plt.plot(x,yp,'-*')
        pass
    x = np.linspace(min(x), max(x), 100)#just to check interpolation isnt going insane
    yp = np.polyval(n_p[i], x)
    if do_i_plot:
        plt.plot(x,yp, label = '8^th degree interpolation')
        plt.legend()
        plt.title('Ray-Traced and Interpolated dispersion law for order '+str(i+1))
        plt.xlabel('Position x on detector [mm]')
        plt.ylabel('Wavelength [nm]')
        plt.grid('on')

del rootpath
#del data
del i
del order
del orders
del p
del x
del y
del yp
#plt.close('all')
# p = np.polyfit(x, y, 8)

# yp = np.polyval(p,x)

# plt.plot(x,(yp-y)/y, '-*')

#%%Kernel generation & centroid calculation
psf_x_collapsed = []
for i in range(n_c*n_w):
    psf_x_collapsed.append(np.sum(psf_all[i], axis=0))
    
psf_x_centroids = []
for i in range(n_c*n_w):
    weightedx = 0
    for j in range(len(psf_x_collapsed[i])):
        weightedx = weightedx + j*psf_x_collapsed[i][j]
    psf_x_centroids.append(weightedx/(psf_x_collapsed[i].sum()))
    
del i
del j
del weightedx
#♣just to see whats happening
x = [i for i in range(100)]
y = psf_x_collapsed[107]

if do_i_plot:
    
    plt.figure()
    plt.plot(x,y)
    plt.grid('on')
    plt.xlabel('Pixel Position x')
    plt.ylabel('Relative luminal intensity')
    plt.title('Kernel for order 12, wavelength 9')
del x
del y
#plt.close('all')

#%%adjusting the kernels(centering them), and inputting them into the correct coordinates

#mm per pixel
px = 0.3/100
psf_x_centroids_mm = psf_x_centroids
for i in range(n_w*n_c):
    psf_x_centroids_mm[i] = np.arange(1- psf_x_centroids[i],101- psf_x_centroids[i],1)
    psf_x_centroids_mm[i] = psf_x_centroids_mm[i] * px
    psf_x_centroids_mm[i] += data[i][2]
    
#psf_x_centroids_mm is the de-meaned version of psf_x_centroids where each element is an array, and all converted into mm
#we use this with the diespersion law to turn these mm length into nm wavelengths(since thats what wwe need for the kernel)
#therefore, we add the offset from the original data(x) to find the correct x coordinate
psf_w_collapsed = []

for i in range(n_c):
    p = n_p[i]
    
    for j in range(n_w):
        x = psf_x_centroids_mm[i*9+j]
        yp = np.polyval(p,x)       
        psf_w_collapsed.append(yp)
del i
del j
del p
del x
del yp

if do_i_plot:
    plt.figure()
    plt.plot(psf_w_collapsed[-1]*(1e3), psf_x_collapsed[-1])
    plt.grid('on')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Relative luminal intensity')
    plt.title('Kernel for order 12, wavelength 9')

#%%Formatting the data to use it in other modules(avoiding to save as text file)
ordered_data = []

for i in range(n_c):
    for j in range(n_w):
        
        wavelength = psf_w_collapsed[i*9+j]
        intensity =  psf_x_collapsed[i*9+j]
        
        weightedx = 0
        for k in range(100):
            weightedx = weightedx + intensity[k]*wavelength[k]
            
        centroid_w = weightedx/(intensity.sum())
        
        ordered_data.append([i+1, centroid_w , wavelength, intensity])
        #                       data[i*9+j][1]
        #•print(centroid_w, data[i*9+j][1])



del i
del j
del k
del weightedx
del centroid_w

#%%clean up useless data when this module is exported
del data
del do_i_plot
del intensity
#del n_c
del n_p
#del n_w
del px
del wavelength
del psf_all
del psf_w_collapsed
del psf_x_collapsed
del psf_x_centroids
del psf_x_centroids_mm

#%% convert to m and adjusting colapsed psf to 1
for count, case in enumerate(ordered_data):
    ordered_data[count][1] *= 1e-6#convert micrometers to meters 
    ordered_data[count][2] *= 1e-6#convert micrometers to meters 
    ordered_data[count][3] /= max( ordered_data[count][3])#scaling the gaussians to a peak value of 1

del count
del case
