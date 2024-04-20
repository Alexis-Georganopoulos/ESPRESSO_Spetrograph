# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 13:33:54 2021

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.signal as sci

from tqdm import tqdm
import os

import psf_analyzer as PSF

ordered_data = PSF.ordered_data#units in m
n_wavelengths = PSF.n_w#the number of wavelengths in a single order
n_configurations = PSF.n_c#the number of orders

#%%
do_i_plot = True

virtual_steps = 1000000 #just how many increments we want for some range of numbers

gaussian_sample_space = 130000#also temporary, depending on what we want
gauss_cuttoff = 0.001 #the percentage height at which ignore the rest of the gaussian

etalon_spacing = 7.6e-3 ##+/- 5e-7 m

c = 3e8

misalignment_container = np.zeros((n_configurations,n_wavelengths))

def plotter_function(xvalues, yvalues, title = '', xlabel = '', ylabel = '', labels = ''):
    global do_i_plot
    if do_i_plot == False:
        return
    if len(xvalues) != len(yvalues):
        raise Warning('Not enough x/y pairs')
        return
    n  = len(xvalues)
    plt.figure()
    
    if type(xvalues[0]) in [int, float, np.float64, np.int64]:
        plt.plot(xvalues, yvalues)
    else:
        for i in range(n):
            temp = labels if len(labels) == 0 else labels[i]
            plt.plot(xvalues[i], yvalues[i], label = temp)
        if labels:
            plt.legend()
    
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.show()
 
def white_light_generator(lambda_min, lambda_max, lambda_divisions):
    difference = lambda_max - lambda_min
    increment = difference/lambda_divisions
    return np.arange(lambda_min, lambda_max+ increment, increment)

def gaussian(cuttoff, sample_space, wdivisions, wincr, target):
    
    gaussian_resolution = target/sample_space #the FWHM we want
    sigma_guassian = gaussian_resolution/(2*np.sqrt(2*np.log(2)))

    xrange = np.arange(-int(wdivisions/2)*wincr, \
                       (int(wdivisions/2)+1)*wincr, wincr)
    #mag_gauss = (1/(np.sqrt(2*np.pi)*sigma_guassian))
    exp_gauss = np.exp(-(xrange**2)/(2*sigma_guassian**2))
    indexes = exp_gauss > cuttoff*np.max(exp_gauss)
    exp_gauss = exp_gauss[indexes]
    #gauss = mag_gauss*exp_gauss
    xrange = xrange[indexes]
    # plt.figure()
    # plt.plot(xrange, exp_gauss)
    # plt.show()
    return exp_gauss, xrange#return gauss, xrange

##We assume for now that the finesse is constant for all lambda
## in fact finesse changes by the incedent angle for each wavelength(not reflectance)
def fabry_perot_transmitance(lambda_range, etalon_thickness, \
                             n = 1, theta = 0, F = 12.27):
    #F = 4*R/(1-R**2)
    delta = (2*np.pi/lambda_range)*2*n*etalon_thickness*np.cos(theta)
    product = F*(np.sin(delta/2))**2
    return 1/(1 + product) ##The transmitance

def peakfinder(signal, axis):#signal and axis MUST be the same dimenstion and correspond to eachother(order matters)
    peak_idx = sci.find_peaks(signal)[0]
    axis_values_of_peaks = axis[peak_idx]
    return axis_values_of_peaks

def peak_allignment(sample, reference):
    global c
    aligned_peaks = sample - reference
    aligned_peaks /= reference
    aligned_peaks *= c
    return aligned_peaks

def reject_false_peaks(array1, array2):#tol is the tolerance below which we reject a peak

    temp1 = array1
    temp2 = array2
    #tol = 5e-11 #this is very arbitrary and happens to work in this application
    # temp = []
    
    # for i in range(len(temp1)-1):
    #     if (temp1[i+1] - temp1[i])
        
    a = temp1[0] - temp2[0]
    b = temp1[-1] - temp2[-1]
    if abs(a) > abs(b):
        if a > 0:
            return[temp1, temp2[1:]]
        else:
            return[temp1[1:], temp2]
    else:
        if b > 0:
            return[temp1[0:-1], temp2]
        else:
            return[temp1, temp2[0:-1]]

def median_span(array, below = 0, above = 0):
    medianspan = []
    temp = []
    
    n = len(array)
    
    if n%2 == 1:
        temp = array
        n = (n+1)/2
    else:
        counter = 0
        while(counter < n - 1):
            temp2 = (array[counter] + array[counter + 1])/2
            temp.append(temp2) 
            counter += 1
            
        n = n/2
   
    n = int(n) -1#python indexing starts at 0, not 1(so we sub 1)
    counter = n - below 
    while(counter <= n + above):
        medianspan.append(temp[counter])
        counter += 1
        
    return medianspan

try:
    
    print('Loading misalignment matrix')
    #misalignment_container = np.loadtxt('misalignment_matrix.txt')
    #misalignment_container = np.loadtxt('laser_matrix.txt')
    misalignment_container = np.loadtxt('test.txt')
    print('Succesefully loaded')
except:
    print("Misalignment matrix not found in current folder, redo from scratch?(y/n)")
    redo = input()
    if redo not in ['y','Y']:
        raise Warning('aborted')

    for i in [107]:#tqdm(range(0,n_wavelengths*n_configurations,1)):#debug cases [106]=c12_w8 weird for lazer
        
        ordered_data_idx = i#i#55   
        
        #%%
        wtest = ordered_data[ordered_data_idx][1]
        xtest = ordered_data[ordered_data_idx][2]
        ytest = ordered_data[ordered_data_idx][3]
        
        wtest_string = str(round(1e9*wtest, 0))
        
        plotter_function(xtest, ytest, "True Kernel(Not normalised) @ "+wtest_string+" nm", "Centered wavelength [nm]")
        
        #%%
        
        lambda_target = wtest #[m]
        lambda_deviation = 0.0005*wtest#0.3e-9 #[m]this varies depending on application    0.3e-9
        
        lambda_min = lambda_target - lambda_deviation
        lambda_max = lambda_target + lambda_deviation
        
        increment = (lambda_max-lambda_min)/virtual_steps
        
        whitelight = white_light_generator(lambda_min, lambda_max, virtual_steps)
        
        #%%
        
        
        gauss, dummy_x = gaussian(gauss_cuttoff, gaussian_sample_space, virtual_steps, increment, lambda_target)
        
        plotter_function([1e9*(dummy_x + wtest), 1e9*xtest], [gauss, ytest],\
                         title = "Unit Gaussian(Not normalised) @ "+wtest_string+" nm",\
                             xlabel = "Centered wavelength [nm]")
        
        f = interp.interp1d(xtest, ytest)
        compatible_res_axis = dummy_x + wtest
        #interpolating the new gaussian in order to compensate for the different resolution
        true_kernel = f(compatible_res_axis)    
        
        
        
        #%%
        
        #this also acts as the reference transmitance, not just for composing systems
        fabry_perot = fabry_perot_transmitance(whitelight, etalon_spacing)
        
        laser_comb = fabry_perot_transmitance(whitelight, etalon_spacing, F = 10000)
        
        plotter_function(1e9*whitelight, fabry_perot, \
                         title = "Transmitance of the Fabry-Perot @ "+wtest_string+" nm",\
                         ylabel = "Transmitance", xlabel = "Wavelength [nm]")
            
        plotter_function(1e9*whitelight, laser_comb, \
             title = "Transmitance of the Laser-Comb @ "+wtest_string+" nm",\
             ylabel = "Transmitance", xlabel = "Wavelength [nm]")
        
        #%% convolutions
        # ideal_convolution = np.convolve(gauss, fabry_perot, mode = 'same')
        # ideal_convolution /= ideal_convolution.max()
        
        # true_convolution = np.convolve(true_kernel, fabry_perot, mode = 'same')
        # true_convolution /= true_convolution.max()
        
        # plotter_function([1e9*whitelight,1e9*whitelight],[ideal_convolution, true_convolution],\
        #                       labels = ['Ideal convolution','Real convolution'],\
        #                       xlabel="Wavelength [nm]", ylabel = "Transmitance", \
        #                       title = "Resultant convolutions @ "+wtest_string+" nm")
            
        fb_convolution = np.convolve(true_kernel, fabry_perot, mode = 'same')
        fb_convolution /= fb_convolution.max()

        las_convolution = np.convolve(true_kernel, laser_comb, mode = 'same')
        las_convolution /= las_convolution.max()

        plotter_function([1e9*whitelight,1e9*whitelight],[fb_convolution, las_convolution],\
                              labels = ['Fabry-Perot convolution','Laser-Comb convolution'],\
                              xlabel="Wavelength [nm]", ylabel = "Transmitance", \
                              title = "Resultant convolutions @ "+wtest_string+" nm")
        
        #%%finding the peaks and alligning them(we reject the extrema due to convolutionnuisances)
    
        
        # ideal_peaks = peakfinder(ideal_convolution, whitelight)
        # true_peaks = peakfinder(true_convolution, whitelight)
       
        # # plt.plot(1e9*ideal_peaks, 0.6*np.ones(len(ideal_peaks)), '*')
        # # plt.plot(1e9*true_peaks, 0.4*np.ones(len(true_peaks)), '*') 
       
        # if ideal_peaks.size != true_peaks.size:
        #     ideal_peaks, true_peaks = reject_false_peaks(ideal_peaks, true_peaks)
        # # plt.plot(1e9*ideal_peaks, np.ones(len(ideal_peaks)), '*')
        # # plt.plot(1e9*true_peaks, 0.8*np.ones(len(true_peaks)), '*')
        # misalignment = peak_allignment(true_peaks, ideal_peaks)
        
        # plotter_function(1e9*ideal_peaks, misalignment, \
        #                  title = "Radial speed @ "+wtest_string+" nm",\
        #                  ylabel = "Error in speed[m/s]", xlabel = "Wavelength [nm]")
        ideal_peaks = peakfinder(fb_convolution, whitelight)
        true_peaks = peakfinder(las_convolution, whitelight)
       
        # plt.plot(1e9*ideal_peaks, 0.6*np.ones(len(ideal_peaks)), '*')
        # plt.plot(1e9*true_peaks, 0.4*np.ones(len(true_peaks)), '*') 
       
        if ideal_peaks.size != true_peaks.size:
            ideal_peaks, true_peaks = reject_false_peaks(ideal_peaks, true_peaks)
        # plt.plot(1e9*ideal_peaks, np.ones(len(ideal_peaks)), '*')
        # plt.plot(1e9*true_peaks, 0.8*np.ones(len(true_peaks)), '*')
        misalignment = peak_allignment(true_peaks, ideal_peaks)
        
        plotter_function(1e9*ideal_peaks, misalignment, \
                         title = "Radial speed @ "+wtest_string+" nm",\
                         ylabel = "Error in speed[m/s]", xlabel = "Wavelength [nm]")
        
        
        
        #%% Taking the center 3 misalignment points and averaging them
        
        focused_average_misalignment = np.mean(median_span(misalignment,1,1))
        
        misalignment_container[i//n_wavelengths][i%n_wavelengths] = focused_average_misalignment


del c
del gauss_cuttoff
del gaussian_sample_space
del etalon_spacing
del virtual_steps
try:
    del dummy_x
    del ordered_data_idx
    del misalignment
    del compatible_res_axis
    del fabry_perot
    del laser_comb ####################
    del focused_average_misalignment
    del gauss
    del i
    del fb_convolution
    #del ideal_convolution
    del ideal_peaks
    del increment
    del lambda_deviation
    del lambda_max
    del lambda_min
    del lambda_target
    del las_convolution
    #del true_convolution
    del true_kernel
    del true_peaks
    del whitelight
    del wtest
    del wtest_string
    del xtest
    del ytest
    del redo
    
except:
    pass

#%%Save the result
#np.savetxt('laser_matrix.txt', misalignment_container)

#%%misalignment by order plots

wavelength_container = np.zeros((n_configurations,n_wavelengths))

for i in range(n_wavelengths*n_configurations):
    wavelength_container[i//n_wavelengths][i%n_wavelengths] = ordered_data[i][1]

plt.figure()
for i in range(n_configurations):
    plt.plot(wavelength_container[i]*(1e9), misalignment_container[i], label = 'order'+str(i+1))
    
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radial Speed [m/s]')
plt.grid()
plt.title('Radial Speed between the true Kernel Fabry-Perot and Laser-Comb')

del i

