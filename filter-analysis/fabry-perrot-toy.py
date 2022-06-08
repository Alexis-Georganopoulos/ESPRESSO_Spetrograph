# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:13:34 2021

@author: Alexis Philip George-Georganopoulos
"""
##All units are in SI in base units, conversions are done just for plotting
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sci
import OpticalSystems as opsys


virtual_steps = 1000000 #just how many increments we want for some range of numbers
lambda_target = 600e-9 #[m]
lambda_deviation = 0.3e-9 #[m]this varies depending on application    0.3e-9

lambda_min = lambda_target - lambda_deviation
lambda_max = lambda_target + lambda_deviation

increment = (lambda_max-lambda_min)/virtual_steps

whitelight = opsys.white_light_generator(lambda_min, lambda_max, virtual_steps)

#%%Light filter before entering fabry_perot

filter_thickness = 2e-3#1e-3 #temporary, we would like to change this
my_tilt = 2# +ve ->counter-clockwise tilt, -ve ->clockwise tilt

transmitance_filter = opsys.filter_transmitance(whitelight, filter_thickness, \
                                           n = 1.5, percentage=0.05)
accurate_transmitance_filter = opsys.accurate_filter_transmitance(whitelight, filter_thickness, path_length=43.718e-3,\
                                                            semi_diameter=0.993e-3, n = 1.5, R1= 0.63947, tilt_deg= my_tilt)    
    
plt.figure()
#plt.plot(1e9*whitelight, transmitance_filter, label = 'Simplified Filter')
plt.plot(1e9*whitelight, accurate_transmitance_filter, label = 'Realistic Filter')
plt.legend()
plt.title("Transmitance of the light filter with "+str(my_tilt)+" degrees tilt")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmitance")
plt.grid()
plt.show()

del filter_thickness
del my_tilt
#%% fabry-perrot
etalon_spacing = 7.6e-3 ##+/- 5e-7 m

#this also acts as the reference transmitance, not just for composing systems
transmitance_fabry_perot = opsys.fabry_perot_transmitance(whitelight, etalon_spacing)

plt.figure()
plt.plot(1e9*whitelight, transmitance_fabry_perot)
plt.title("Transmitance of the Fabry-Perot")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmitance")
plt.grid()
plt.show()

perfect_lambda = opsys.analytic_fabry_perot_peaks(lambda_min, lambda_max, etalon_spacing)

del etalon_spacing
#%%Neutral Density Filter

#Assuming normal incidence, air->Ni-Fe
R1 = abs((1-2.2714)/(1+2.2714))**2
neutral_density_transmitance = opsys.neutral_density(whitelight, R1)
#%%filter + fabry perot
#(system composition)
fabry_perot = transmitance_fabry_perot#redundant, but the wording is consistent
fabry_perot_filter = transmitance_fabry_perot*transmitance_filter
accurate_fabry_perot_filter = transmitance_fabry_perot*accurate_transmitance_filter

fig, axs = plt.subplots(2, 1, sharex = 'all', sharey='all')
axs[0].plot(1e9*whitelight, fabry_perot, color = 'r')
axs[0].set_title("Transmitance of the Fabry-Perot")
# axs[1].plot(1e9*whitelight, fabry_perot_filter)
# axs[1].set_title("Transmitance of the Fabry-Perot + simple filter")
axs[1].plot(1e9*whitelight, accurate_fabry_perot_filter, color = 'g')
axs[1].set_title("Transmitance of the Fabry-Perot + realistic filter")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmitance")
axs[0].grid()
axs[1].grid()
#axs[2].grid()
plt.show()

#%%gaussian generation

gaussian_sample_space = 130000#also temporary, depending on what we want
gauss_cuttoff = 0.001 #the percentage height at which ignore the rest of the gaussian


gauss, dummy_x = opsys.gaussian(gauss_cuttoff, gaussian_sample_space, lambda_target, \
                          lambda_max, lambda_min, virtual_steps)
plt.figure()
plt.plot(1e9*dummy_x, gauss)
plt.title("Unit Gaussian(Not normalised)")
plt.xlabel("Centered wavelength [nm]")
plt.ylabel("")
plt.grid()
plt.show()

del gauss_cuttoff
del gaussian_sample_space
del dummy_x

#%%Convolution
pure_convolution = np.convolve(gauss, fabry_perot, mode = 'same')
pure_convolution /= pure_convolution.max()

simple_convolution = np.convolve(gauss, fabry_perot_filter, mode = 'same')
simple_convolution /= simple_convolution.max()

accurate_convolution = np.convolve(gauss, accurate_fabry_perot_filter, mode = 'same')
accurate_convolution /= accurate_convolution.max()

fig, axs = plt.subplots(2, 1, sharex = 'all', sharey='all')
axs[0].plot(1e9*whitelight, pure_convolution, color = 'r')
axs[0].set_title("Transmitance of the Fabry-Perot after convolution")
# axs[1].plot(1e9*whitelight, simple_convolution)
# axs[1].set_title("Transmitance of the Fabry-Perot + simple filter after convolution")
axs[1].plot(1e9*whitelight, accurate_convolution, color = 'g')
axs[1].set_title("Transmitance of the Fabry-Perot + realistic filter after convolution")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmitance")
axs[0].grid()
axs[1].grid()
#axs[2].grid()
plt.show()

#%%Finding peaks before discretization

#python peak finder
fb_peaks = sci.find_peaks(fabry_perot)[0]
#analytic peaks
perfect_lambda;


pure_peaks = sci.find_peaks(pure_convolution)[0] 
simple_peaks = sci.find_peaks(simple_convolution)[0]
acc_conv_peaks = sci.find_peaks(accurate_convolution)[0]

#just a quick check to see things are ok
plt.figure()
plt.plot(1e9*whitelight, fabry_perot)
plt.plot(1e9*whitelight[fb_peaks], fabry_perot [fb_peaks],'*')
#plt.plot(1e9*whitelight, convolution)
plt.plot(1e9*whitelight, accurate_convolution)
plt.plot(1e9*whitelight[acc_conv_peaks], accurate_convolution[acc_conv_peaks],'*')
plt.title("Convoluted Fabry-perot + realistic filter, naive peak-finding")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmitance")
plt.grid()
plt.show()

# fig, axs = plt.subplots(3, 2, sharex = 'all', sharey='all')
# axs[0,0].plot(1e9*whitelight[fb_peaks[0]], fabry_perot [fb_peaks[0]],'*')
# axs[0,0].plot(1e9*whitelight, pure_convolution)
# axs[0,0].plot(1e9*whitelight[pure_peaks[0]], pure_convolution [pure_peaks[0]],'*')
# axs[0,0].set_title("Transmitance of the Fabry-Perot after convolution, naive peak-finding")
# axs[1].plot(1e9*whitelight, simple_convolution)
# axs[1].set_title("Transmitance of the Fabry-Perot + simple filter after convolution")
# axs[2].plot(1e9*whitelight, accurate_convolution, color = 'g')
# axs[2].set_title("Transmitance of the Fabry-Perot + realistic filter after convolution")
# plt.xlabel("Wavelength [nm]")
# plt.ylabel("Transmitance")
# for i in range(3):
#     for j in range (2):
#         axs[i,j].grid()
# plt.show()

# del i
# del j
#%%plotting the difference between fabry-perot and the convolution
c = 3e8

pure_aligned_peaks = (pure_peaks-fb_peaks)*((lambda_max-lambda_min)/virtual_steps)
pure_aligned_peaks[1:-1] /= whitelight[fb_peaks[1:-1]]
pure_aligned_peaks *= c

pure_aligned_peaks_analytic = whitelight[pure_peaks] - perfect_lambda
pure_aligned_peaks_analytic /= perfect_lambda
pure_aligned_peaks_analytic *= c

simple_aligned_peaks = (simple_peaks-fb_peaks)*((lambda_max-lambda_min)/virtual_steps)
simple_aligned_peaks[1:-1] /= whitelight[fb_peaks[1:-1]]
simple_aligned_peaks *= c

simple_aligned_peaks_analytic = whitelight[simple_peaks] - perfect_lambda
simple_aligned_peaks_analytic /= perfect_lambda
simple_aligned_peaks_analytic *= c

acc_aligned_peaks = (acc_conv_peaks-fb_peaks)*((lambda_max-lambda_min)/virtual_steps)
acc_aligned_peaks[1:-1] /= whitelight[fb_peaks[1:-1]]
acc_aligned_peaks *= c

acc_aligned_peaks_analytic = whitelight[acc_conv_peaks] - perfect_lambda
acc_aligned_peaks_analytic /= perfect_lambda
#TODO
#FABRY-PEROT+INFRARED REFERENCE
acc_aligned_peaks_analytic *= c

# plt.figure()
# plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_aligned_peaks[1:-1])
# plt.plot(1e9*whitelight[fb_peaks[1:-1]], simple_aligned_peaks[1:-1])
# plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_aligned_peaks[1:-1])
# plt.title("Error between fb and conv")
# plt.xlabel("Wavelength [nm]")
# plt.ylabel("Error in speed[m/s]")
# plt.grid()
# plt.show()

#in orange are the analytic peak estimations, in blue are the python peakfinders
fig, axs = plt.subplots(2, 1, sharex = 'all', sharey='all')
#axs[0].plot(1e9*whitelight[fb_peaks[1:-1]], pure_aligned_peaks[1:-1])
axs[0].plot(1e9*perfect_lambda[1:-1], pure_aligned_peaks_analytic[1:-1])
axs[0].set_title("Error of convoluted Fabry-Perot before discretisation")
#axs[1].plot(1e9*whitelight[simple_peaks[1:-1]], simple_aligned_peaks[1:-1])
# axs[1].plot(1e9*perfect_lambda[1:-1], simple_aligned_peaks_analytic[1:-1])
# axs[1].set_title("Error of convoluted Fabry-Perot + simple filter before discretisation")
#axs[2].plot(1e9*whitelight[acc_conv_peaks[1:-1]], acc_aligned_peaks[1:-1])
axs[1].plot(1e9*perfect_lambda[1:-1], acc_aligned_peaks_analytic[1:-1])
axs[1].set_title("Error of convoluted Fabry-Perot + realistic filter before discretisation")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
axs[0].grid()
axs[1].grid()
#axs[2].grid()
plt.show()

# del c
# del pure_aligned_peaks
# del pure_aligned_peaks_analytic
# del simple_aligned_peaks
# del simple_aligned_peaks_analytic
# del acc_aligned_peaks
# del acc_aligned_peaks_analytic

#%%Discretisation
wavelength_resolution = lambda_target/130000/3

pure_discretized_convolution, pure_edges = opsys.discretize(whitelight, pure_convolution, \
                                            wavelength_resolution, lambda_max, lambda_min)

discretized_convolution, edges = opsys.discretize(whitelight, simple_convolution, \
                                            wavelength_resolution, lambda_max, lambda_min)
acc_discretized_convolution, acc_edges = opsys.discretize(whitelight, accurate_convolution, \
                                            wavelength_resolution, lambda_max, lambda_min)

plt.figure()
plt.plot(1e9*whitelight, pure_convolution)
plt.plot(1e9*pure_edges, pure_discretized_convolution, 'ro')
plt.bar(1e9*pure_edges,pure_discretized_convolution, alpha = 0.5, color ='orange', \
        width = 1e9*wavelength_resolution, align = 'center')
plt.title("Convoluted Fabry-Perot against discretized points")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmittance")
plt.grid()
plt.show()
    
plt.figure()
plt.plot(1e9*whitelight, simple_convolution)
plt.plot(1e9*edges, discretized_convolution, 'ro')
plt.bar(1e9*edges,discretized_convolution, alpha = 0.5, color ='orange', \
        width = 1e9*wavelength_resolution, align = 'center')
plt.title("Convoluted Fabry-Perot + simple filter against discretized points")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmision")
plt.grid()
plt.show()
#plt.close()

plt.figure()
plt.plot(1e9*whitelight, accurate_convolution)
plt.plot(1e9*acc_edges, acc_discretized_convolution, 'ro')
plt.bar(1e9*acc_edges,acc_discretized_convolution, alpha = 0.5, color ='orange', \
        width = 1e9*wavelength_resolution, align = 'center')
plt.title("Convoluted Fabry-Perot + realistic filter against discretized points")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Transmision")
plt.grid()
plt.show()
#plt.close()

del wavelength_resolution
#%%Gaussian fitting
pure_disc_peaks = sci.find_peaks(pure_discretized_convolution)[0]
pure_disc_valleys = sci.find_peaks(-pure_discretized_convolution)[0]

disc_peaks = sci.find_peaks(discretized_convolution)[0]
disc_valleys = sci.find_peaks(-discretized_convolution)[0]

acc_disc_peaks = sci.find_peaks(acc_discretized_convolution)[0]
acc_disc_valleys = sci.find_peaks(-acc_discretized_convolution)[0]


#simple visual inspection to see if satisfactory
dummy = [1 for i in range(len(perfect_lambda))]
plt.figure()
#plt.plot(1e9*whitelight, convolution)
plt.plot(1e9*whitelight, accurate_convolution)
plt.plot(1e9*perfect_lambda, dummy, 'yo')
# plt.plot(1e9*edges, discretized_convolution, 'ro')
# plt.plot(1e9*edges[disc_peaks], discretized_convolution[disc_peaks], 'bo')
# plt.plot(1e9*edges[disc_valleys], discretized_convolution[disc_valleys], 'go')
plt.plot(1e9*acc_edges, acc_discretized_convolution, 'ro')
plt.plot(1e9*acc_edges[acc_disc_peaks], acc_discretized_convolution[acc_disc_peaks], 'bo')
plt.plot(1e9*acc_edges[acc_disc_valleys], acc_discretized_convolution[acc_disc_valleys], 'go')


        
pure_gmeans, pure_gstd = opsys.gauss_model(pure_disc_peaks, pure_disc_valleys, pure_discretized_convolution, \
                                     pure_edges, vtol=0.05, show = True, perfect = perfect_lambda[1:-1])
pure_gmeans_w, pure_gstd_w = opsys.gauss_model(pure_disc_peaks, pure_disc_valleys, pure_discretized_convolution, \
                                     pure_edges, vtol=0.05, show = True, perfect = perfect_lambda[1:-1], weights=1)
   
gmeans, gstd = opsys.gauss_model(disc_peaks, disc_valleys, discretized_convolution, \
                          edges, vtol = 0.05, show = True, perfect = perfect_lambda[1:-1])
gmeans_w, gstd_w = opsys.gauss_model(disc_peaks, disc_valleys, discretized_convolution, \
                          edges, vtol = 0.05, show = True, perfect = perfect_lambda[1:-1], weights=1)

acc_gmean, acc_gstd = opsys.gauss_model(acc_disc_peaks, acc_disc_valleys, acc_discretized_convolution,\
                                  acc_edges, vtol = 0.05, show = True, perfect = perfect_lambda[1:-1])
acc_gmean_w, acc_gstd_w = opsys.gauss_model(acc_disc_peaks, acc_disc_valleys, acc_discretized_convolution,\
                                  acc_edges, vtol = 0.05, show = True, perfect = perfect_lambda[1:-1], weights=1)
del dummy
#%%ERF model
#@custom_model




#acc_gmean_erf = erf_model(acc_disc_peaks, acc_disc_valleys, acc_discretized_convolution, acc_edges, vtol = 0.2, show = True, perfect = perfect_lambda[1:-1])

pure_gmeans_erf= opsys.erf_model(pure_disc_peaks, pure_disc_valleys, \
                                 pure_discretized_convolution, pure_edges,\
                                 pure_gmeans, pure_gstd, vtol = 0.3, \
                                 show = True, incr = increment)
    
pure_gmeans_erf_w = opsys.erf_model(pure_disc_peaks, pure_disc_valleys, \
                                    pure_discretized_convolution, pure_edges,\
                                    pure_gmeans_w, pure_gstd_w, vtol = 0.3, \
                                    show = True, weights=1, incr = increment)
    

acc_gmean_erf_w = opsys.erf_model(acc_disc_peaks, acc_disc_valleys, \
                                  acc_discretized_convolution, acc_edges,\
                                  acc_gmean_w, acc_gstd_w, vtol = 0.3, \
                                  show = True, weights=1, incr = increment)

    
                                      
# gmean_erf= erf_model(disc_peaks, disc_valleys, discretized_convolution, edges,\
#                       gmeans, gstd, vtol = 0.3, show = True, perfect = perfect_lambda[1:-1])
    
# gmean_erf_w = erf_model(disc_peaks, disc_valleys, discretized_convolution, edges,\
#                      gmeans_w, gstd_w, vtol = 0.3, show = True, weights = True, perfect = perfect_lambda[1:-1])
#%%Gaussian estimation of peaks vs fabry-perot
c = 3e8

def peak_allignment(means, reference):
    global c
    aligned_peaks = means - reference
    aligned_peaks /= reference
    aligned_peaks *= c
    
    return aligned_peaks

#def peak_permutator


#For gaussian modelling
pure_g_aligned_peaks = peak_allignment(pure_gmeans, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_analytic = peak_allignment(pure_gmeans, perfect_lambda[1:-1])
pure_g_aligned_peaks_w = peak_allignment(pure_gmeans_w, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_w_analytic = peak_allignment(pure_gmeans_w, perfect_lambda[1:-1])

plt.figure()
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_peaks, label = 'Python peaks, unweighted')
plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_analytic, label = 'Analytic peaks, unweighted')
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_peaks_w, label = 'Python peaks, weighted')
plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_w_analytic, label = 'Analytic peaks, weighted')
plt.legend()
plt.title("Error between Fabry-Perot and discretized Fabry-Perot, guaussian fit")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.grid()
plt.show()

acc_g_aligned_peaks = peak_allignment(acc_gmean, whitelight[fb_peaks[1:-1]])
acc_g_aligned_peaks_analytic = peak_allignment(acc_gmean, perfect_lambda[1:-1])
acc_g_aligned_peaks_w = peak_allignment(acc_gmean_w, whitelight[fb_peaks[1:-1]])
#TODO
#FABRY-PEROT+INFRA GAUSSIAN
acc_g_aligned_peaks_w_analytic = peak_allignment(acc_gmean_w, perfect_lambda[1:-1])

plt.figure()
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_peaks-acc_g_aligned_peaks.mean(), label = 'Python peaks, unweighted')
plt.plot(1e9*perfect_lambda[1:-1], acc_g_aligned_peaks_analytic, label = 'Analytic peaks, unweighted')
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_peaks_w-acc_g_aligned_peaks_w.mean(), label = 'Python peaks, weighted')
plt.plot(1e9*perfect_lambda[1:-1], acc_g_aligned_peaks_w_analytic, label = 'Analytic peaks, weighted')
plt.legend()
plt.title("Error between Fabry-Perot and discretized Fabry-Perot + Realistic filter, guaussian fit")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.grid()
plt.show()

#for ERF modeling
pure_g_aligned_erf_peaks = peak_allignment(pure_gmeans_erf, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_erf_analytic = peak_allignment(pure_gmeans_erf, perfect_lambda[1:-1])
pure_g_aligned_peaks_erf_w = peak_allignment(pure_gmeans_erf_w, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_erf_w_analytic = peak_allignment(pure_gmeans_erf_w, perfect_lambda[1:-1])

plt.figure()
plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_erf_peaks - pure_g_aligned_erf_peaks.mean(), label = 'Python peaks, unweighted')
#plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_erf_analytic, label = 'Analytic peaks, unweighted')
plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_peaks_erf_w - pure_g_aligned_peaks_erf_w.mean(), label = 'Python peaks, weighted')
#plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_erf_w_analytic, label = 'Analytic peaks, weighted')
plt.legend()
plt.title("Error between Fabry-Perot and discretized Fabry-Perot, erf fit")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.grid()
plt.show()

#TODO
#FABRY-PEROT+INFRA ERF
acc_g_aligned_erf_peaks_w_analytic = peak_allignment(acc_gmean_erf_w, perfect_lambda[1:-1])

# g_alligned_peaks = gmeans - whitelight[fb_peaks[0][1:-1]]
# g_alligned_peaks /= whitelight[fb_peaks[0][1:-1]]
# g_alligned_peaks *= c

# p_alligned_peaks = gmeans - perfect_lambda[1:-1]
# p_alligned_peaks /= perfect_lambda[1:-1]
# p_alligned_peaks *= c

# acc_g_aligned_peaks = acc_gmean - whitelight[fb_peaks[1:-1]]
# acc_g_aligned_peaks /= whitelight[fb_peaks[1:-1]]
# acc_g_aligned_peaks *= c

# plt.figure()
# # plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], g_alligned_peaks, label = 'Python peaks(Ideal filter)')
# #plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], p_alligned_peaks, label = 'Analytic peaks(Ideal filter)')
# plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_peaks, label = 'Realistic filter')
# plt.legend()
# plt.title("Error between fb and discretized conv, guaussian fit")
# plt.xlabel("Wavelength [nm]")
# plt.ylabel("Error in speed[m/s]")
# plt.grid()
# plt.show()

# acc_g_aligned_erf_peaks = acc_gmean_erf - whitelight[fb_peaks[0][1:-1]]
# acc_g_aligned_erf_peaks /= whitelight[fb_peaks[0][1:-1]]
# acc_g_aligned_erf_peaks *= c

# plt.figure()
# plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], acc_g_aligned_erf_peaks, label = 'Realistic filter')
# plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], acc_g_aligned_erf_peaks - acc_g_aligned_erf_peaks.mean(), label = 'Realistic filter, de-meaned')
# plt.legend()
# plt.title("Error between fb and discretized conv, erf fit")
# plt.xlabel("Wavelength [nm]")
# plt.ylabel("Error in speed[m/s]")
# plt.grid()
# plt.show()

# g_aligned_erf_peaks = gmean_erf - whitelight[fb_peaks[0][1:-1]]
# g_aligned_erf_peaks /= whitelight[fb_peaks[0][1:-1]]
# g_aligned_erf_peaks *= c

# plt.figure()
# plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], g_aligned_erf_peaks, label = 'Realistic filter')
# plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], g_aligned_erf_peaks - g_aligned_erf_peaks.mean(), label = 'Realistic filter, de-meaned')
# plt.legend()
# plt.title("Error between fb and discretized conv, erf fit")
# plt.xlabel("Wavelength [nm]")
# plt.ylabel("Error in speed[m/s]")
# plt.grid()
# plt.show()

#del c
#%%
#plt.close('all')

#%%plots for thesis

plt.figure()
plt.plot(1e9*perfect_lambda[1:-1], acc_aligned_peaks_analytic[1:-1], label = 'No Discretisation')
plt.plot(1e9*perfect_lambda[1:-1], acc_g_aligned_peaks_w_analytic, label = 'Gaussian Estimation')
plt.plot(1e9*perfect_lambda[1:-1], acc_g_aligned_erf_peaks_w_analytic, label = 'ERF Estimation')
plt.legend()
plt.title("Radial Speed of Fabry-Perot+Neutral Density Filter, weight = x^2, tilt = 2")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Radial Speed[m/s]")
plt.grid()
plt.show()

np.mean(acc_aligned_peaks_analytic[1:-1])
np.var(acc_aligned_peaks_analytic[1:-1])

np.mean(acc_g_aligned_peaks_w_analytic)
np.var(acc_g_aligned_peaks_w_analytic)
       
np.mean(acc_g_aligned_erf_peaks_w_analytic)
np.var(acc_g_aligned_erf_peaks_w_analytic)
#%%Error of errors

# diff_models = (g_alligned_peaks - acc_g_aligned_peaks)#/g_alligned_peaks

# plt.figure()
# plt.plot(1e9*whitelight[fb_peaks[0][1:-1]], diff_models)
# plt.title("Error of errors")
# plt.xlabel("Wavelength [nm]")
# plt.ylabel("Error in m/s")
# plt.grid()
# plt.show()
my_weight = '1/sqrt(x)'

acc_g_aligned_peaks = peak_allignment(acc_gmean, whitelight[fb_peaks[1:-1]])
acc_g_aligned_peaks_analytic = peak_allignment(acc_gmean, perfect_lambda[1:-1])
acc_g_aligned_peaks_w = peak_allignment(acc_gmean_w, whitelight[fb_peaks[1:-1]])
acc_g_aligned_peaks_w_analytic = peak_allignment(acc_gmean_w, perfect_lambda[1:-1])
plt.figure()
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_peaks-acc_g_aligned_peaks.mean(), label = 'Python peaks, unweighted')
#plt.plot(1e9*perfect_lambda[1:-1], acc_g_aligned_peaks_analytic, label = 'Analytic peaks, unweighted')
plt.plot(1e9*whitelight[acc_conv_peaks[1:-1]], acc_aligned_peaks[1:-1], label = 'Not discretized')#Python peaks, 
plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_peaks_w-acc_g_aligned_peaks_w.mean(), label = 'Weight = '+my_weight)#Python peaks, 
plt.title("")
plt.legend()
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.title("Error between Fabry-Perot and discretized Fabry-Perot + realistic filter, guaussian fit")
plt.grid()
plt.show()

acc_g_aligned_erf_peaks_analytic = peak_allignment(acc_gmean_erf_w, perfect_lambda[1:-1])
plt.figure()
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_peaks-acc_g_aligned_peaks.mean(), label = 'Python peaks, unweighted')
#plt.plot(1e9*perfect_lambda[1:-1], acc_g_aligned_peaks_analytic, label = 'Analytic peaks, unweighted')
plt.plot(1e9*whitelight[acc_conv_peaks[1:-1]], acc_aligned_peaks[1:-1], label = 'Not discretized')#Python peaks, 
plt.plot(1e9*whitelight[fb_peaks[1:-1]], acc_g_aligned_erf_peaks_w_analytic-acc_g_aligned_erf_peaks.mean(), label = 'Weight = '+my_weight)#Python peaks, 
plt.title("")
plt.legend()
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.title("Error between Fabry-Perot and discretized Fabry-Perot + realistic filter, erf fit")
plt.grid()
plt.show()


pure_g_aligned_peaks = peak_allignment(pure_gmeans, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_analytic = peak_allignment(pure_gmeans, perfect_lambda[1:-1])
pure_g_aligned_peaks_w = peak_allignment(pure_gmeans_w, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_w_analytic = peak_allignment(pure_gmeans_w, perfect_lambda[1:-1])
plt.figure()
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_peaks, label = 'Python peaks, unweighted')
plt.plot(1e9*perfect_lambda[1:-1], pure_aligned_peaks[1:-1], label = 'Not discretized')
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_peaks_w, label = 'Python peaks, weighted')
plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_w_analytic, label = 'Weight = '+my_weight)
plt.legend()
plt.title("Error between Fabry-Perot and discretized Fabry-Perot, guaussian fit")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.grid()
plt.show()

pure_g_aligned_erf_peaks = peak_allignment(pure_gmeans_erf, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_erf_analytic = peak_allignment(pure_gmeans_erf, perfect_lambda[1:-1])
pure_g_aligned_peaks_erf_w = peak_allignment(pure_gmeans_erf_w, whitelight[fb_peaks[1:-1]])
pure_g_aligned_peaks_erf_w_analytic = peak_allignment(pure_gmeans_erf_w, perfect_lambda[1:-1])

plt.figure()
#plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_erf_peaks - pure_g_aligned_erf_peaks.mean(), label = 'Python peaks, unweighted')
#plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_erf_analytic, label = 'Analytic peaks, unweighted')
plt.plot(1e9*perfect_lambda[1:-1], pure_aligned_peaks[1:-1], label = 'Not discretized')
plt.plot(1e9*whitelight[fb_peaks[1:-1]], pure_g_aligned_peaks_erf_w - pure_g_aligned_peaks_erf_w.mean(), label = 'Weight = '+my_weight)

#plt.plot(1e9*perfect_lambda[1:-1], pure_g_aligned_peaks_erf_w_analytic, label = 'Analytic peaks, weighted')
plt.legend()
plt.title("Error between Fabry-Perot and discretized Fabry-Perot, erf fit")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Error in speed[m/s]")
plt.grid()
plt.show()