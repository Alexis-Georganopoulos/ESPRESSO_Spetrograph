 # -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:09:42 2022

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat
from astropy.modeling import models, fitting
from scipy.special import erf
from scipy.optimize import curve_fit
from tqdm import tqdm

###################### OPTICAL SYSYTEM FUNCTIONS ##############################


def white_light_generator(lambda_min, lambda_max, lambda_divisions):
    difference = lambda_max - lambda_min
    increment = difference/lambda_divisions
    return np.arange(lambda_min, lambda_max+ increment, increment)

##We assume for now that the finesse is constant for all lambda
## in fact finesse changes by the incedent angle for each wavelength(not reflectance)
def fabry_perot_transmitance(lambda_range, etalon_thickness, \
                             n = 1, theta = 0, F = 12.27):
    #F_coeff = 
    delta = (2*np.pi/lambda_range)*2*n*etalon_thickness*np.cos(theta)
    product = F*(np.sin(delta/2))**2
    return 1/(1 + product) ##The transmitance

##we model the transmitance directly as a sinusoid for a simple filter
def filter_transmitance(lambda_range,filter_thickness, \
                        theta = 0, n = 1, percentage = 0.05, peak = 1):

    delta = (2*np.pi/lambda_range)*2*n*filter_thickness*np.cos(theta)
    amplitude = peak*percentage
    wave = (np.sin(delta/2))**2
    return amplitude*wave+1-amplitude

#realistic model of an infrared fiter
# From data sheet of 800FL07-12.5
#Transmitance @600nm ~ 0.9
#reflectance of glass ~ 0.04(= R2) -> T_glass ~ 0.96
#-> T_filter = 0.9 = T_coating * T_glass
#-> T_coating ~ 0.9/0.96 = 0.9375 -> R_caoting(= R1) = 1 - 0.9375 = 0.0625
def accurate_filter_transmitance(lambda_range,filter_thickness, \
                                 path_length = 56.547e-3, semi_diameter = 1.129e-3, \
                                 n = 1, R1 = 0.0625, R2 = 0.04, tilt_deg = 0):
    
    def _mini_fabry_perot(theta, one_lambda):
        delta = (2*np.pi/one_lambda)*2*n*filter_thickness*np.cos(theta)
        geomean = np.sqrt(R1*R2)
        
        numerator = (1 - R1)*(1 - R2)
        denominator = (1-geomean)**2 + 4*geomean*(np.sin(delta/2)**2)
        transmitance = numerator/denominator
        return transmitance
    
    def _deg2rad (deg):
        return deg*2 *np.pi/180
    
    transmitance_array = []
    theta_upper = np.arctan(semi_diameter/path_length)
    theta_lower = -theta_upper
    theta_upper, theta_lower = (_deg2rad(tilt_deg) + x for x in (theta_upper, theta_lower))
    #####
    divisions = 100
    theta_increment = (theta_upper - theta_lower)/divisions
    theta_integral_vector = np.arange(theta_lower, theta_upper+theta_increment, theta_increment)
    for wavelength in tqdm(lambda_range):
        #####Trapezoidal Integration
        I = _mini_fabry_perot(theta_integral_vector, wavelength)
        I = (np.sum(I) - (I[0]+I[-1])/2)*theta_increment
        transmitance_array.append(I)

    transmitance_array /= np.max(transmitance_array)
    return transmitance_array

# realistic model of a neutral density filter
def neutral_density(lambda_range, R1, filter_thickness = 2e-3,\
                    n = 1.5, R2 = 0.04, theta = 0):
    
    delta = (2*np.pi/lambda_range)*2*n*filter_thickness*np.cos(theta)
    geomean = np.sqrt(R1*R2)
    
    numerator = (1 - R1)*(1 - R2)
    denominator = (1-geomean)**2 + 4*geomean*(np.sin(delta/2)**2)
    transmitance = numerator/denominator
    return transmitance

#for reference, the analytic peaks:
def analytic_fabry_perot_peaks(lower, upper, etalon, n = 1,theta = 0):
    k_high = int((2*n*etalon*np.cos(theta))/(upper))+1
    k = k_high - 1
    perfect_lambda = []
    while True:
        k = k + 1
        new_lamnda = (2*n*etalon*np.cos(theta))/(k)
        if new_lamnda < lower:
            perfect_lambda.reverse()
            perfect_lambda = np.asarray(perfect_lambda)
            return perfect_lambda
        
        perfect_lambda.append(new_lamnda)
        
###############################################################################
        
#################### UTILITY AND MODELING FUNCTIONS ###########################        
def gaussian(cuttoff, sample_space, target, upper, lower, steps):

    gaussian_resolution = target/sample_space #the FWHM we want
    sigma_guassian = gaussian_resolution/(2*np.sqrt(2*np.log(2)))

    increment = (upper-lower)/steps
    xrange = np.arange(-int(steps/2)*increment, \
                       (int(steps/2)+1)*increment, increment)
    #mag_gauss = (1/(np.sqrt(2*np.pi)*sigma_guassian))
    exp_gauss = np.exp(-(xrange**2)/(2*sigma_guassian**2))
    indexes = exp_gauss > cuttoff*np.max(exp_gauss)
    exp_gauss = exp_gauss[indexes]
    #gauss = mag_gauss*exp_gauss
    xrange = xrange[indexes]
    # plt.figure()
    # plt.plot(xrange, exp_gauss)
    # plt.show()
    return exp_gauss, xrange
    #return gauss, xrange
    
def discretize(lambda_range, y_value, resolution, upper, lower):

    statistic, edges, _ = stat.binned_statistic(lambda_range, y_value, \
                                bins = int((upper-lower)/resolution))

    for i in range(len(edges)-1):
        edges[i] = (edges[i]+edges[i+1])/2
    edges = np.delete(edges, -1)
    
    return statistic, edges

def gauss_model(peaks, valleys, d_conv, edg, perfect, vtol = 0, show = True , weights = 0):
    #"show", "perfect" variables mainly used for debugging
    mean_package = []
    std_package = []
    if show: plt.figure()
    for i in range(len(valleys)-1):
        indices = np.arange(valleys[i],valleys[i+1]+1,1)
        indices = indices[d_conv[indices] > vtol*(d_conv[indices].max())]
        x = edg[indices]
        y = d_conv[indices]
        mu = (x[-1]-x[0])/2 + x[0]
        sigma = abs(x[0]-x[-1])/10 #this denominator is arbitrary
        
        #carefull how you initialise! it is sensitive
        g_init = models.Gaussian1D(amplitude=1., mean=mu, stddev=sigma)
        
        fit_g = fitting.LevMarLSQFitter()

        if weights == 0:
            g = fit_g(g_init, x, y)
        elif weights == 1:
            g = fit_g(g_init, x, y, weights=np.array(y)**2)
        elif weights == 2:
            g = fit_g(g_init, x, y, weights=1/np.sqrt(np.array(y)))
        else:
            raise Warning("invalid weight in gauss_model")
            
        dummy_x = np.arange(x[0],x[-1], 1e-14)
        if show: 
            #dummy_perfect = [1 for i in range(len(perfect))]
            plt.plot(x,y,'*', dummy_x,g(dummy_x), g.mean.value, g.amplitude.value, 'bo')# perfect, dummy_perfect,'ro',
        mean_package.append(g.mean.value)
        std_package.append(g.stddev.value)
    if show:
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Transmission")
        plt.grid()
        plt.show()
    return mean_package, std_package

def erf_model(peaks, valleys, d_conv, edg, means, stds, incr, vtol = 0, show = True, weights = False):
    
    mean_package = []
    
    def _erfunc(x, mFL =0, a=0, b=1,c=0):
        return mFL*erf((x-a)/(b*np.sqrt(2))) + c
    
    if show: plt.figure()
    for i in range(len(valleys)-1):
        indices = np.arange(valleys[i],valleys[i+1]+1,1)
        indices = indices[d_conv[indices] > vtol*(d_conv[indices].max())]
        x = edg[indices]
        y = d_conv[indices]
        #x_adj = x + incr/2 
        my_erf = [0]
        for j in range(len(y)-1):
            temp = incr*(y[j]+y[j+1])/2 + my_erf[j]
            my_erf.append(temp)
        del my_erf[0]
        x = np.delete(x,0)
        
        # print(len(my_erf))
        # print(len(x))
        
        if weights == 0:
            my_sig = np.array(y[1:])
        elif weights == 1:
            my_sig = np.array(y[1:])**2
        elif weights == 2:
            my_sig = (1/np.sqrt(np.array(y[1:])))
        else:
            raise Warning("invalid weight in erf_model")
            
        params, extras = curve_fit(_erfunc, x, my_erf, \
                                   p0 = [1e-16,means[i],stds[i], 0], \
                                   sigma = my_sig, method='lm')#std = 0.4e-6
        plt.plot(x,my_erf,'*',x,_erfunc(x, *params))#params[1] are the means
        mean_package.append(params[1])

            # erf_init = _erfunc(mFL = 1e-16, a = means[i], b = stds[i], c = 0 )
            # fit_erf = fitting.LevMarLSQFitter()
            # erf = fit_erf(erf_init, x, my_erf)
            # plt.plot(x,my_erf,'*',x,erf(x))
            #mean_package.append(erf.)

    if show:
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Sampled erf(not normalised)")
        plt.grid()
        plt.show()
    return mean_package

###############################################################################

####################### USEFUL PLOTTING DEFINITIONS ###########################
