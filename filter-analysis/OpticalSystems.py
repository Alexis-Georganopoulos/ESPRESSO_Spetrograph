 # -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:09:42 2022

@author: User
"""
import numpy as np
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
    
    def deg2rad (deg):
        return deg*2 *np.pi/180
    
    transmitance_array = []
    theta_upper = np.arctan(semi_diameter/path_length)
    theta_lower = -theta_upper
    theta_upper, theta_lower = (deg2rad(tilt_deg) + x for x in (theta_upper, theta_lower))
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
def neutral_density(lambda_range, R1, filter_thickness = 2e-3, n = 1.5, R2 = 0.04, theta = 0):
    
    delta = (2*np.pi/lambda_range)*2*n*filter_thickness*np.cos(theta)
    geomean = np.sqrt(R1*R2)
    
    numerator = (1 - R1)*(1 - R2)
    denominator = (1-geomean)**2 + 4*geomean*(np.sin(delta/2)**2)
    transmitance = numerator/denominator
    return transmitance

###############################################################################

####################### USEFUL PLOTTING DEFINITIONS ###########################
