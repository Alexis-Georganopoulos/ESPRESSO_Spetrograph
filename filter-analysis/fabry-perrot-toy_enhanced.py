# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 02:05:23 2021

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sci
import scipy.stats as stat
import scipy.integrate as intgr
from astropy.modeling import models, fitting
import time
from Espresso import Perfect_Fabry_Perot, Simple_Filter, Realistic_Filter, Optical_System, white_light_generator

virtual_steps = 100000 #just how many increments we want for some range of numbers
lambda_target = 600e-9 #[m]
lambda_deviation = 0.3e-9 #[m]this varies depending on application    0.3e-9
lambda_min = lambda_target - lambda_deviation
lambda_max = lambda_target + lambda_deviation
increment = (lambda_max-lambda_min)/virtual_steps

#general parameters of the filters & fabry_perot(the rest are default arguments)
etalon_spacing = 7.6e-3
filter_thickness = 2e-3
tilt = 2 *np.pi/180

#reflectance of the netrual density filter coating
R_nd = 0.05#abs((1-(2.2714+3.8077j))/(1+(2.2714+3.8077j)))**2
path_nd = 43.718e-3
semi_nd = 0.993e-3

#for generating the guassian needed for convolution
sample_space = 130000#also temporary, depending on what we want
cuttoff = 0.001 #the percentage height at which ignore the rest of the gaussian

#for discretization
resolution = lambda_target/130000/3

#for beat pattern examination
c = 3e8

whitelight = white_light_generator(lambda_min, lambda_max, virtual_steps)

fabry_perot = Perfect_Fabry_Perot(whitelight, etalon_spacing)
fabry_perot.plot_transmitance()
simple_filter = Simple_Filter(whitelight,filter_thickness)
#simple_filter.plot_transmitance()
real_filter = Realistic_Filter(whitelight, filter_thickness)
#real_filter.plot_transmitance()

optical_system = Optical_System(whitelight, fabry_perot, real_filter)
optical_system.plot_transmitance()
optical_system.generate_gaussian(increment, virtual_steps, lambda_target, sample_space, cuttoff, show = True)
optical_system.convolve_system(show = True)
optical_system.beat_pattern_pre_discretization(fabry_perot, c)
optical_system.discretize_system(resolution, increment, virtual_steps, show = True)
optical_system.gauss_peaks(show = True)
optical_system.beat_pattern(fabry_perot, c, show = True)

# optical_system_pure = Optical_System(whitelight, fabry_perot)
# optical_system_pure.plot_transmitance()
# optical_system_pure.generate_gaussian(increment, virtual_steps, lambda_target, sample_space, cuttoff, show = True)
# optical_system_pure.convolve_system(show = True)
# optical_system_pure.beat_pattern_pre_discretization(fabry_perot, c)
# optical_system_pure.discretize_system(resolution, increment, virtual_steps, show = True)
# optical_system_pure.gauss_peaks(show = True)
# optical_system_pure.beat_pattern(fabry_perot, c, show = True)

# we assume the reflectance does not vary with incident angle
# neutral_density = Realistic_Filter(whitelight,filter_thickness, R1 = R_nd, \
#                                     path_length = path_nd, semi_diameter = semi_nd)
# neutral_density.plot_transmitance()
    
# optical_system_extended = Optical_System(whitelight,  fabry_perot,real_filter, neutral_density)
# optical_system_extended.plot_transmitance()
# optical_system_extended.generate_gaussian(increment, virtual_steps, lambda_target, sample_space, cuttoff, show = True)
# optical_system_extended.convolve_system(show = True)
# optical_system_extended.discretize_system(resolution, increment, virtual_steps, show = True)
# optical_system_extended.gauss_peaks(show = True)
# optical_system_extended.beat_pattern(fabry_perot, c, show = True)

# optical_system_nd = Optical_System(whitelight, fabry_perot,neutral_density)
# optical_system_nd.plot_transmitance()
# optical_system_nd.generate_gaussian(increment, virtual_steps, lambda_target, sample_space, cuttoff, show = True)
# optical_system_nd.convolve_system(show = True)
# optical_system_nd.discretize_system(resolution, increment, virtual_steps, show = True)
# optical_system_nd.gauss_peaks(show = True)
# optical_system_nd.beat_pattern(fabry_perot, c, show = True)
    
