# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 16:02:27 2021

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat
from astropy.modeling import models, fitting
import scipy.signal as sci

def white_light_generator(lambda_min, lambda_max, lambda_divisions):
    difference = lambda_max - lambda_min
    increment = difference/lambda_divisions
    return np.arange(lambda_min, lambda_max+ increment, increment)

class Fabry_Perot:
    
    def __init__(self, theta, wavelength, n, etalon_spacing, R1, R2):
        self.theta = np.asmatrix(theta)
        self.wavelength = np.asmatrix(wavelength)
        self.n = n
        self.etalon_spacing= etalon_spacing
        self.R1 = R1
        self.R2 = R2
        self._transmitance = False
        self._peaks = False
        #resulting matrix is:
        #------->theta increasing
        #|
        #|wavelength increasing
        #\/
    def _sin_delta_square(self):
            scalar_delta = 2*np.pi*2*self.n*self.etalon_spacing
            matrix_delta = ((1/self.wavelength).transpose()) @ np.cos(self.theta)
            delta = scalar_delta*matrix_delta
            s_delta = np.sin(delta/2)
            s_delta_2 = np.multiply(s_delta,s_delta)
            return s_delta_2
        
    def _return_form(self,result):
        if result.shape[1] == 1:
            self._transmitance = np.squeeze(np.asarray(result))
            return self._transmitance
        else:
            self._transmitance = result
            return self._transmitance
        
    #cannot handle matrix transmitance    
    def peaks(self):
        if type(self._transmitance) == bool:
            raise Warning('Cannot find peaks before transmitance calculation')
            return
        if type(self._peaks) == bool:
            self._peaks = sci.find_peaks(self._transmitance)[0]
        return self._peaks
    
    def transmitance(self):
        if type(self._transmitance) == bool:
            geomean = np.sqrt(self.R1*self.R2)
            sin_delta_2 = self._sin_delta_square()
            numerator = (1 - self.R1)*(1 - self.R2)
            denominator = (1 - geomean)**2 + 4*geomean*(sin_delta_2)
            result = numerator/denominator #The transmitance
            return self._return_form(result)
        else:
            return self._transmitance
        
    def plot_transmitance(self):
        y = self.transmitance()
        x = np.squeeze(np.asarray(self.wavelength))
        plt.figure()
        plt.plot(1e9*x,y)
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Transmitance")
        plt.grid()
        plt.show()
    
class Perfect_Fabry_Perot(Fabry_Perot):
    def __init__(self, wavelength, etalon_spacing, n = 1, theta = 0, F = 12.27):
        # b = (4+2*F)/F
        # R = (b-np.sqrt(b**2-4))/2
        super().__init__(theta, wavelength, n, etalon_spacing, 0, 0)
        self.F = F
        self._transmitance = False
        
    def transmitance(self):
        if type(self._transmitance) == bool:
            sin_delta_2 = self._sin_delta_square()
            product = self.F*sin_delta_2
            result = 1/(1 + product)
            return self._return_form(result)
        else:
            return self._transmitance

class Simple_Filter(Fabry_Perot):
    def __init__(self, wavelength,filter_thickness, theta = 0, n = 1.5, \
                 percentage = 0.05, peak = 1):
        
        super().__init__(theta, wavelength, n, filter_thickness, 0, 0)
        self.percentage = percentage
        self.peak = peak
        self._transmitance = False
         
    def transmitance(self):
        if type(self._transmitance) == bool:
            sin_delta_2 = self._sin_delta_square()
            amplitude = self.peak*self.percentage
            result = amplitude*sin_delta_2+1-amplitude
            return self._return_form(result)
        else:
            return self._transmitance
        
class Realistic_Filter(Fabry_Perot):
    ##reflecatance of generic glass is around 0.04
    def __init__(self, wavelength, filter_thickness, n = 1.5, tilt = 0, R1 = 0.0625, \
                 R2 = 0.04, path_length = 56.547e-3, semi_diameter = 1.129e-3, divisions = 100):
        theta_upper = np.arctan(semi_diameter/path_length)
        theta_lower = -theta_upper
        theta_upper, theta_lower = (tilt + x for x in (theta_upper, theta_lower))
        theta_increment = (theta_upper - theta_lower)/divisions
        theta_integral_vector = np.arange(theta_lower, theta_upper+theta_increment, theta_increment)
        super().__init__(theta_integral_vector, wavelength, n, filter_thickness, R1, R2)
        self._transmitance = False
        
    def transmitance(self):
        if type(self._transmitance) == bool:
            
            integration_matrix = super().transmitance()
            #return integration_matrix
            I = (integration_matrix[:,0] + integration_matrix[:,-1])/2
            I = integration_matrix.sum(axis = 1) - I
            I /= I.max() 
            return self._return_form(I)
        else:
            return self._transmitance
    
class Optical_System:
    
    def __init__(self, wavelength, *args):
        if type(args) != tuple:
            self.objects = [args]
        elif len(args) >= 0:
            temp = []
            for i in args:
                temp.append(i)
            self.objects = temp
        else:
            raise Warning('unkown error in Optical_System: __init__')
        self.wavelength = wavelength
        self._reset_hidden()
        
    def _reset_hidden(self):
        self._transmitance = False
        self._gaussian = False
        self._gaussian_x = False 
        self._convolution = False
        self._discretized = False
        self._edges = False
        self._disc_res = False
        self._means = False
        self._i_computational_peaks = False
        self._i_computational_valeys = False
        self._gauss_beats = False
        
    def append(self,*args):
        if type(args) != tuple:
            self.objects.append(args)
        elif len(args) >= 0:
            for i in args:
                self.objects.append(i)
        else:
            raise Warning('unkown error in Optical_System: append')
        self._reset_hidden()
            
    def transmitance(self):
        if type(self._transmitance) == bool:
            result = 1
            if len(self.objects) == 0:
                raise Warning('Empty system')
                return
            else:
                for i in self.objects:
                   result *= i.transmitance()
                
                self._transmitance = result
                return self._transmitance
        else:
            return self._transmitance
        
    def generate_gaussian(self, increment, virtual_steps, lambda_target, sample_space, cuttoff, show = False):
        if type(self._gaussian) == bool:
            gaussian_resolution = lambda_target/sample_space #the FWHM we want
            sigma_guassian = gaussian_resolution/(2*np.sqrt(2*np.log(2)))
        
            xrange = np.arange(-int(virtual_steps/2)*increment, \
                               (int(virtual_steps/2)+1)*increment, increment)
                
            exp_gauss = np.exp(-(xrange**2)/(2*sigma_guassian**2))
            indexes = exp_gauss > cuttoff*np.max(exp_gauss)
            exp_gauss = exp_gauss[indexes]
            xrange = xrange[indexes]
            self._gaussian = exp_gauss
            self._gaussian_x = xrange
        else:
            print('Gaussian already generated')
        if show:
            plt.figure()
            plt.plot(1e9*self._gaussian_x, self._gaussian)
            plt.title("Unit Gaussian(Not normalised)")
            plt.xlabel("Centered wavelength [nm]")
            plt.ylabel("")
            plt.grid()
            plt.show()
            
     # def external_gaussian(self, )
            
    def convolve_system(self, show = False):
        if type(self._gaussian) == bool:
            raise Warning('No gaussian to convolute with')
            return
        if type(self._convolution) == bool:
            convolution = np.convolve(self._gaussian, self.transmitance(), mode = 'same')
            convolution /= convolution.max()
            self._convolution = convolution
        else:
            print('Convolution already completed')
        if show:
            plt.figure()
            plt.plot(1e9*self.wavelength,self._convolution)
            plt.title("Fabry-perot + filter after gaussian convolution")
            plt.xlabel("Wavelength [nm]")
            plt.ylabel("Transmitance")
            plt.grid()
            plt.show()
            
    ##############################################################questionable function
    def beat_pattern_pre_discretization(self,comparator, c):
        if type(self._convolution) == bool:
            raise Warning('No convolution product to compare to')
            return
        idx_peaks = sci.find_peaks(self._convolution)[0]  
        beats = self.wavelength[idx_peaks[1:-1]] - self.wavelength[comparator.peaks()[1:-1]]
        beats /= self.wavelength[comparator.peaks()[1:-1]]
        beats *= c
        
        plt.figure()
        plt.plot(1e9*self.wavelength[comparator.peaks()[1:-1]],beats)
        plt.title('Beat pattern pre-discretization')
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Error in speed [m/s]")
        plt.grid()
        plt.show()
    ###############################################################
            
    def discretize_system(self,resolution, increment, virtual_steps, show = False):
        if type(self._convolution) == bool:
            raise Warning('No convolution to discretize')
            return
        if type(self._discretized) == bool:
            statistic, edges, _ = stat.binned_statistic(self.wavelength, self._convolution, \
                                bins = int(increment*virtual_steps/resolution))
    
            for i in range(len(edges)-1):
                edges[i] = (edges[i]+edges[i+1])/2
            edges = np.delete(edges, -1)
            self._discretized = statistic
            self._edges = edges
            self._disc_res = resolution
        else:
            print('System already discretized')
        if show:
            plt.figure()
            plt.plot(1e9*self.wavelength,self._convolution)
            plt.plot(1e9*self._edges, self._discretized, 'ro')
            plt.bar(1e9*self._edges,self._discretized, alpha = 0.5, \
                    color ='orange', width = 1e9*self._disc_res, align = 'center')
            plt.title("Convoluted Fabry-Perot against discretized points")
            plt.xlabel("Wavelength [nm]")
            plt.ylabel("Transmitance")
            plt.grid()
            plt.show()
    
    #plotting functionality is inneficient/inelegant(have to rerun loop)
    def gauss_peaks(self, vtol = 0.2, show = False):
        if type(self._discretized) == bool:
            raise Warning('No discretization to fit')
            return
        if True:#type(self._means) == bool:
            self._i_computational_peaks = sci.find_peaks(self._discretized)[0]
            self._i_computational_valleys = sci.find_peaks(-self._discretized)[0]
            #peaks = self._i_computational_peaks
            valleys = self._i_computational_valleys
            d_conv = self._discretized
            edg = self._edges
            mean_package = []
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
                g = fit_g(g_init, x, y)
                mean_package.append(g.mean.value)
                
                dummy_x = np.arange(x[0],x[-1], 1e-14)
                if show: 
                    plt.plot(x,y,'*', dummy_x,g(dummy_x), g.mean.value, g.amplitude.value, 'bo')
            plt.title("Gaussian peak estimations vs discretized points")
            plt.xlabel("Wavelength [nm]")
            plt.ylabel("Transmision")
            plt.grid()
            plt.show()

            self._means = mean_package
            
    def efr_peaks(self):
        if type(self._discretized) == bool:
            raise Warning('No discretization to fit')
            return
        pass
            
    #CAREFULL! .PEAKS() IS A METHOD IN THE FABRY-PEROT CLASS AND IS THEREFORE A TIGHT DEPENDANCE
    def beat_pattern(self,comparator, c, show = False):
        if type(self._means) == bool:
            raise Warning('No peak estimations to compare against')
            return
        if type(self._gauss_beats) == bool:
            beats = self._means - self.wavelength[comparator.peaks()[1:-1]]
            beats /= self.wavelength[comparator.peaks()[1:-1]]
            beats *= c
            self._gauss_beats = beats
        if show:
            plt.figure()
            plt.plot(1e9*self.wavelength[comparator.peaks()[1:-1]],self._gauss_beats)
            plt.xlabel("Wavelength [nm]")
            plt.ylabel("Error in speed [m/s]")
            plt.grid()
            plt.show() 
            
    def plot_transmitance(self):
        y = self.transmitance()
        x = self.wavelength
        plt.figure()
        plt.plot(1e9*x,y)
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Transmitance")
        plt.grid()
        plt.show()
        

        
        
