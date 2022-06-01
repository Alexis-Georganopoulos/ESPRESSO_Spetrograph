# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:23:09 2021

@author: Alexis Philip George-Georganopoulos
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sci
import scipy.stats as stat
from astropy.modeling import models, fitting
import scipy.integrate as intgr

from PIL import Image, ImageFilter

import psf_analyzer as PSF
# a = np.asmatrix(np.arange(1,11,1))
# b = np.asmatrix(np.arange(10,110,10))
# increment = np.pi/1000
# testx = np.arange(0, 3*2*np.pi, increment)

# testy = np.sin(testx)

# plt.figure()
# plt.plot(testx,testy)

# def test_func(a):
#     return a[-1]#np.median(a)#
# bin_width = 100
# statistic, edges, _ = stat.binned_statistic(testx, testy, statistic = 'mean', bins = bin_width)

# for i in range(len(edges)-1):
#     edges[i] = (edges[i]+edges[i+1])/2

# edges = np.delete(edges, -1)

# plt.plot(edges, statistic, 'ro')

# plt.bar(edges,statistic, alpha = 0.5, color ='orange', width = -(testx[-1]-testx[0])/bin_width, align = 'center')
# x = np.arange(20,30+1,1)
# for i in list(range(len(x)-1)):
#     print(x[i])
#     print(x[i+1])
#     print("BREAK")    

# def some_func(x,y):
#     return np.sin(x) + np.cos(y)

# dummytest = intgr.quad(some_func, -np.pi/2, np.pi/2, args=(np.pi/2))

# #Error term for integration(we take the high-quality value as reference)
# 100*np.max(np.abs((np.array(transmitance_filter) - np.array(reference_transmition))/np.array(reference_transmition)))
# np.arctan(2.577e-3/294.475e-3)*180/np.pi

# First create some toy data:
# x = np.linspace(0, 2*np.pi, 400)
# y = np.sin(x**2)

# # Create just a figure and only one subplot
# fig, ax = plt.subplots()
# ax.plot(x, y)
# ax.set_title('Simple plot')

# # Create two subplots and unpack the output array immediately
# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
# ax1.plot(x, y)
# ax1.set_title('Sharing Y axis')
# ax2.scatter(x, y)

# # Create four polar axes and access them through the returned array
# fig, axs = plt.subplots(2, 2, subplot_kw=dict(projection="polar"))
# axs[0, 0].plot(x, y)
# axs[1, 1].scatter(x, y)

# # Share a X axis with each column of subplots
# plt.subplots(2, 2, sharex='col')

# # Share a Y axis with each row of subplots
# plt.subplots(2, 2, sharey='row')

# # Share both X and Y axes with all subplots
# plt.subplots(2, 2, sharex='all', sharey='all')

# # Note that this is the same as
# plt.subplots(2, 2, sharex=True, sharey=True)

# # Create figure number 10 with a single subplot
# # and clears it if it already exists.
# fig, ax = plt.subplots(num=10, clear=True)

# def A(m,n):
#     if m == 0:
#         return n+1
#     elif n == 0:
#         return A(m-1, 1)
#     else:
#         return A(m-1,B(m,n-1))
    
# def B(m,n):
#     if n == 0:
#         return m+1
#     elif m==0:
#         return B(1, n-1)
#     else:
#         return B(A(m-1,n),n-1)
    
# print(A(1,1))
# n_w = 9
# rootpath = 'C:/Users/User/Documents/EPFL/Robotics-2019_2020-Year 3/Masters Thesis/zemax_psf/'
# w_range = [i+1 for i in range(n_w)]
# lst = ['C1_W1.txt']
# lst = []

# for i in range(n_w):
#     lst.append(rootpath+'C1_W'+str(i+1)+'.txt')
    
# del n_w
# del i
# del w_range
# del rootpath

# plt.figure(figsize = (18,6))
# for i,f in enumerate(lst):
#     data = np.loadtxt(f,skiprows=17)
#     plt.subplot(1,8,1+i)
#     plt.imshow(data,extent = [-150,150,-150,150])

# test_const_1 = 5
# test_const_2 = 10
# test_var = 0
# for i in range(10):
#     test_var += i



