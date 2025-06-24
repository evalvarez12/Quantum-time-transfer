#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:26:10 2023

@author: eduardo
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt




# Gaussian function for fitting 
def Gaussian(t, amp, mean, std):
    return amp*np.exp(-1/2*((t-mean)/std)**2)



def fit(x, y, plot=False):

    param_init = [np.max(y), np.average(x), 0.01]
    
    my = np.mean(y)
    y2 = y - my

    popt, pocv = curve_fit(Gaussian, x, y2, p0=param_init)
    
    if plot:
        plt.plot(x,y,color='black',label='Data',alpha=0.6)
        plt.plot(x, Gaussian(x ,popt[0],popt[1],popt[2]) + my,label='Fitted Gaussian',color='red',linestyle='dotted')
        plt.xlabel(r"$\tau$")
        plt.ylabel(r"$c(\tau)$")
        
    return popt

   

def fit_w_hist(data, bins, plot=True):
    hist_data = np.histogram(data, bins)
    
    x = hist_data[1][:-1]
    y = hist_data[0]
    
    plt.plot(x,y)
    
    return fit(x, y, plot)
    
    