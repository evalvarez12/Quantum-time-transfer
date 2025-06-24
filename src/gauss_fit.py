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
def Gaussian(t, amp, mean, std, floor):
    return amp*np.exp(-1/2*((t-mean)/std)**2) + floor



def fit(x, y, plot=False):

    param_init = [np.max(y), np.average(x), (np.max(x) - np.min(x))/100, 0]
    print(param_init)
    
    
    popt, pcov = curve_fit(Gaussian, x, y, p0=param_init)
    
    perr = sum(np.sqrt(np.diag(pcov)))
    print(popt)
    # print(pcov)
    print(perr)
    cont = 0
    
    while perr > 1 and cont < 100:
        param_init[1] += param_init[2]
        popt, pcov = curve_fit(Gaussian, x, y, p0=param_init)
        
        perr = sum(np.sqrt(np.diag(pcov)))
        cont += 1
        
    print('Fitting attempts: ', cont)
    
    if plot:
        plt.figure()
        plt.plot(x,y,color='black',label='Data',alpha=0.6)
        plt.plot(x, Gaussian(x ,popt[0],popt[1],popt[2], popt[3]) ,label='Fitted Gaussian',color='red',linestyle='dotted')
        plt.xlabel(r"$\tau$")
        plt.ylabel(r"$c(\tau)$")
        
        # plt.plot(x, Gaussian(x ,param_init[0], param_init[1], param_init[2], param_init[3]),label='Fitted Gaussian',color='green',linestyle='dotted')
        
    return popt

   

def fit_w_hist(data, bins):
    hist_data = np.histogram(data, bins)
    
    x = hist_data[1][:-1]
    y = hist_data[0]
    
    # plt.figure()
    # plt.plot(x,y)
    
    return fit(x, y, True)
    
    