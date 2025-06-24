# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:40:06 2023

@author: vil034
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt



def get_overlap(ta, tb, qa, qb, shift, margin):
    ao = []
    bo = []
    
    aoq = []
    boq = []
    
    for i in range(len(tb)):
        b = tb[i]
        idx = find_nearest_idx(ta, b - shift)
        
        a = ta[idx]
                
        # c = np.abs(a - b + shift)
        if np.abs(a - b + shift) < margin:
            ao += [a]
            bo += [b]
            
            aoq += [qa[idx]]
            boq += [qb[i]]
    
    return np.array(ao), np.array(bo), np.array(aoq), np.array(boq)


def get_overlap_old(ta, tb, qa, qb, shift, margin):
    ao = []
    bo = []
    
    aoq = []
    boq = []
    for i in range(len(tb)):
        b = tb[i]
        
        m = np.where(np.abs(ta + shift - b) < 2*margin)
        
        
        # m1 = a + shift > tb - margin
        # m2 = a + shift < tb + margin
        # m = np.logical_and(m1, m2)

        if any(m):
            ao += list(ta[m])
            bo += [b]
            
            aoq += list(qa[m])
            boq += [qb[i]]
            
    return np.array(ao), np.array(bo), np.array(aoq), np.array(boq)

# Gaussian function for fitting 
def Gaussian2(t, amp, mean, std, floor):
    return amp*np.exp(-1/2*((t-mean)/std)**2) + floor

def Gaussian(t, amp, mean, std):
    return amp*np.exp(-1/2*((t-mean)/std)**2)





def fit(x, y, plot=False):
    my = np.max(y)
    param_init = [my, np.average(x), (np.max(x) - np.min(x))/10]
    # print(param_init)
    
    popt, pcov = curve_fit(Gaussian, x, y, p0=param_init, maxfev=2000)
  
    # perr = sum(np.sqrt(np.diag(pcov)))
    # print(popt)
    # # print(pcov)
    # print(perr)
    # cont = 0
    
    # while perr > 1 and cont < 100:
    #     param_init[1] += param_init[2]
    #     popt, pcov = curve_fit(Gaussian, x, y, p0=param_init, maxfev=1000)
    #     plt.plot(x, Gaussian(x ,popt[0],popt[1],popt[2], popt[3]) ,label='Fitted Gaussian',linestyle='dotted')
        
    #     perr = sum(np.sqrt(np.diag(pcov)))
    #     cont += 1
        
    # print('Fitting attempts: ', cont)
    
    if plot:
        plt.figure()
        plt.plot(x,y,color='black',label='Data',alpha=0.6)
        plt.plot(x, Gaussian(x, popt[0], popt[1], popt[2]) ,label='Fitted Gaussian',color='red',linestyle='dotted')
        plt.xlabel(r"$\tau (\mu s)$", fontsize=14)
        plt.ylabel(r"Counts",  fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.plot(x, Gaussian(x , param_init[0],  param_init[1], param_init[2]),label='Fitted Gaussian',color='green',linestyle='dotted')
        
    return popt[1], popt[2], popt[0]




def fit2(x, y, plot=False):
    my = np.max(y)
    param_init = [my, np.average(x), (np.max(x) - np.min(x))/10, 0]
    # print(param_init)
    
    popt, pcov = curve_fit(Gaussian2, x, y, p0=param_init, maxfev=2000)
  
    # perr = sum(np.sqrt(np.diag(pcov)))
    # print(popt)
    # # print(pcov)
    # print(perr)
    # cont = 0
    
    # while perr > 1 and cont < 100:
    #     param_init[1] += param_init[2]
    #     popt, pcov = curve_fit(Gaussian, x, y, p0=param_init, maxfev=1000)
    #     plt.plot(x, Gaussian(x ,popt[0],popt[1],popt[2], popt[3]) ,label='Fitted Gaussian',linestyle='dotted')
        
    #     perr = sum(np.sqrt(np.diag(pcov)))
    #     cont += 1
        
    # print('Fitting attempts: ', cont)
    
    if plot:
        plt.figure()
        plt.plot(x,y,color='black',label='Data',alpha=0.6)
        plt.plot(x, Gaussian2(x, popt[0], popt[1], popt[2], popt[3]) ,label='Fitted Gaussian',color='red',linestyle='dotted')
        plt.xlabel(r"$\tau (\mu s)$", fontsize=14)
        plt.ylabel(r"Counts",  fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.plot(x, Gaussian(x , param_init[0],  param_init[1], param_init[2]),label='Fitted Gaussian',color='green',linestyle='dotted')
        
    return popt[1], popt[2], popt[0], popt[3]

   

def fit_w_hist(data, bins):
    hist_data = np.histogram(data, bins)
    
    x = hist_data[1][:-1]
    y = hist_data[0]
    
    # plt.figure()
    # plt.plot(x,y)
    
    return fit2(x, y, True)

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    