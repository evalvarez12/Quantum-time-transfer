# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:04:00 2023

@author: vil034
"""

import numpy as np
from scipy.fft import fft, ifft
from scipy.signal import correlate
import matplotlib.pyplot as plt

def to_analog(time_data, N, time_res):

    analog_data = np.zeros(int(N), dtype=bool)
        
    for ti in time_data:
         idx = np.floor(np.mod(ti/time_res, N))
         analog_data[int(idx)] = 1
    
    return analog_data



def cross_corr_mine(a, b):
    return ifft(np.conj(fft(a))*fft(b))


def analysis(ta, tb, n, total_time, plot=False, verbose=False):
    N = 2**n
    
    time_res = total_time/N 
    
    if verbose:
        print('Time resolution: ', time_res)
        print('Desired resolution: ', 10e-6)
    
    analog_alice = to_analog(ta, N, time_res)
    analog_bob = to_analog(tb, N, time_res)

    # cc = cross_corr_mine(analog_bob, analog_alice)
    cc = correlate(analog_bob, analog_alice, 'same')


    xt = np.arange(-total_time/2, total_time/2, time_res)
    offset_ind = np.argmax(np.real(cc)) 
    offset = xt[offset_ind]
    
    print('avg:::: ', np.average(np.real(cc)))
    print('std:::: ', np.std(np.real(cc)))
    
    
    # plt.figure()
    # plt.hist(np.real(cc), 500)
    
    if verbose:
        print('Offset:', offset)
       
    if plot:
        plt.figure()

        plt.plot(xt, np.real(cc), 'k-', alpha=0.5)
        plt.xlabel(r'$ \tau ~(\mu {s})$',  fontsize=14)
        plt.ylabel('$C$',  fontsize=14)
        plt.xlim([30, 150])
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        
    return offset, time_res


