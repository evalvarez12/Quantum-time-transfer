# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 11:18:58 2023

@author: vil034
"""

import numpy as np
from scipy.fft import fft, ifft
from scipy.signal import correlate

def to_analog2(time_data, time_res, time_start, time_stop):
    total_time = time_stop - time_start
    
    # Number of bins in analog data
    N = np.ceil(total_time/time_res)
    
    analog_data = np.zeros(int(N))
    
    for ti in time_data:
         idx = np.floor(np.mod(ti/time_res, N))
         analog_data[int(idx)] = 1
    
    return analog_data


def to_analog(time_data, N, total_time):

    time_res = total_time/N
    
    analog_data = np.zeros(int(N))
        
    for ti in time_data:
         idx = np.floor(np.mod(ti/time_res, N))
         analog_data[int(idx)] = 1
    
    return analog_data



def cross_corr_mine(a, b):
    return ifft(np.conj(fft(a))*fft(b))



# t1 = np.array([0,1, 4, 5, 6, 7, 9, 10, 12])
# analog_t1 = to_analog(t1, 1, 0, 20)
# print(analog_t1)


# t2 = np.array([0, 1, 4, 5, 6, 7, 9, 10, 12,])
# analog_t2 = to_analog(t2, 1, 0, 20)
# print(analog_t2)


# cc = cross_corr(analog_t1, analog_t2)
# print(cc)
# # 
# plt.hist(np.real(cc))
