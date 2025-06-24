# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 11:46:36 2023

@author: vil034
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time_tagging
import detector
import analog

import scipy.signal as sg
from scipy.io import savemat

import kanalysis as ka

####### Generate time tagging data
clock_shift = 46.123
total_protocol_time_sec = 1
alice_times, bob_times, Ns, Nb = time_tagging.generate_time_data(clock_shift, total_protocol_time_sec, loss=True, background_arg=1, jitter=True, verbose=True)

# np.save('alice_times', alice_times)
# np.save('bob_times', bob_times)

m_data = {'label':total_protocol_time_sec,
          'offset': clock_shift,
          'alice_times': alice_times,
          'bob_times': bob_times,
          'Ns': Ns,
          'Nb': Nb}


data_string = 'matlab/data/time_data.mat'
savemat(data_string, m_data)





# Partitioning time data
n = 9
dt = total_protocol_time_sec*1e6/n
alice_part = [[] for i in range(n)]
bob_part = x = [[] for i in range(n)]

for i in range(n):
    alice_part[i] = list(alice_times[np.logical_and(alice_times < dt*(i+1), alice_times >= dt*i)])
    bob_part[i] = list(bob_times[np.logical_and(bob_times < dt*(i+1), bob_times >= dt*i)])



# n = 6
# time_diffs = ka.find_peak_reverse(alice_times, bob_times, n)
# plt.hist(time_diffs, 500)

########################## iterative data generation
# for i in range(10): 
#     clock_shift = 46.123
#     total_protocol_time_sec = .615
#     alice_times, bob_times, Ns, Nb = time_tagging.generate_time_data(clock_shift, total_protocol_time_sec, loss=True, background_arg=0, jitter=True, verbose=True)
    
#     # np.save('alice_times', alice_times)
#     # np.save('bob_times', bob_times)
    
#     m_data = {'label':total_protocol_time_sec,
#               'offset': clock_shift,
#               'alice_times': alice_times,
#               'bob_times': bob_times,
#               'Ns': Ns,
#               'Nb': Nb}
    
    
#     data_string = 'matlab/data/time_data' + str(i+1) + '.mat'
#     savemat(data_string, m_data)






# alice_times = np.load('alice_times.npy')
# bob_times = np.load('bob_times.npy')
# 

####### Testing statistics
# et = list(entangled_pair_rate)
# aa = list(alice_times)


# plt.figure('Alice')
# plt.plot(entangled_pairs_time, 'r,')
# plt.plot(alice_times, 'o')

# plt.figure('Bob')
# plt.plot(entangled_pairs_time, 'r,')
# plt.plot(bob_times, 'o')




####################### Time data analysis
###### Time diferences parameters
# n0 = 0
# n = 2


# time_deltas0 = detector.time_diff_reverse(bob_times, alice_times, n, n0)
# ddt = 1.5
# idx = np.logical_and(time_deltas0 < ddt, time_deltas0 > -ddt)

# time_deltas2 = time_deltas0[idx]
# time_deltas = time_deltas0

# plt.close('all')

# plt.figure()
# bins = 100

# plt.hist(time_deltas0, bins)


#################### Analog analysis
# N = 2**22

# total_time = total_protocol_time_sec * 1e6
# time_res = total_time/N
# print('Time resolution: ', time_res)
# # alice_a = analog.to_analog(alice_times, time_res, start_time, total_protocol_time + start_time)
# # bob_a = analog.to_analog(bob_times, time_res, start_time, total_protocol_time + start_time)

# analog_alice = analog.to_analog2(alice_times, N, total_time)
# analog_bob = analog.to_analog2(bob_times, N, total_time)

# cc_dist = int(len(analog_alice)/1e3)

# cc = analog.cross_corr2(analog_alice, analog_bob)

# xt = np.arange(-total_time/2,total_time/2, time_res)
# # xt = np.arange(0,total_time, time_res)


# print('Offset:', xt[np.argmax(np.real(cc))])



# plt.figure()

# plt.plot(xt, np.real(cc))
# plt.xlabel('$c(t)~ \mu s$')


############## Analog analysis 2



# cc2 = sg.correlate(alice_a, bob_a, 'same')
# xt2 = np.arange(start_time - total_protocol_time,total_protocol_time + start_time, time_res)

# plt.figure()

# # plt.plot(xt[:-1], np.real(cc2))
# plt.plot(xt, np.real(cc2))

# plt.xlabel('$c(t)~ \mu s$')
# print('Offset:', xt[np.argmax(np.real(cc2))])







######### Fitting a gaussian curve to hist



# bins = 1000
# hist_data=np.histogram(time_deltas,bins)


# # Gaussian function for fitting 
# def Gaussian(t, amp, mean, std):
#     return amp*np.exp(-1/2*((t-mean)/std)**2)


# x_data = hist_data[1][:-1]
# y_data = hist_data[0]


# coarse_shift_guess = 1

# param_init = [np.max(y_data), coarse_shift_guess, (np.max(time_deltas) - np.min(time_deltas))/20]

# popt, pocv = curve_fit(Gaussian, x_data, y_data, p0=param_init)


# x=np.linspace(x_data[0], x_data[-1],1000)


# plt.figure()


# # plt.plot([popt[1],popt[1]],[0,popt[0]],color='red',linestyle='dotted')

# plt.legend()

# print('fitted_amplitude= {}counts'.format(popt[0]))

# print('time difference = {}s'.format(popt[1]))

# print('time resolution = {}s'.format(popt[2]))

# print('n(adjacent photons)={}'.format(n))


# plt.xlabel('Time Difference [mus]')

# plt.ylabel('Coincidences')



# plt.plot(x_data,y_data,color='black',label='Simulation',alpha=0.6)

# plt.plot(x, Gaussian(x ,popt[0],popt[1],popt[2]),label='Fitted Gaussian',color='red',linestyle='dotted')
# # plt.plot(x, Gaussian(x , param_init[0], param_init[1], param_init[2]),label='Fitted Gaussian',color='red',linestyle='dotted')

# # plt.xlim(-7.5e-9,7.5e-9)

# plt.legend()

# plt.savefig('histogram.jpg')

# plt.show()








