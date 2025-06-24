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
c = 299792458 # m/s
clock_shift = 46.123
total_protocol_time_sec = .2
alice_times, bob_times, Ns, Nb = time_tagging.generate_time_data(clock_shift, total_protocol_time_sec, loss=True, background_arg=0, jitter=True, verbose=True)

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
n = 1
dt = total_protocol_time_sec*1e6/n
alice_part = [[] for i in range(n)]
bob_part = x = [[] for i in range(n)]

for i in range(n):
    alice_part[i] = list(alice_times[np.logical_and(alice_times < dt*(i+1), alice_times >= dt*i)])
    bob_part[i] = list(bob_times[np.logical_and(bob_times < dt*(i+1), bob_times >= dt*i)])


# m_data_part = {'label':total_protocol_time_sec,
#                'offset': clock_shift,
#                'alice_times': alice_part,
#                'bob_times': bob_part,
#                'Ns': Ns,
#                'Nb': Nb}


# data_string = 'matlab/data/time_data_part.mat'
# savemat(data_string, m_data_part)



#################### Analog analysis
N = 2**24

time_res = dt/N
print('Time partition: ', dt)
print('Time resolution: ', time_res)
print('Desired resolution: ', 10e-6)
print('Probability pair lost: ', clock_shift/dt)
print('Pairs lost: ', clock_shift/dt *Ns)
cc = np.zeros(N)

for i in range(n):
    
    analog_alice = analog.to_analog(alice_part[i], N, dt)
    analog_bob = analog.to_analog(bob_part[i], N, dt)

    cc += sg.correlate(analog_bob, analog_alice, 'same')

    # xt = np.arange(0,total_time, time_res)

xt = np.arange(-dt/2,dt/2, time_res)
offset_ind = np.argmax(np.real(cc)) 
offset = xt[offset_ind]
print('Offset:', offset)

 # bob_times = bob_times offset + time_res


# plt.figure()

# plt.plot(xt, np.real(cc))
# plt.xlabel('$c(t)~ \mu s$')




def get_overlap(va, vb, shift, margin):
    ao = []
    bo = []
    for a in va:
        m1 = a > vb + shift - margin
        m2 = a < vb + shift + margin
        
        m = np.logical_and(m1, m2)
        if any(m):
            ao += [a]
            bo += list(vb[m])
            
    return np.array(ao), np.array(bo)

# def get_overlap(va, vb, margin=0):
#     return [a for a in set(va) if 
#         any(a > (b - margin) and a < (b + margin)
#         for b in set(vb))]



aa, bb = get_overlap(bob_times, alice_times, offset, time_res)

np.save('alice_points', aa)
np.save('bb_points', bb)
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








