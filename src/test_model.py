# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:10:31 2023

@author: vil034
"""

import model as mod
import numpy as np
import cross_corr as cc
import matplotlib.pyplot as plt
import fine_analysis as fa
from scipy import special

plt.close('all')

new_data = True

if new_data:

    setup = mod.Model()
    
    background_rate = 20000
    qber = 0.01
    loss_channel = 25
    loss_detector = 5
    signal_rate = 3
    
    setup.set_background_rate(background_rate)
    setup.set_qber(qber)
    setup.set_loss(loss_channel, loss_detector)
    # setup.set_jitter(0)
    setup.set_signal_rate(signal_rate)
    
    
    time_shift = 46.123
    start_time = 0
    total_protocol_time_sec = .1
    
    
    
    total_protocol_time = total_protocol_time_sec * 1e6
    ta, tb, qa, qb = setup.generate_time_data(time_shift, start_time, total_protocol_time_sec, verbose=True)
    
    
    # Number of bins 2**n
    n = 25
    offset, time_res = cc.analysis(ta, tb, n, total_protocol_time, True, True)
    
    
    
    ####### CC analysis again
    # time2 = max(tb2 - ta2) - min(tb2 - ta2)
    
    # offset2, time_res2 = cc.analysis(ta2, tb2, n, time2, True, True)
    
    
    print('----------- Fine analysis')
    
    ta2, tb2, qa2, qb2 = fa.get_overlap(ta, tb, qa, qb, offset, time_res)
    print('Data points: ', len(ta2))
    
    np.save('alice_time_points', ta2)
    np.save('bob_time_points', tb2)
    np.save('alice_q_points', qa2)
    np.save('bob_q_points', qb2)




alice_times = np.load('alice_time_points.npy')
bob_times = np.load('bob_time_points.npy')
alice_q = np.load('alice_q_points.npy')
bob_q = np.load('bob_q_points.npy')

diff = bob_times - alice_times

pfit = fa.fit_w_hist(diff, 200)


idx = np.where(np.abs(diff - pfit[0]) < 2*np.abs(pfit[1]) )

qber = sum(alice_q[idx] != bob_q[idx])/len(alice_q[idx])


# plt.hist(diff[idx])

print('Offset: ', pfit[0])
print('Precision: ', pfit[1])
print('Measured QBER: ', qber)

print('-------------------- Theory')

time_jitter = 500e-6
N = total_protocol_time/(time_jitter)
dt = total_protocol_time/N
# each bin has the width of time jitter - time_jitter is the time unit
Ns = signal_rate*total_protocol_time*(10**(-(loss_channel + 2*loss_detector)/10))
ra = signal_rate*total_protocol_time*(10**(-loss_detector/10))/N
rb = background_rate*total_protocol_time_sec/N
p = ra*rb
mu_b = N*p
sig_b2 = N*p*(1-p)
print('Avg noise: ', mu_b)
print('Std noise: ', sig_b2)

P_failure = .5*special.erfc((Ns - mu_b)/np.sqrt(2*sig_b2))
print('Probability of large bin: ', P_failure)
print('Probability of failure: ', P_failure*N)    


n_background_in_peak = ra*N*background_rate*1e-6*2*time_jitter 
qber_est = (Ns*qber + n_background_in_peak/2)/(Ns + n_background_in_peak)
print('QBER: ', qber_est)

# ###############################################################
# print('----------- Fine analysis 2')


# # m1 = (diff - pfit[1]) >  - 2*pfit[2]
# # m2  = (diff - pfit[1]) < 2*pfit[2]

# # m = np.logical_and(m1, m2)

# m = np.where(np.abs(diff - pfit[1]) < pfit[2])

# diff2 = diff[m]
# qa3 = alice_q[m]
# qb3 = bob_q[m]

# # if len()

# print('Data points: ', len(diff2))



# pfit2 = fa.fit_w_hist(diff2, 200)
# qber = sum(qa3 != qb3)/len(qa3)

# print('Offset: ', pfit2[1])
# print('Precision: ', pfit2[2])
# print('Floor: ', pfit2[3])
# print('Measured QBER: ', qber)

# print('Mean: ', np.average(diff2))

# print('Std: ', np.std(diff2))

