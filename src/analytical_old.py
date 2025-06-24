# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:12:14 2023

@author: vil034
"""

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

def P_fail(signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time):
    N = protocol_time/(time_jitter)
    # each bin has the width of time jitter - time_jitter is the time unit
    Ns = signal_rate*protocol_time*(10**(-(loss_channel + 2*loss_detector)/10))
    ra = signal_rate*protocol_time*(10**(-loss_detector/10))/N
    rb = background_rate*total_protocol_time_sec/N
    p = ra*rb
    mu_b = N*p
    sig_b2 = N*p*(1-p)
    # print('Avg noise: ', mu_b)
    # print('Std noise: ', sig_b2)

    P_failure = .5*special.erfc((Ns - mu_b)/np.sqrt(2*sig_b2))
    # print('Probability of large bin: ', P_failure)
    # print('Probability of failure: ', P_failure*N)    
    
    P_failure = P_failure*500*1e3
    P_failure[P_failure > 1] = 1  
    # P_failure[P_failure == 0] = 1e-50    

    return P_failure

def QBER(qber_dev, signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time):  
    Ns = signal_rate*protocol_time*(10**(-(loss_channel + 2*loss_detector)/10))
    Na = signal_rate*protocol_time*(10**(-loss_detector/10))

    time_filter = 2*time_jitter
    n_background_in_peak = Na*background_rate*1e-6*time_filter
    qber_est = (Ns*qber_dev + n_background_in_peak/2)/(Ns + n_background_in_peak)
    return qber_est, Ns

def QBER_BB84(qber_dev, signal_rate, mu, loss_channel, loss_detector, time_jitter, background_rate, protocol_time):  
    transmissivity = 10**(-(loss_channel + loss_detector)/10)
    Ns = signal_rate*protocol_time*(1 - np.exp(-transmissivity*mu))
    Na = signal_rate*protocol_time

    time_filter = 2*time_jitter
    n_background_in_peak = Na*background_rate*1e-6*time_filter
    qber_est = (Ns*qber_dev + n_background_in_peak/2)/(Ns + n_background_in_peak)
    return qber_est, Ns


plt.close('all')

background_rate = 100000
qber_dev = 0.01
loss_channel = 25
loss_detector = 5
signal_rate = 3


total_protocol_time_sec = .1
protocol_time = total_protocol_time_sec * 1e6

time_jitter = 500e-6


background_rate = np.linspace(1000, 1.4e6)
ps = P_fail(signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time)
# ps = -10*np.log10(ps)

qs, ns = QBER(qber_dev, signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time)
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax2.plot(background_rate, ps)
ax1.plot(background_rate, qs, 'r')

ax1.set_ylim(0, 0.26)
# ax1.set_xscale('log')

ax1.set_xlabel('Background counts per second')
ax1.set_ylabel('QBER')
ax2.set_ylabel('Prob failure')
ax1.yaxis.label.set_color('red')
ax2.yaxis.label.set_color('b')
#################################################
# background_rate = 10000
# protocol_time = np.linspace(0.01, 1)
# qs, ns = QBER(qber_dev, signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time*1e6)

# fig, ax1 = plt.subplots()

# ax2 = ax1.twinx()
# ax1.plot(protocol_time, qs, 'r')
# ax2.plot(protocol_time, ns, 'b')

#################################################
background_rate = np.linspace(1000, 600000)
protocol_time = 0.1
mu = 0.5
mu2 = 0.1

loss_detector = 1

qs, ns = QBER(qber_dev, signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time*1e6)
qs_bb84, ns_bb84 = QBER_BB84(qber_dev, 3, mu, loss_channel, loss_detector, time_jitter, background_rate, protocol_time*1e6)
qs_bb842, ns_bb842 = QBER_BB84(qber_dev, 3, mu2, loss_channel, loss_detector, time_jitter, background_rate, protocol_time*1e6)

fig, ax1 = plt.subplots()

bg_night = np.ones(15)*10000
bg_day = np.ones(15)*500000



# ax2 = ax1.twinx()
ax1.plot(background_rate, qs, '-', label='entangled source')

ax1.plot(background_rate, qs_bb84, '--', label='BB84 $\mu=0.5$')
ax1.plot(background_rate, qs_bb842, '--', label='BB84 $\mu=0.1$')

ax1.plot(bg_night,np.linspace(0,max(qs_bb842), 15), 'k--')
ax1.plot(bg_day,np.linspace(0,max(qs_bb842), 15), 'k--')

ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(background_rate), max(background_rate))

ax1.text(11000, 0.1,'Night')
ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Background counts per second')
ax1.set_ylabel('QBER')
ax1.legend()
# ax2.plot(background_rate, ns, '-')
# ax2.plot(background_rate, ns_bb84, '-')