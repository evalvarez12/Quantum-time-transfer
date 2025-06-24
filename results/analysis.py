# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 12:34:30 2023

@author: vil034
"""


import numpy as np
from scipy import special
import matplotlib.pyplot as plt

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

def QBER_2sig(qd, ra, rb, eta, sig, T):
    time_filter = 2*sig
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*ra*T
    Nb = rb*ra*T*time_filter
    
    return (Ns*qd + Nb/2)/(Ns + Nb)

def QBER_4sig(qd, ra, rb, eta, sig, T):
    time_filter = 4*sig
    time_filter_keep = 0.95
    
    Ns = time_filter_keep*eta*ra*T
    Nb = rb*ra*T*time_filter
    
    return (Ns*qd + Nb/2)/(Ns + Nb)

    

def P_fail_4jitter(ra, rb, eta, sig, T):
    dt = 4*sig
    time_filter_keep = 0.95
    # dt = 2*sig
    # time_filter_keep = 0.68
    # number of bins in the coarse guess
    # N0 = 500
    N0 = T/dt
    
    Ns = eta*ra*T*time_filter_keep
    mub = ra*rb*dt*T
    
    
    
    P_failure = N0*.5*special.erfc((Ns)/np.sqrt(2*mub))
    P_failure[P_failure > 1] = 1 
    return  P_failure

def P_fail_2jitter(ra, rb, eta, sig, T):
    dt = 2*sig
    time_filter_keep = 0.68
    # dt = 2*sig
    # time_filter_keep = 0.68
    # number of bins in the coarse guess
    # N0 = 500
    N0 = T/dt
    # N0 = 1e2
    
    Ns = eta*ra*T*time_filter_keep
    mub = ra*rb*dt*T
    
    
    
    P_failure = N0*.5*special.erfc((Ns)/np.sqrt(2*mub))
    P_failure[P_failure > 1] = 1 
    return  P_failure


#############################################################################
### Loss plot

loss = np.arange(10,50)
eta = 10**(-loss/10)
qber_dev = 0.05
signal_rate = 3e6 # 3 MHz 
T = .1 # in seconds
sig = 500e-12
# sig = 2e-9

background_rate = np.array([10000, 80000, 20000, 350000, 500000])

background_rate_labs = ['10k', '80k', '200k', '350k', '500k']


ps = np.zeros((len(background_rate), len(loss)))
qs = np.zeros((len(background_rate), len(loss)))

for i in range(len(background_rate)):
    ps[i, :] = P_fail_2jitter(signal_rate, background_rate[i], eta, sig, T)

    qs[i, :] = QBER_2sig(qber_dev, signal_rate, background_rate[i], eta, sig, T)
    
    
fig = plt.figure(figsize=(1.6*4, 1*4))

for i in range(len(background_rate)):
    plt.plot(loss, ps[i], linewidth=2, label=background_rate_labs[i])
    
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.legend(fontsize=20)
# plt.yscale('log')
# plt.ylim([0, 1])

fig = plt.figure(figsize=(1.6*4, 1*4))

for i in range(len(background_rate)):
    plt.plot(loss, qs[i], linewidth=2, label=str(int(background_rate[i])))
    
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
# plt.legend(fontsize=15)
plt.ylim([0.045, 0.255])
