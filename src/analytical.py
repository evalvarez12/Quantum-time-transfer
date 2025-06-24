# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:18:21 2023

@author: vil034
"""

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

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
    
    Ns = eta*ra*T*time_filter_keep
    mub = ra*rb*dt*T
    
    
    
    P_failure = N0*.5*special.erfc((Ns)/np.sqrt(2*mub))
    P_failure[P_failure > 1] = 1 
    return  P_failure

#############################################################################
### Background counts plot

plt.close('all')

loss = 30
loss2 = 20
eta = 10**(-loss/10)
eta2 = 10**(-loss2/10)

qber_dev = 0.01
signal_rate = 3e6 # 3 MHz 
T = .1 # in seconds
sig = 500e-12
# sig = 2e-9

background_rate = np.logspace(3, 7)
ps = P_fail_4jitter(signal_rate, background_rate, eta, sig, T)
# ps = -10*np.log10(ps)
# ps2 = P_fail_2jitter(signal_rate, background_rate, eta, sig, T)
ps2 = P_fail_4jitter(signal_rate, background_rate, eta2, sig, T)

# qs = QBER_4sig(qber_dev, signal_rate, background_rate, eta, sig, T)
qs = QBER_2sig(qber_dev, signal_rate, background_rate, eta, sig, T)
qs2 = QBER_2sig(qber_dev, signal_rate, background_rate, eta2, sig, T)


fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))

ax2 = ax1.twinx()
ax2.plot(background_rate, ps, '--', color='tab:blue')
ax2.plot(background_rate, ps2, '-', color='tab:blue')

ax1.plot(background_rate, qs, '--', color='tab:orange', label= '30 dB', linewidth=2)
ax1.plot(background_rate, qs2, '-', color='tab:orange', label =  '20 dB', linewidth=2)

ax1.set_ylim(0.005, 0.26)
ax1.set_xscale('log')

bg_night = np.ones(15)*10000
bg_day = np.ones(15)*500000

ax1.plot(bg_night,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)
ax1.plot(bg_day,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)

ax1.text(11000, 0.06,'Night')
ax1.text(550000, 0.03,'Day')

ax1.set_xlabel('Background counts per second')
ax1.set_ylabel('QBER')
ax1.legend()
ax2.set_ylabel('Prob failure')
ax1.yaxis.label.set_color('tab:orange')
ax2.yaxis.label.set_color('tab:blue')

#############################################################################
### Loss plot

loss = np.arange(10,60)
eta = 10**(-loss/10)
qber_dev = 0.01
signal_rate = 3e6 # 3 MHz 
T = .1 # in seconds
sig = 500e-12
# sig = 2e-9

background_rate = 10000
background_rate2 = 500000

ps = P_fail_4jitter(signal_rate, background_rate, eta, sig, T)
# ps = -10*np.log10(ps)
# ps2 = P_fail_2jitter(signal_rate, background_rate, eta, sig, T)
ps2 = P_fail_4jitter(signal_rate, background_rate2, eta, sig, T)

qs = QBER_2sig(qber_dev, signal_rate, background_rate, eta, sig, T)
qs2 = QBER_2sig(qber_dev, signal_rate, background_rate2, eta, sig, T)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))

ax2 = ax1.twinx()
ax2.plot(loss, ps, '-', color='tab:blue', linewidth=2)
ax2.plot(loss, ps2, '--', color='tab:blue', linewidth=2)

ax1.plot(loss, qs, '-', color='tab:orange', label= '10 kHz bcps', linewidth=2)
ax1.plot(loss, qs2, '--', color='tab:orange', label='500 kHz bcps', linewidth=2)


loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)
ax1.plot(loss_downlink,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)
ax1.plot(loss_uplink,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)

ax1.text(26, 0.06,'10 km')
ax1.text(36,0.11,'Downlink')
ax1.text(56, 0.04,'Uplink')

ax1.set_ylim(0.005, 0.26)
# ax1.set_xscale('log')

ax1.set_xlabel('Loss (dB)')
ax1.set_ylabel('QBER')
ax1.legend()
ax2.set_ylabel('Prob failure')
ax1.yaxis.label.set_color('tab:orange')
ax2.yaxis.label.set_color('tab:blue')



#############################################################################
### Loss plot 2

loss = np.arange(10,60)
eta = 10**(-loss/10)
qber_dev = 0.01
signal_rate = 3e6 # 3 MHz 
T = 10 # in seconds
sig = 50e-12
# sig = 2e-9

background_rate = 10000
background_rate2 = 500000

ps = P_fail_4jitter(signal_rate, background_rate, eta, sig, T)
# ps = -10*np.log10(ps)
# ps2 = P_fail_2jitter(signal_rate, background_rate, eta, sig, T)
ps2 = P_fail_4jitter(signal_rate, background_rate2, eta, sig, T)

qs = QBER_2sig(qber_dev, signal_rate, background_rate, eta, sig, T)
qs2 = QBER_2sig(qber_dev, signal_rate, background_rate2, eta, sig, T)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))

ax2 = ax1.twinx()
ax2.plot(loss, ps, '-', color='tab:blue', linewidth=2)
ax2.plot(loss, ps2, '--', color='tab:blue', linewidth=2)

ax1.plot(loss, qs, '-', color='tab:orange', label= '10 kHz bcps', linewidth=2)
ax1.plot(loss, qs2, '--', color='tab:orange', label='500 kHz bcps', linewidth=2)


loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)
ax1.plot(loss_downlink,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)
ax1.plot(loss_uplink,np.linspace(0,max(qs2), 15), 'k--', linewidth=.7)

ax1.text(26, 0.06,'10 km')
ax1.text(36,0.11,'Downlink')
ax1.text(56, 0.04,'Uplink')

ax1.set_ylim(0.005, 0.26)
# ax1.set_xscale('log')

ax1.set_xlabel('Loss (dB)')
ax1.set_ylabel('QBER')
ax1.legend()
ax2.set_ylabel('Prob failure')
ax1.yaxis.label.set_color('tab:orange')
ax2.yaxis.label.set_color('tab:blue')

#############################################################
## Time jitter plot
# loss = 30
# eta = 10**(-loss/10)
# qber_dev = 0.01
# signal_rate = 3e6 # 3 MHz 
# T = .1 # in seconds
# sig = np.linspace(200e-12, 3e-9)
# background_rate = 10000

# ps = P_fail_4jitter(signal_rate, background_rate, eta, sig, T)
# # ps = -10*np.log10(ps)
# ps2 = P_fail_2jitter(signal_rate, background_rate, eta, sig, T)

# qs = QBER_4sig(qber_dev, signal_rate, background_rate, eta, sig, T)
# qs2 = QBER_2sig(qber_dev, signal_rate, background_rate, eta, sig, T)

# fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))

# ax2 = ax1.twinx()
# ax2.plot(sig, ps, '-', color='tab:blue')
# ax2.plot(sig, ps2, '--', color='tab:blue')

# ax1.plot(sig, qs, '-', color='tab:orange', label= '10 kHz cps')
# ax1.plot(sig, qs2, '--', color='tab:orange', label='500 kHz cps')

# # ax1.set_ylim(0, 0.26)
# # ax1.set_xscale('log')

# ax1.set_xlabel('Background counts per second')
# ax1.set_ylabel('QBER')
# ax1.legend()
# ax2.set_ylabel('Prob failure')
# ax1.yaxis.label.set_color('red')
# ax2.yaxis.label.set_color('b')

#################################################################


