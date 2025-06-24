# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:12:14 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import analytical_decoy as decoy
from scipy import special



def P_fail_ES_4jitter(ra, rb, channel_loss, detector_loss, sig, T):
    eta = 10**(-(channel_loss + detector_loss)/10)
    ra = ra*10**(-detector_loss/10)
    return P_fail_4jitter(ra, rb, eta, sig, T)
    

def P_fail_ES_2jitter(ra, rb, channel_loss, detector_loss, sig, T):
    eta = 10**(-(channel_loss + detector_loss)/10)
    ra = ra*10**(-detector_loss/10)
    return P_fail_2jitter(ra, rb, eta, sig, T)

    
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


def QBER_ES(qber_dev, loss, time_jitter, background_rate):
    eta = 10**(-loss/10)
    return QBER_2sig(qber_dev, background_rate, eta, time_jitter)


def QBER_2sig(qd, rb, eta, sig):
    time_filter = 2*sig
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta
    Nb = rb*time_filter
    
    return (Ns*qd + Nb/2)/(Ns + Nb)








def QBER_BB84(qber_dev, mu, loss, time_jitter, background_rate):  
    eta = 10**(-loss/10)
    eta = (1 - np.exp(-eta*mu))
    
    return QBER_2sig(qber_dev, background_rate, eta, time_jitter)
    


def QBER_BB84_decoy(qber_dev, mu, nu, loss, time_jitter, background_rate):
    eta = 10**(-loss/10)
    eta = (1 - np.exp(-eta*mu))
    
    Y0 = background_rate*2*time_jitter
    # Ev = QBER_BB84(qber_dev, mu, loss, time_jitter, background_rate)
    Ev = (Y0/2 + qber_dev*(1-np.exp(-eta*nu)))/(Y0 + 1-np.exp(-eta*nu))
    # Time filter
    
    return np.abs(decoy.e1(nu, mu, eta, Y0, Ev))
    

plt.close('all')

#################################################
# background_rate = 10000
# protocol_time = np.linspace(0.01, 1)
# qs, ns = QBER(qber_dev, signal_rate, loss_channel, loss_detector, time_jitter, background_rate, protocol_time*1e6)

# fig, ax1 = plt.subplots()

# ax2 = ax1.twinx()
# ax1.plot(protocol_time, qs, 'r')
# ax2.plot(protocol_time, ns, 'b')

#################################################
background_rate = np.logspace(3, 6)
protocol_time = 0.1
mu = 0.5
nu = 0.18
mu2 = 0.1
qber_dev =  0.01
qber_dev_es = 2/3 * (0.05)
qber_dev_es2 = 2/3 * (0.01)

loss = 25
time_jitter = 500e-12

qs = QBER_ES(qber_dev_es, loss, time_jitter, background_rate)
qs2 = QBER_ES(qber_dev_es2, loss, time_jitter, background_rate)

qs_bb84 = QBER_BB84_decoy(qber_dev, mu, nu, loss, time_jitter, background_rate)
qs_bb842 = QBER_BB84(qber_dev, mu2, loss, time_jitter, background_rate)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))





# ax2 = ax1.twinx()
ax1.plot(background_rate, qs, '-', label='entangled source - $F=95\%$', linewidth=2)
ax1.plot(background_rate, qs2, '-', label='entangled source - $F=99\%$', linewidth=2)

ax1.plot(background_rate, qs_bb84, '--', label='BB84 $\mu=0.5 - decoy ~state$', linewidth=2)
ax1.plot(background_rate, qs_bb842, '--', label='BB84 $\mu=0.1$', linewidth=2)

bg_night = np.ones(15)*10000
bg_day = np.ones(15)*500000

ax1.plot(bg_night,np.linspace(0,max(qs_bb842), 15), 'k--', linewidth=.7)
ax1.plot(bg_day,np.linspace(0,max(qs_bb842), 15), 'k--', linewidth=.7)

ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(background_rate), max(background_rate))

ax1.text(11000, 0.06,'Night')
ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Background counts per second')
ax1.set_ylabel('QBER')
ax1.legend()
# ax2.plot(background_rate, ns, '-')
# ax2.plot(background_rate, ns_bb84, '-')

###############################################
background_rate = 1e4
protocol_time = 0.1
mu = 0.5
nu = 0.18
mu2 = 0.1
qber_dev =  0.01
qber_dev_es = 2/3 * (0.05)
qber_dev_es2 = 2/3 * (0.01)
loss = np.arange(10,60)
time_jitter = 500e-12

qs = QBER_ES(qber_dev_es, loss, time_jitter, background_rate)
qs2 = QBER_ES(qber_dev_es2, loss, time_jitter, background_rate)

qs_bb84 = QBER_BB84_decoy(qber_dev, mu, nu, loss, time_jitter, background_rate)
qs_bb842 = QBER_BB84(qber_dev, mu2, loss, time_jitter, background_rate)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))





# ax2 = ax1.twinx()
ax1.plot(loss, qs, '-', label='entangled source - $F=95\%$', linewidth=2)
ax1.plot(loss, qs2, '-', label='entangled source - $F=99\%$', linewidth=2)

ax1.plot(loss, qs_bb84, '--', label='BB84 $\mu=0.5 - decoy ~state$', linewidth=2)
ax1.plot(loss, qs_bb842, '--', label='BB84 $\mu=0.1$', linewidth=2)

loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
ax1.plot(loss_downlink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
ax1.plot(loss_uplink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)

ax1.text(26, 0.06,'10 km')
ax1.text(36,0.11,'Downlink')
ax1.text(56, 0.04,'Uplink')

# ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(loss), max(loss))

ax1.text(11000, 0.06,'Night')
ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Loss (dB)')
ax1.set_ylabel('QBER')
ax1.legend()


###############################
background_rate1 = 1e4
background_rate2 = 5e5

protocol_time = 0.1
qber_dev_es = 2/3 * (0.05)
qber_dev_es2 = 2/3 * (0.01)
loss = np.arange(10,60)
time_jitter = 500e-12

ra = 3e6
ra2 = 25e6

qs = QBER_ES(qber_dev_es, loss+3, time_jitter, background_rate2)
qs2 = QBER_ES(qber_dev_es2, loss+3, time_jitter, background_rate2)

pfail1 = P_fail_ES_2jitter(ra, background_rate2, loss, 3, time_jitter, protocol_time)
pfail2 = P_fail_ES_2jitter(ra2, background_rate2, loss, 3, time_jitter, protocol_time)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))





ax2 = ax1.twinx()
ax1.plot(loss, qs, '-', color='tab:orange', label='$F=95\%$', linewidth=2)
ax1.plot(loss, qs2, '--', color='tab:orange', label='$F=99\%$', linewidth=2)

ax2.plot(loss, pfail1, '-', color='tab:blue', label='$ r_a= 3 MHz $', linewidth=2)
ax2.plot(loss, pfail2, '--', color='tab:blue', label='$r_a = 25 MHz $', linewidth=2)

loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
ax1.plot(loss_downlink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
ax1.plot(loss_uplink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)

ax1.text(26, 0.06,'10 km')
ax1.text(36,0.11,'Downlink')
ax1.text(56, 0.04,'Uplink')

# ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(loss), max(loss))

ax1.text(11000, 0.06,'Night')
ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Channel loss')
ax1.set_ylabel('QBER')
ax1.legend()
ax2.legend(loc='lower right')

ax2.set_ylabel('Prob failure')
ax1.yaxis.label.set_color('tab:orange')
ax2.yaxis.label.set_color('tab:blue')



###############################################
background_rate = 203000*2.4
F1 = 0.9
F2 = 0.95
F3 = 0.99

qber1 = 2/3 * (1-F1)
qber2 = 2/3 * (1-F2)
qber3 = 2/3 * (1-F3)

protocol_time = 0.1
loss = np.arange(10,40)
time_jitter = 500e-12

qs1 = QBER_ES(qber1, loss, time_jitter, background_rate)
qs2 = QBER_ES(qber2, loss, time_jitter, background_rate)
qs3 = QBER_ES(qber3, loss, time_jitter, background_rate)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))


# ax2 = ax1.twinx()
ax1.plot(loss, qs1, '-', label='entangled source - $F=90\%$', linewidth=2)
ax1.plot(loss, qs2, '-', label='entangled source - $F=95\%$', linewidth=2)
ax1.plot(loss, qs3, '-', label='entangled source - $F=99\%$', linewidth=2)

loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
# ax1.plot(loss_downlink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
# ax1.plot(loss_uplink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)

ax1.text(26, 0.06,'10 km')
# ax1.text(33,0.11,'Downlink')
# ax1.text(56, 0.04,'Uplink')

# ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(loss), max(loss))

# ax1.text(11000, 0.06,'Night')
# ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Loss (dB)')
ax1.set_ylabel('QBER')
ax1.set_title('Entangled pairs')

ax1.legend()


###############################################
background_rate = 203000
protocol_time = 0.1
mu = 0.5
nu = 0.18
mu2 = 0.1
mu3 = 1
qber_dev =  0.01
loss = np.arange(10,40)
time_jitter = 500e-12


qs_bb84 = QBER_BB84_decoy(qber_dev, mu, nu, loss, time_jitter, background_rate)
qs_bb843 = QBER_BB84_decoy(qber_dev, mu3, nu, loss, time_jitter, background_rate)

qs_bb842 = QBER_BB84(qber_dev, mu2, loss, time_jitter, background_rate)

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))





# ax2 = ax1.twinx()
ax1.plot(loss, qs_bb842, '--', label='BB84 $\mu=0.1$', linewidth=2)

ax1.plot(loss, qs_bb84, '--', label='BB84 $\mu=0.5 - decoy ~state$', linewidth=2)
ax1.plot(loss, qs_bb843, '--', label='BB84 $\mu=1 - decoy ~state$', linewidth=2)

loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
# ax1.plot(loss_downlink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
# ax1.plot(loss_uplink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)

ax1.text(26, 0.06,'10 km')
# ax1.text(36,0.11,'Downlink')
# ax1.text(56, 0.04,'Uplink')

# ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(loss), max(loss))

# ax1.text(11000, 0.06,'Night')
# ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Loss (dB)')
ax1.set_ylabel('QBER')
ax1.set_title('BB84 source')
ax1.legend()

###############################################
background_rate1 = 5100*2.4
background_rate2 = 203000*2.4
F1 = 0.97


qber1 = 2/3 * (1-F1)
qber2 = 2/3 * (1-F2)
qber3 = 2/3 * (1-F3)

protocol_time = 0.1
loss = np.arange(10,40)
time_jitter = 500e-12

qs1 = QBER_ES(qber1, loss, time_jitter, background_rate1)
qs2 = QBER_ES(qber1, loss, time_jitter, background_rate2)


plt.rc('font', size=18) 

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))


# ax2 = ax1.twinx()
ax1.plot(loss, qs1, '-', label='Nighttime', linewidth=3)
ax1.plot(loss, qs2, '-', label='Daytime', linewidth=3)

loss_10k = np.ones(15)*25
loss_downlink = np.ones(15)*35
loss_uplink = np.ones(15)*55


ax1.plot(loss_10k,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
ax1.plot(loss_downlink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
ax1.plot(loss_uplink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)

# ax1.text(25.2, 0.06,'10 km')
# ax1.text(35.2,0.11,'Downlink')
# ax1.text(56, 0.04,'Uplink')

# ax1.set_xscale('log')

ax1.set_ylim(0, 0.26)
ax1.set_xlim(min(loss), max(loss))

# ax1.text(11000, 0.06,'Night')
# ax1.text(300000, 0.03,'Day')

# ax1.set_xlabel('Loss (dB)')
# ax1.set_ylabel('QBER')

ax1.legend()


###############################################
loss810 = np.array([10.7, 11.8, 12.8, 13, 13.8, 14.5]) + 7.7
loss1550 = np.array([19, 20, 20.9, 21.8, 22.6, 23.2]) + 9.5
ds = [0.7, 0.8, 0.9, 1, 1.1, 1.2]

background_rate810 = 187000 
background_rate1550 = 89000 
F = 0.95


qber = 2/3 * (1-F)
protocol_time = 0.1
time_jitter = 500e-12

qs810 = QBER_ES(qber, loss810, time_jitter, background_rate810)
qs1550 = QBER_ES(qber, loss1550, time_jitter, background_rate1550)


plt.rc('font', size=18) 

fig, ax1 = plt.subplots(figsize=(1.6*6, 1*6))


# ax2 = ax1.twinx()
ax1.plot(ds, qs810, '-', label='810', linewidth=3)
ax1.plot(ds, qs1550, '-', label='1550', linewidth=3)


# ax1.plot(loss_10k,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
# ax1.plot(loss_downlink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)
# ax1.plot(loss_uplink,np.linspace(0,max(qs), 15), 'k--', linewidth=.7)

# ax1.text(25.2, 0.06,'10 km')
# ax1.text(35.2,0.11,'Downlink')
# ax1.text(56, 0.04,'Uplink')

# ax1.set_xscale('log')

# ax1.set_ylim(0, 0.26)
# ax1.set_xlim(min(loss), max(loss))

# ax1.text(11000, 0.06,'Night')
# ax1.text(300000, 0.03,'Day')

ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('QBER')

ax1.legend()
