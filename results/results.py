# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:07:05 2023

@author: vil034
"""

###################################################
###################################################
###################################################
## PLOTS FROM THIS SCRIPT ARE NOT USED IN THE PAPER

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

def QBER_2sig(qber_dev, sig_rate, back_rate, eta, jitter, time):
    time_filter = 2*jitter
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*sig_rate*time
    Nb = sig_rate*back_rate*time*time_filter
    
    return (Ns*qber_dev + Nb/2)/(Ns + Nb)


background_day = 500000
background_night = 10000

loss_det = 10**(-16/10)

jitter = 500e-12

time = 0.1

sig_rate = 20e6

qber_dev = 0.05

file_day = 'data/HORIZONTAL_DISTANCE_trx=15_L=%s_wl=0.81'
file_night = 'data/HORIZONTAL_DISTANCE_NIGHT_trx=15_L=%s_wl=0.81'

dists = ['7000', '8000', '9000', '10000', '11000', '12000', '13000', '14000', '15000', '16000']
distsi = [7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000]
m_d = []
std_d = []
m_n = []
std_n = []

for d in dists:
    data_day = sp.io.loadmat(file_day %d)
    data_day = data_day['res']
    # data_day = -10*np.log10(data_day)
    
    m_d += [np.average(data_day)]
    std_d += [np.std(data_day)]
    
    data_n = sp.io.loadmat(file_night %d)
    data_n = data_n['res']
    # data_n = -10*np.log10(data_n)
    
    m_n += [np.average(data_n)]
    std_n += [np.std(data_n)]
    
    
loss_day = np.array(m_d)*loss_det
loss_night = np.array(m_n)*loss_det


qber_day = QBER_2sig(qber_dev, sig_rate, background_day, loss_day, jitter, time)
qber_night = QBER_2sig(qber_dev, sig_rate, background_night, loss_night, jitter, time)

fig = plt.figure()
plt.plot(distsi, qber_day, label='Day', linewidth=3)
plt.plot(distsi, qber_night, label='Night', linewidth=3)
# plt.xlabel('Distance (seg)', fontsize=15)
# plt.ylabel('QBER', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

fig.set_size_inches(11.5*.6, 10.5*.6)
fig.savefig('qber_main.png', dpi=200)

####################################################

def avg_std(arr):
    return np.average(arr), np.std(arr)

qber_real = qber_day[3]
# qber_real = qber_day[5]

qber_limit = 0.25


time = 10

rate = int(sig_rate*loss_day[3]) 
# rate  = int(sig_rate*loss_day[5])

updates_every_time = 0.1
cumulative_time = 1
N_update = int(rate*updates_every_time)

N_update_block = int(cumulative_time/updates_every_time)

time_x = np.arange(0, time, updates_every_time)
N = rate * time

roll = np.random.rand(N)

qbits_err = roll < qber_real
qbits_err = qbits_err.astype(int)

cumulative = False

avgs = []
stds = []
for i in range(int(time/updates_every_time)):
    Ni = N_update*(i+1)
    Ni0 = int(i/N_update_block)*N_update*N_update_block
    
    if cumulative:
        a, s = avg_std(qbits_err[:Ni])
        Ni0 = 0
    else:
        a, s = avg_std(qbits_err[Ni0:Ni])
        # print(Ni0,Ni, Ni - Ni0)
        
    avgs += [a]
    stds += [s/np.sqrt(Ni - Ni0)]
    
    
avgs = np.array(avgs)
stds = np.array(stds)    
confidence = 1-sp.special.erfc((qber_limit-avgs)/(np.sqrt(2)*stds))/2
    
  

fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

ax1.errorbar(time_x,avgs, stds)
ax1.plot(time_x, np.ones_like(time_x)*qber_limit, 'r--', alpha = 0.2)    
# ax1.text(6.1, .251, 'Security limit', color='r', fontsize=15)

# ax1.set_xlabel('Time (seg)', fontsize=15)
# ax1.set_ylabel('QBER', fontsize=15)
ax1.set_ylim([0.18,0.26])

ax1.tick_params(axis='both', which='major', labelsize=13, labelrotation=0)
# ax2.set_ylabel('Confidence')
# ax2.set_ylim([0.9, 1.005])

# ax2.plot(time_x, confidence, 'r--', alpha = 0.3)
# ax2.plot(time_x[9::10], confidence[9::10], 'ro', alpha = 0.8)

# ax2.set_yscale("log")
# plt.ylim([0.2, 0.26])

fig.set_size_inches(18.5*.3, 10.5*.3)
fig.savefig('qber2.png', dpi=200)

plt.show()    
    