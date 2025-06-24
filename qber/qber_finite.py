# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 10:25:59 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def avg_std(arr):
    return np.average(arr), np.std(arr)

qber_real = 0.22

qber_limit = 0.25


time = 10

rate = 1000 # coincidence rate per second

updates_every_time = 0.1
cumulative_time = 1
N_update = int(rate*updates_every_time)

N_update_block = int(cumulative_time/updates_every_time)

time_x = np.arange(0, time, updates_every_time)
N = rate * time

roll = np.random.rand(N)

qbits_err = roll < qber_real
qbits_err = qbits_err.astype(int)

cumulative = True

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
ax1.text(6.1, .251, 'Security limit', color='r', fontsize=18)

ax1.set_xlabel('Time (seg)', fontsize=18)
ax1.set_ylabel('QBER', fontsize=18)
ax1.set_ylim([0.18,0.26])

# ax2.set_ylabel('Confidence')
# ax2.set_ylim([0.9, 1.005])

# ax2.plot(time_x, confidence, 'r--', alpha = 0.3)
# ax2.plot(time_x[9::10], confidence[9::10], 'ro', alpha = 0.8)

# ax2.set_yscale("log")
# plt.ylim([0.2, 0.26])

plt.show()    
    