# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 09:00:17 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as spc

def QBER_2sig(qber_dev, sig_rate, back_rate, eta, jitter, time):
    time_filter = 2*jitter
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*sig_rate*time
    Nb = sig_rate*back_rate*time*time_filter
    
    return (Ns*qber_dev + Nb/2)/(Ns + Nb)


def QBER_2sig(qber_dev, sig_rate, back_rate, eta, jitter, time):
    time_filter = 2*jitter
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*sig_rate*time
    Nb = sig_rate*back_rate*time*time_filter
    
    return (Ns*qber_dev + Nb/2)/(Ns + Nb)

def avg_std(arr):
    return np.average(arr), np.std(arr)

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"



##################################################################
# CHANNEL LOSS PLOT


file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
file_night = 'data/HORIZONTAL_DISTANCE_NIGHT2_trx=15_L=%i_wl=0.81'

# dists = [7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000]
dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]


m_d = []
std_d = []
m_n = []
std_n = []


scattering_loss = 0.45 # db/km for 810

for d in dists:
    data_day = sp.io.loadmat(file_day %d)
    data_day = data_day['res']
    data_day = -10*np.log10(data_day)+scattering_loss*d*1e-3
    # data_day = -10*np.log10(data_day) 

    m_d += [np.average(data_day)]
    std_d += [np.std(data_day)]
    
    data_n = sp.io.loadmat(file_night %d)
    data_n = data_n['res']
    data_n = -10*np.log10(data_n)+scattering_loss*d*1e-3
    # data_n = -10*np.log10(data_n)
    
    m_n += [np.average(data_n)]
    std_n += [np.std(data_n)]
    
# plt.errorbar(distsi, m_d, std_d)
# plt.errorbar(distsi, m_n, std_n)

fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.plot(dists, m_d, linewidth=3, label='Day')
plt.plot(dists, m_n, linewidth=3, label='Night')

plt.xlabel('Distance (m)',fontsize=13)
plt.ylabel('Loss (dB)', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)


plt.legend()
fig.savefig('loss.pdf', dpi=200)


############################################
# CALCULATION OF THE BACKGROUND COUNTS

# cn2_day = 2.47e-13
cn2_day = 1.52e-14
cn2_night = 1.25e-15
wl = 810e-9

det_loss = 5
fov_angle_AO = 30e-6
fov_angle_nAO = 182e-6
fov = 4*np.pi*np.sin(fov_angle_AO/2)**2

spectral_width = 2.3
# spectral_width = 0.17
sky_radiance = 0.004
a = 0.075 # aperture radious 

counts_d = spectral_width*sky_radiance*fov*np.pi*(a**2)*wl/(sp.constants.c* sp.constants.h)*10**(-det_loss/10)

background_day = counts_d
background_night = spectral_width*(sky_radiance/50)*fov*np.pi*(a**2)*wl/(sp.constants.c* sp.constants.h)*10**(-det_loss/10) 

print('Background day: ', background_day)
print('Background night: ', background_night)


#####################################################################################################



jitter = np.sqrt(2)*350e-12

time = 1

sig_rate = 20e6*10**(-det_loss/10)
# sig_rate = 500e3*10**(-det_loss/10)

qber_dev = 0.05

# loss_filter = 4.3 # 4.3 dB of extra loss due to 1 nm filter on 2.4 nm FWHM source
loss_filter = 0 

# dists = [7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000]
dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]

m_d = []
std_d = []
m_n = []
std_n = []

for d in dists:
    data_day = sp.io.loadmat(file_day %d)
    data_day = data_day['res']
    # data_day = -10*np.log10(data_day)
    
    data_day = data_day*10**(-(scattering_loss*d*1e-3 + det_loss + loss_filter)/10)
    
    m_d += [np.average(data_day)]
    std_d += [np.std(data_day)]
    
    data_night = sp.io.loadmat(file_night %d)
    data_night = data_night['res']
    
    data_night = data_night*10**(-(scattering_loss*d*1e-3 + det_loss + loss_filter)/10)
    # data_n = -10*np.log10(data_n)
    
    m_n += [np.average(data_night)]
    std_n += [np.std(data_night)]
    
    
loss_day = np.array(m_d)
loss_night = np.array(m_n)


qber_day = QBER_2sig(qber_dev, sig_rate, background_day, loss_day, jitter, time)
qber_night = QBER_2sig(qber_dev, sig_rate, background_night, loss_night, jitter, time)

fig = plt.figure()
plt.plot(dists, qber_day, 'o', label='Day', linewidth=3)
plt.plot(dists, qber_night, 'o', label='Night', linewidth=3)
# plt.xlabel('Distance (seg)', fontsize=15)
# plt.ylabel('QBER', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)
plt.ylim([0.0, 0.27])

fig.set_size_inches(18.5*.3, 10.5*.3)
fig.savefig('qber_main.pdf', dpi=200)

####################################################

def avg_std(arr):
    return np.average(arr), np.std(arr)

# .21 - 11000
# .24 - 12000

qber_real = 0.21
# qber_real = qber_day[5]


qber_limit = 0.25


time = 10

rate = int(sig_rate*loss_day[5]/2) 
# rate  = int(sig_rate*0.0011/2)

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
ax1.plot(time_x, np.ones_like(time_x)*qber_limit, 'r--', alpha = 1)    
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

# print(confidence[::10])

fig.set_size_inches(18.5*.3, 10.5*.3)
fig.savefig('qber2.pdf', dpi=200)

plt.show()    
    
############### PRECISION PLOT

alpha_precision = 2

coincidences_d = .1 * sig_rate * np.array(m_d)
coincidences_n = .1 * sig_rate * np.array(m_n)

sigs_d = jitter*alpha_precision/np.sqrt(coincidences_d)
sigs_n = jitter*alpha_precision/np.sqrt(coincidences_n)

fig = plt.figure()
plt.plot(dists, sigs_d*1e12, '-', label='Day', linewidth=3)
plt.plot(dists, sigs_n*1e12, '-', label='Night', linewidth=3)
# plt.xlabel('Distance (seg)', fontsize=15)
plt.ylabel('Precision (picoseg)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

fig.set_size_inches(18.5*.3, 10.5*.3)
plt.ylim([0.0, 70])


fig.savefig('time_precision.pdf', dpi=200)