# -*- coding: utf-8 -*-

import numpy as np


def QBER_2sig(ed, ra, rb, eta, sig, T):
    time_filter = 2*sig
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*ra*T
    Nb = rb*ra*T*time_filter
    
    return (Ns*ed + Nb/2)/(Ns + Nb)


source_rate = 20e6

source = 'entangled'
# source = 'bb84'

det_loss = 5 # detector loss - dB
channel_loss = 12.5 # in dB
time_QBER = 1 # integration time QBER estimation
time = .1 # integration time PNT
ed = 0.05 # device QBER contribution 

if source == 'entangled':
    eta_db = channel_loss + 2*det_loss # total loss
    total_jitter = 500e-12 # entangled source
elif source == 'bb84':
    eta_db = channel_loss + det_loss # total loss
    total_jitter = 364e-12 # BB84 source



eta = 10**(-eta_db/10) # loss





# Background counts calculation 

hc = 1.98644586e-25
w = 24e-6 # half angle
fov = 4*np.pi*np.sin((w/2))**2 # half angle field of view
area = np.pi * (.15/2)**2  # area of aperture - 150 mm
L = 0.004 # sky radiance
filt = 2.4

lamb = 810e-9

energy = L*fov*area*filt

N = energy * lamb/ hc

detector_eff = 10**(-det_loss/10)

bg_rate = N * detector_eff

print('Background rate: ', bg_rate)
bg_counts = time*source_rate* bg_rate * 2 * total_jitter
print('Background counts: ', bg_counts)
QBER = QBER_2sig(ed, source_rate, bg_rate, eta, total_jitter, time_QBER)

print('Asymtotic QBER: ', QBER)



# Signal calculation
coincidences = source_rate*eta*time + time*source_rate* bg_rate * 2 * total_jitter
# coincidences = source_rate*eta*time
print('Coincidences: ', coincidences)
precision = 2*total_jitter/np.sqrt(coincidences/2) # factor of 2 after sifting
print('Precision (picosec): ', precision*1e12)

