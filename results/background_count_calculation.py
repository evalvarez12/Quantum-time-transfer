# -*- coding: utf-8 -*-

"""
Created on Thu Jan  4 11:17:38 2024

@author: vil034
"""

import numpy as np
import scipy as sp

c = 3e8
h = 6.62607015e-34


time = 1
apperture = 0.15
w = 100e-6 # half angle of the cone that creates the solid angle
fov = 4*np.pi*np.sin((w/2))**2 # field of view - solid angle
area = np.pi * (apperture/2)**2  # area of aperture - 150 mm
# L = 0.004 # sky radiance
L = 0.004/50

filt = 3

lamb = 810e-9

energy = L*fov*area*filt

N = energy * lamb/ (h*c)

detector_eff = 10**(-6/10)

print('Background counts')
print(N*detector_eff)


# ################# Second calculations with smaller aperture
# area = np.pi * (0.044)**2
# energy = L*fov*area*filt
# N = energy * lamb/ hc
# print()
# # print('Background counts')
# # print(N*detector_eff)

# N = energy * lamb/ hc




# ############### Theory of fail probability

# background_rate = 37000 # background counts per second 
# qber = 0.05 # device qber not considering backgroud radiation
# loss_channel = 20 # loss of the channel - only afects Bob
# loss_detector = 5 # loss of the detection system - affects Alice and Bob
# signal_rate = 200000 # average entangled pair generation - 20 MHz
# jitter = 364e-12
# # jitter = 10e-6 
# total_protocol_time = .1

# print('-------------------- Theory')
# # Final step compare with theory to se if results are similar to what is expected
# # Results might not be perfect as theory does not accounts for small effects such as detector dead time

# # each bin has the width of time jitter - jitter is the time unit
# total_jitter = np.sqrt(2)*jitter
# N = total_protocol_time/(total_jitter)
# dt = total_protocol_time/N

# Ns = signal_rate*total_protocol_time*(10**(-(loss_channel + 2*loss_detector)/10)) # number of coincidences
# ra = signal_rate*total_protocol_time*(10**(-loss_detector/10))/N # number of Alice data points
# rb = background_rate*total_protocol_time/N # Number of background photons for Bob
# p = ra*rb
# mu_b = N*p
# # sig_b2 = N*p*(1-p)
# # print('Avg noise: ', mu_b)
# # print('Std noise: ', sig_b2)

# # Probability of misidentifiying the peak due to low SNR
# P_failure = .5*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
# print('Probability of large bin: ', P_failure)
# print('Probability of failure: ', P_failure*N)    

# # QBER estimation based on statistics
# print('---> 4 sig time filter')
# print('Num of coincidences: ', Ns)
# Nb = ra*N*background_rate*4*total_jitter 
# print('Num background photons misidentified as signal: ', Nb)
# qber_est = (Ns*qber + Nb/2)/(Ns + Nb)
# print('QBER: ', qber_est)

# print('---> 2 sig time filter')
# Ns2  = 0.68*Ns
# print('Num of coincidences: ', Ns2)
# Nb2 = Nb/2
# print('Num background photons misidentified as signal: ', Nb2)
# print('QBER - 2 sig: ', (Ns2*qber + Nb2/2)/(Ns2 + Nb2))
# # Estimated QBER always a bit lower than measured one for some reason - needs a closer look



