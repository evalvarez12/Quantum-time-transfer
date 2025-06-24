# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 10:19:57 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp



#######################

def binom_prob(n, p):
    binom = np.zeros(n+1)
    for i in range(n+1):
        # array contianing the probabilities of each photon loss combination
        # from binomial expansion, starts from no photons to all photons 
        binom[i] = sp.special.binom(n, i)*((1-p)**(n-i))*(p**(i))
    return binom



#######################

ra = 20e6 # pulse rate
rb = 500e3 # background rate
jitter = 364e-12 # total time jitter
e_d = 0.05 # detector QBER


channel_loss = 12.5 # in dB
detector_loss = 5 
eta_db = channel_loss + detector_loss # total loss
eta = 10**(-eta_db/10) # loss
T = 1 # total time of prototocol



time_filter = 2
bg_prob = time_filter*jitter*rb*10**(-detector_loss/10) # time filter to select a coincidence - units of time jitter sigma
# time filter adds extra loss
if time_filter == 2:
    eta = eta*0.6826
if time_filter == 3:
    eta = eta*0.954

m = 0.4 # signal mean photon number
v = 0.2 # decoy mean photon number  - required v < m
l = 0.001 # vacuum mean photon number


pm = 0.5 # probability of sigal produced
pv = 0.5 # probability of decoy produced
pl = 1 - pm - pv # probability of vacuum

# values from optimizer
m = .5154
v = .1645
pm = .7966
pv = 1 - pm

# T = 0.1
# eta = 0.2
# m = 1.5
# v = 0.5

###### Parameter optimization
print('-------------------------- Theory')
Y0 = bg_prob # bg probability after sifting

# Number of sent pulses
# Factor of 2 dividing for basis sifting
nm_t = ra*T*pm
nv_t = ra*T*pv
# nl_t = ra*T*pl

# Number of received pulses 
Nm_t = nm_t*(bg_prob + 1 - np.exp(-eta*m))
Nv_t = nv_t*(bg_prob + 1 - np.exp(-eta*v))
# Nl_t = nl_t*(bg_prob + 1 - np.exp(-eta*l))/2


print('nm: ', nm_t)
print('nv: ', nv_t)
# print('nl: ', nl_t)

print('Nm: ', Nm_t)
print('Nv: ', Nv_t)
# print('Nl: ', Nl_t)



Qm_t = Nm_t/nm_t
Qv_t = Nv_t/nv_t
# Ql_t = Nl_t/nl_t

std_Qm = np.sqrt(Qm_t*(1- Qm_t)/nm_t)
std_Qv = np.sqrt(Qv_t*(1- Qv_t)/nv_t)
# std_Ql = np.sqrt(Ql_t *(1- Ql_t)/nl_t)


print('Qm: ', Qm_t, std_Qm)
print('Qv: ', Qv_t, std_Qv)
# print('Ql: ', Ql_t, std_Ql)

# Ma et al
qber_m_t = (Y0/2 + e_d*(1-np.exp(-eta*m)))/(1-np.exp(-eta*m) + Y0)
qber_v_t = (Y0/2 + e_d*(1-np.exp(-eta*v)))/(1-np.exp(-eta*v) + Y0)
# qber_l_t = (Y0/2 + e_d*(1-np.exp(-eta*l)))/(1-np.exp(-eta*l) + Y0)

std_qber_m_t = np.sqrt(qber_m_t*(1 - qber_m_t)/(Nm_t/2)) # Factor of 2 dividing for basis sifting
std_qber_v_t = np.sqrt(qber_v_t*(1 - qber_v_t)/(Nv_t/2))
# std_qber_l = np.sqrt(qber_l_t*(1 - qber_l_t)/Nl_t)


print('QBER m: ', qber_m_t, std_qber_m_t)
print('QBER v: ', qber_v_t, std_qber_v_t)
# print('QBER l: ', qber_l_t, std_qber_l)

Y1_t = m/(m*v-v**2) * (Qv_t*np.exp(v) - Qm_t*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * Y0)

e1_t = (Qm_t*qber_m_t*np.exp(m) - Qv_t*qber_v_t*np.exp(v))/ (Y1_t*(m-v))

e1_t2 = (Qv_t*qber_v_t*np.exp(v) - Y0/2)/(Y1_t*v)


e1_t3 = (Qm_t*qber_m_t*np.exp(m) - Y0/2)/(Y1_t*m)

print('Y1: ', Y1_t)
print('e1 - signal & decoy : ', e1_t)

print('e1 - signal : ', e1_t2)

print('e1 - decoy : ', e1_t3)


print('Num bg photons: ', bg_prob*ra*T)

print('Single photon pulses detected: ', nm_t*m*Y1_t*np.exp(-m) + nv_t*v*Y1_t*np.exp(-v))
print('Single photon signal: ',  nm_t*m*Y1_t*np.exp(-m))
print('Single photon decoy: ', nv_t*v*Y1_t*np.exp(-v))


print('---------- With central limit theorem')

prob_success_protocol = 0.95
a = sp.stats.norm.ppf(prob_success_protocol)
Qm_t_u = Qm_t + a*std_Qm
Qv_t_l = Qv_t - a*std_Qv 
Y0_u = Y0 

Y1_t_l = m/(m*v-v**2) * (Qv_t_l*np.exp(v) - Qm_t_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * Y0_u)
print('Single photon pulses detected: ', nm_t*m*Y1_t_l*np.exp(-m) + nv_t*v*Y1_t_l*np.exp(-v))

Qm_t_l = Qm_t - a*std_Qm
Qv_t_u = Qv_t + a*std_Qv 
qber_m_t_l = qber_m_t - a*std_qber_m_t
qber_m_t_u = qber_m_t + a*std_qber_m_t
qber_v_t_u = qber_v_t + a*std_qber_v_t

e1_t_u = (Qv_t_u*qber_v_t_u*np.exp(v) - Y0/2)/(Y1_t_l*v)

e1_t_u_2 = (Qm_t_u*qber_m_t_u*np.exp(m) - Y0/2)/(Y1_t_l*m)

# e1_t_u = (Qm_t_l*qber_m_t_l*np.exp(m) - Qv_t_u*qber_v_t_u*np.exp(v))/ (Y1_t_l*(m-v))
print('e1 - signal: ', e1_t_u_2)
print('e1 - decoy: ', e1_t_u)

print('Y1_l: ', Y1_t_l)

####### Numerical simulation - Statistics

print('------------- Numerical simulation')
pulses = np.random.choice([1, 2, 0], int(T*ra), p=[pm, pv, pl])

# Create initial polarization states - for simplicity one basis
pulses_pol_prepared = np.zeros_like(pulses)
pulses_pol_prepared[np.random.rand(len(pulses)) < 0.5] = 1

pulses_pol = np.copy(pulses_pol_prepared)

# Add random error
mask = np.random.rand(len(pulses)) < e_d
pulses_pol[mask] = (pulses_pol[mask] + 1) % 2


nm = sum(pulses == 1)
nv = sum(pulses == 2)
nl = sum(pulses == 0)

photons = np.copy(pulses)
photons[pulses == 1] = np.random.poisson(m, nm)
photons[pulses == 2] = np.random.poisson(v, nv)
photons[pulses == 0] = np.random.poisson(l, nl)

og_photons = np.copy(photons)


max_ph = np.max(og_photons)

for i in range(1,max_ph):
    l = sum(og_photons == i)
    
    probs = binom_prob(i, eta)
    photons[og_photons == i] = np.random.choice(np.arange(i+1), p = probs, size=l)

# # independent loss of photons
# for i in range(len(photons)):
#     probs = binom_prob(photons[i], eta)
    
#     photons[i] = np.random.choice(np.arange(photons[i]+1), p = probs)


bg = (np.random.rand(int(T*ra)) < bg_prob).astype(int)

pulses_pol[bg.astype(bool)] = np.random.choice([0,1], np.sum(bg.astype(bool)))


photons = photons + bg

clicks = photons != 0

# Number of sent pulses
nm = len(pulses[pulses == 1])
nv = len(pulses[pulses == 2])
# nl = len(pulses[pulses == 0])

# Number of received pulses
Nm = sum(clicks[pulses == 1])
Nv = sum(clicks[pulses == 2])
# Nl = sum(clicks[pulses == 0])

Qm = Nm/nm
Qv = Nv/nv
# Ql = Nl/nl

print('Nm: ', Nm)
print('Nv: ', Nv)
# print('Nl: ', Nl)

print('Qm: ', Qm)
print('Qv: ', Qv)
# print('Ql: ', Ql)

inds_clicks_m = np.logical_and(clicks == 1, pulses == 1)
inds_clicks_v = np.logical_and(clicks == 1, pulses == 2)
# inds_clicks_l = np.logical_and(clicks == 1, pulses == 0)


qber_m = np.sum(pulses_pol_prepared[inds_clicks_m] != pulses_pol[inds_clicks_m])/Nm
qber_v = np.sum(pulses_pol_prepared[inds_clicks_v] != pulses_pol[inds_clicks_v])/Nv
# qber_l = np.sum(pulses_pol_prepared[inds_clicks_l] != pulses_pol[inds_clicks_l])/Nl

std_qber_m = np.std(pulses_pol_prepared[inds_clicks_m] != pulses_pol[inds_clicks_m])/np.sqrt(Nm)
std_qber_v = np.std(pulses_pol_prepared[inds_clicks_v] != pulses_pol[inds_clicks_v])/np.sqrt(Nv)

print('QBER m: ', qber_m, std_qber_m)
print('QBER v: ', qber_v , std_qber_v)
# print('QBER l: ', qber_l)


Y1_bound = m/(m*v-v**2) * (Qv*np.exp(v) - Qm*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * Y0)
e1_bound = (Qv*qber_v*np.exp(v) - Y0/2)/(Y1_bound*v)

e1_bound2 = (Qm*qber_m*np.exp(m) - Y0/2)/(Y1_bound*m)

print('Y1 bound: ', Y1_bound)

print('e1 - signal & decoy: ', (Qm*qber_m*np.exp(m) - Qv*qber_v*np.exp(v))/ (Y1_bound*(m-v)))
print('e1 - signal: ', e1_bound2)
print('e1 - decoy: ', e1_bound)

print('Bound single photon pulses detected: ', nm*m*Y1_bound*np.exp(-m) + nv*v*Y1_bound*np.exp(-v))

print('------ 1 sigma')
qber_v_1s = qber_v + std_qber_v
e1_bound_1s = (Qv*qber_v_1s*np.exp(v) - Y0/2)/(Y1_bound*v)
print('e1 bound: ', e1_bound_1s)

print('------ 2 sigma')
qber_v_2s = qber_v + 2*std_qber_v
e1_bound_2s = (Qv*qber_v_2s*np.exp(v) - Y0/2)/(Y1_bound*v)
print('e1 bound: ', e1_bound_2s)

print('------ 3 sigma')
qber_v_3s = qber_v + 3*std_qber_v
e1_bound_3s = (Qv*qber_v_3s*np.exp(v) - Y0/2)/(Y1_bound*v)
print('e1 bound: ', e1_bound_3s)


print('------ Real')

clicks_1 = np.logical_and(og_photons == 1, clicks)
ns = np.logical_and(clicks_1, bg == 0)

qber1 = np.sum(pulses_pol_prepared[clicks_1] != pulses_pol[clicks_1])/np.sum(clicks_1)
print('e1: ', qber1)

print('Num bg photons: ', np.sum(bg))


print('Single photon pulses detected: ', np.sum(clicks_1))
print('Num true single photons: ', np.sum(ns))
print('Num bg passing as single photons: ', np.sum(np.logical_and(clicks_1, bg == 1)))


print('Single photon signal: ', np.sum(np.logical_and(ns, pulses == 1)))
print('Single photon decoy: ', np.sum(np.logical_and(ns, pulses == 2)))



print('----------------- Comparing theory simulation')
print('-- Theory')
print('Y1: ', Y1_t)
print('e1 - signal & decoy : ', e1_t)

print('e1 - signal : ', e1_t2)

print('e1 - decoy : ', e1_t3)

print('-- Simulation')
print('Y1 bound: ', Y1_bound)

print('e1 - signal & decoy: ', (Qm*qber_m*np.exp(m) - Qv*qber_v*np.exp(v))/ (Y1_bound*(m-v)))
print('e1 - signal: ', e1_bound2)
print('e1 - decoy: ', e1_bound)