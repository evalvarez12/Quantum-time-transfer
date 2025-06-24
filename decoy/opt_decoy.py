# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
import scipy.special as sp

def func_Y1(m, v, Qm, Qv):
    return m/(m*v-v**2) * (Qv*np.exp(v) - Qm*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * BG_PROB)

def func_e1(m, v, Qm, Qv, qber_m, qber_v):
    Y1 = func_Y1(m, v, Qm, Qv)
    return  (Qm*qber_m*np.exp(m) - Qv*qber_v*np.exp(v))/ (Y1*(m-v))



def ser(m, v):
    N = 1e5
    Y = np.random.rand(int(N))
    
    res1 = 0
    for i in range(2, len(Y)):
        res1 += Y[i]*(m**i)/sp.factorial(i)
        
    res2 = 0
    for i in range(2, len(Y)):
        res2 += Y[i]*(v**i)/sp.factorial(i)
    return res2 < (v/m)**2*res1

        
########################### PARAMETERS


ra = 20e6 # pulse rate
rb = 500e3 # background rate
jitter = 500e-12 # total time jitter
e_d = 0.05 # detector QBER


channel_loss = 25 # in dB
detector_loss = 5 
eta_db = channel_loss + detector_loss # total loss
eta = 10**(-eta_db/10) # loss
t = .1 # total time of prototocol

time_filter = 2
bg_prob = time_filter*jitter*rb*10**(-detector_loss/10) # time filter to select a coincidence - units of time jitter sigma
# time filter adds extra loss
if time_filter == 2:
    eta = eta*0.6826
if time_filter == 4:
    eta = eta*0.954

m = 0.5 # signal mean photon number
v = 0.3 # decoy mean photon number  - required v < m
l = 0.001 # vacuum mean photon number


pm = 0.6 # probability of sigal produced
pv = 0.4 # probability of decoy produced



#### GLOBAL VARS
ETA = eta
TIME = t
BG_PROB = bg_prob
RA = ra


print('-------------------------- Theory')
Y0 = bg_prob

# Number of sent pulses
nm_t = ra*t*pm
nv_t = ra*t*pv


# Number of received pulses
Nm_t = nm_t*(bg_prob + 1 - np.exp(-eta*m))
Nv_t = nv_t*(bg_prob + 1 - np.exp(-eta*v))



print('nm: ', nm_t)
print('nv: ', nv_t)


print('Nm: ', Nm_t)
print('Nv: ', Nv_t)


print('Num bg photons: ', bg_prob*ra*t)

Qm_t = Nm_t/nm_t
Qv_t = Nv_t/nv_t


print('Qm: ', Qm_t)
print('Qv: ', Qv_t)


# Ma et al
qber_m_t = (Y0/2 + e_d*(1-np.exp(-eta*m)))/(1-np.exp(-eta*m) + Y0)
qber_v_t = (Y0/2 + e_d*(1-np.exp(-eta*v)))/(1-np.exp(-eta*v) + Y0)

print('QBER m: ', qber_m_t)
print('QBER v: ', qber_v_t)

Y1_t = m/(m*v-v**2) * (Qv_t*np.exp(v) - Qm_t*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * Y0)

e1_t = (Qm_t*qber_m_t*np.exp(m) - Qv_t*qber_v_t*np.exp(v))/ (Y1_t*(m-v))

print('Y0: ', Y0)
print('Q1: ', np.exp(-m)*Y1_t)
print('e1: ', e1_t)

print('Single photon pulses detected: ', nm_t*Y1_t*np.exp(-m)*m)
print('Single photon pulses detected: ', (nm_t*np.exp(-m)*m + nv_t*np.exp(-v)*v)*Y1_t)








########################### OPTIMIZATION
print('----------------- Optimization')

def cost_Y1(x):
    m, v, pm = x
    pv = 1 - pm
    
    Qm = (BG_PROB + 1 - np.exp(-ETA*m))
    Qv = (BG_PROB + 1 - np.exp(-ETA*v))
    
    nm = RA*TIME*pm
    nv = RA*TIME*pv

    Y1 = func_Y1(m, v, Qm, Qv)
    
    n = (nm*np.exp(-m)*m + nv*np.exp(-v)*v)*Y1
    
    qber_m = (Y0/2 + e_d*(1-np.exp(-eta*m)))/(1-np.exp(-eta*m) + Y0)
    qber_v = (Y0/2 + e_d*(1-np.exp(-eta*v)))/(1-np.exp(-eta*v) + Y0)

    e1 = func_e1(m, v, Qm, Qv, qber_m, qber_v)
    
    cost = - n
    
    if e1 > .2:
        order_n = np.floor(np.log10(n)) + 1
        cost = e1*10**(order_n) - n 
    

    return cost



x0 = [m, v, pm]
cons = ({'type': 'ineq', "fun": lambda x: x[0] - x[1] - 0.1},
        {'type': 'ineq', "fun": lambda x: x[0] - 0.0001},
        {'type': 'ineq', "fun": lambda x: x[1] - 0.0001},
        {'type': 'ineq', "fun": lambda x: 1 - x[0]},
        {'type': 'ineq', "fun": lambda x: 1 - x[1]},
        {'type': 'ineq', "fun": lambda x: x[2] - 0.001},
        {'type': 'ineq', "fun": lambda x: 1 - x[2]})


# x0 = [0.5, 0.4]
# cons = ({'type': 'ineq', "fun": lambda x: x[0] - x[1]},
#         {'type': 'ineq', "fun": lambda x: x[0]},
#         {'type': 'ineq', "fun": lambda x: x[1]},
#         {'type': 'ineq', "fun": lambda x: 1 - x[0]},
#         {'type': 'ineq', "fun": lambda x: 1 - x[1]})

# cons = ({'type': 'ineq', "fun": lambda x: 2 - x[0]},
#         {'type': 'ineq', "fun": lambda x: x[0] - 0.0001},
#         {'type': 'ineq', "fun": lambda x: x[1] - 0.0001},
#         {'type': 'ineq', "fun": lambda x: x[0] - x[1]})

res = opt.minimize(cost_Y1, x0, constraints=cons)

print(res)


m, v, pm = res['x']
pv = 1 - pm

Qm = (BG_PROB + 1 - np.exp(-ETA*m))
Qv = (BG_PROB + 1 - np.exp(-ETA*v))

qber_m = (Y0/2 + e_d*(1-np.exp(-eta*m)))/(1-np.exp(-eta*m) + Y0)
qber_v = (Y0/2 + e_d*(1-np.exp(-eta*v)))/(1-np.exp(-eta*v) + Y0)

Y1 = func_Y1(m, v, Qm, Qv)
e1 = func_e1(m, v, Qm, Qv, qber_m, qber_v)

print('Q1: ', np.exp(-m)*Y1)
print('e1: ', e1)

print('QBER m: ', qber_m)
print('QBER v: ', qber_v)

print('Single photon pulses: ', np.exp(-m)*Y1*m*ra*t)
print('Single photon pulses: ', (RA*TIME*pm*np.exp(-m)*m + RA*TIME*pv*np.exp(-v)*v)*Y1)

# m, v, pv, pm = res['x']
# print(cost_Y1(res['x']))
####################################################################################

# LOW_COST = 10.
# MID_COST = 150.
# HIGH_COST = 400.

# def weight(a, b, c, d):
#     return "calculated weight of structure"

# def frequency(a, b, c, d):
#     return "calculated resonant frequency"

# def freq_penalty(freq):
#     # Example linear piecewise penalty function -
#     #   increasing cost for frequencies below 205 or above 395
#     if freq < 205:
#         return MID_COST * (205 - freq)
#     elif freq < 395:
#         return 0.
#     else:
#         return MID_COST * (freq - 395)

# def stress_fraction(a, b, c, d):
#     return "calculated stress / failure criteria"

# def stress_penalty(stress_frac):
#     # Example linear piecewise penalty function -
#     #   low extra cost for stress fraction below 0.85,
#     #   high extra cost for stress fraction over 0.98
#     if stress_frac < 0.85:
#         return LOW_COST * (0.85 - stress_frac)
#     elif stress_frac < 0.98:
#         return 0.
#     else:
#         return HIGH_COST * (stress_frac - 0.98)

# def overall_fitness(parameter_vector):
#     a, b, c, d = parameter_vector
#     return (
#         # D'oh! it took me a while to get this right -
#         # we want _minimum_ weight and _minimum_ penalty
#         # to get _maximum_ fitness.
#        -weight(a, b, c, d)
#       - freq_penalty(frequency(a, b, c, d))
#       - stress_penalty(stress_fraction(a, b, c, d)
#     )


# from scipy.optimize import fmin

# initial_guess = [29., 45., 8., 0.06]
# result = fmin(lambda x: -overall_fitness(x), initial_guess, maxfun=100000, full_output=True)

