# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 16:01:49 2023

@author: vil034
"""

import numpy as np


def e1(v, m, eta, Y0, Ev):
    Qv = Y0 + 1 - np.exp(-eta*v)
    Qm = Y0 + 1 - np.exp(-eta*m)
    
    Y1 = (m/(m*v -v**2)) * (Qv*np.exp(v)- Qm*np.exp(m)*(v/m)**2 - ((m**2-v**2)/m**2) * Y0)
    # print(Y1)
    return (Ev*Qv*np.exp(v) - Y0/2) / (Y1*v)

def QBER_2sig(qd, ra, rb, eta, sig, T):
    time_filter = 2*sig
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*ra*T
    Nb = rb*ra*T*time_filter
    
    return (Ns*qd + Nb/2)/(Ns + Nb)


loss = 26
eta = 10**(-loss/10)
qber_dev = 0.01
signal_rate = 20e6 # 3 MHz 
T = .1 # in seconds
sig = 500e-12
# sig = 2e-9
background_rate2 = 500000

Y0 = 2*sig*background_rate2
v = 0.5 # decoy mean photon number
m = 1 # signal mean photon number


# print(e1(v,m,eta, Y0, Ev))


em = np.exp(m)
ev = np.exp(v)

Qv = Y0 + 1 - np.exp(-eta*v)
Qm = Y0 + 1 - np.exp(-eta*m)

# Ev = 0.20
T = 0.1
Ev = QBER_2sig(qber_dev, signal_rate, background_rate2, eta, sig, T)


Y1 = (m/(m*v -v**2)) * (Qv*np.exp(v)- Qm*np.exp(m)*(v/m)**2 - ((m**2-v**2)/m**2) * Y0)
Q1 = m**2/((m*v-v**2)*em) * (Qv*ev - Qm*em * (v/m)**2 - (m**2 - v**2)/m**2 * Y0)

print('probability of single photon pulse:',  Y1)
print(Q1*signal_rate*T)
# print((Ev*Qv*np.exp(v) - Y0/2) / (Y1*v))

print('Upper bound to single photon pulse qber: ', e1(v,m, eta, Y0, Ev))