# -*- coding: utf-8 -*-


import numpy as np
import scipy as sp

ra = 20e6 # pulse rate
bg_rate = 500e3 # background rate
jitter = 364e-12 # total time jitter
e_d = 0.05 # detector QBER
time_filter = 2 # 2 sigma time filter on the coincidence window

channel_loss = 12.5 # in dB
detector_loss = 5 
eta_db = channel_loss + detector_loss # total loss
eta = 10**(-eta_db/10) # loss
time = 1 # total time of prototocol

eta = eta*0.6826
eta_es = 10**(- (channel_loss + 2*detector_loss)/10)*0.6826


bg_prob = time_filter*jitter*bg_rate*10**(-detector_loss/10) 

conf = 0.99

def restriction(x):
    s = signal(x)
    n1, _ = theoretical_vals(x, conf)
    
    # print('--->', s, n1)
    
    return (2*n1 - s)

    
def signal(x):
    m, v, pm = x
    pv = 1 - pm
    s = ra*time*eta*(pm*(1-np.exp(-m)) + pv*(1-np.exp(-v)))
    return s


def theoretical_vals(x, confidence):
    m, v, pm = x
    pv = 1 - pm
    
    a = sp.stats.norm.ppf(confidence)
    
    # Number of sent pulses
    nm = ra*time*pm
    nv = ra*time*pv


    # Number of received pulses
    Nm = nm*(bg_prob + 1 - np.exp(-eta*m))
    Nv = nv*(bg_prob + 1 - np.exp(-eta*v))
    
    Qm = Nm/nm
    Qv = Nv/nv
    
    std_Qm = np.sqrt(Qm*(1- Qm)/nm)
    std_Qv = np.sqrt(Qv*(1- Qv)/nv)
    
    qber_m = (bg_prob/2 + e_d*(1-np.exp(-eta*m)))/(1-np.exp(-eta*m) + bg_prob)
    qber_v = (bg_prob/2 + e_d*(1-np.exp(-eta*v)))/(1-np.exp(-eta*v) + bg_prob)
    
    # Number of coincidences divided by two after sifting
    std_qber_m = np.sqrt(qber_m*(1- qber_m)/(Nm/2))
    std_qber_v = np.sqrt(qber_v*(1- qber_v)/(Nv/2))
    
    
    qber_m_u = qber_m + a*std_qber_m
    qber_v_u = qber_v + a*std_qber_v
    
    Qm_u = Qm + a*std_Qm
    Qv_l = Qv - a*std_Qv 
    Qv_u = Qv + a*std_Qv 
    
    Y1_l = m/(m*v-v**2) * (Qv_l*np.exp(v) - Qm_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * bg_prob)
    ns = ra*time*Y1_l*(pm*np.exp(-m)*m + pv*np.exp(-v)*v)
    
    e1_u_1 = (Qv_u*qber_v_u*np.exp(v) - bg_prob/2)/(Y1_l*v)
    e1_u_2 = (Qm_u*qber_m_u*np.exp(m) - bg_prob/2)/(Y1_l*m)
    
    # e1_u = min(e1_u_1, e1_u_2)
    # Decoy pulses give better QBER bounds
    e1_u = e1_u_1

    # print('QBERS:',e1_u_1, e1_u_2)
    # print('Nm: ', Nm)

    return ns, e1_u




mr = 1
vr = 0.4
pmr = 0.5


def pfun(x):
    xv = [x, vr, pmr]
    return restriction(xv)



roots = sp.optimize.fsolve(pfun, mr)
print('root: ', roots)

print(restriction([roots[0], vr, pmr]))