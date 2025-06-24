# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize as opt
import scipy.special as sp

class Optimizer:
    
    def __init__(self, ra, eta, bg_rate, jitter, time_filter, time, e_d):
        self.ra = ra
        
        self.time = time
        self.e_d = e_d
        
        # time filter to select a coincidence - units of time jitter sigma
        # time filter adds extra loss
        bg_prob = time_filter*jitter*rb*10**(-detector_loss/10) 
        if time_filter == 2:
            eta = eta*0.6826
        if time_filter == 3:
            eta = eta*0.954        

        self.eta = eta
        self.bg_prob = bg_prob
        
        
    def cost_Y1(self, x):
        m, v, pm = x
        pv = 1 - pm
                
        Qm = (self.bg_prob + 1 - np.exp(-self.eta*m))
        Qv = (self.bg_prob + 1 - np.exp(-self.eta*v))
        
        nm = self.ra*self.time*pm
        nv = self.ra*self.time*pv

        Y1 = self.func_Y1(m, v, Qm, Qv)
        
        n = (nm*np.exp(-m)*m + nv*np.exp(-v)*v)*Y1
        
        qber_m = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*m)))/(1-np.exp(-self.eta*m) + self.bg_prob)
        qber_v = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*v)))/(1-np.exp(-self.eta*v) + self.bg_prob)

        e1 = self.func_e1(m, v, Qm, Qv, qber_m, qber_v)
        
        cost = - n
        
        if e1 > .2:
            order_n = np.floor(np.log10(n)) + 1
            cost = e1*10**(order_n) - n 
        

        return cost
    
    def func_Y1(self, m, v, Qm, Qv):
        return m/(m*v-v**2) * (Qv*np.exp(v) - Qm*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * self.bg_prob)

    def func_e1(self, m, v, Qm, Qv, qber_m, qber_v):
        Y1 = self.func_Y1(m, v, Qm, Qv)
        return  (Qm*qber_m*np.exp(m) - Qv*qber_v*np.exp(v))/ (Y1*(m-v))

    
    def theoretical_vals(self, confidence):
        m, v, pm = self.x
        pv = 1 - pm
        
        a = sp.stats.norm.ppf(confidence)
        
        # Number of sent pulses
        nm = self.ra*self.t*pm
        nv = self.ra*self.t*pv


        # Number of received pulses
        Nm = nm*(self.bg_prob + 1 - np.exp(-self.eta*m))
        Nv = nv*(self.bg_prob + 1 - np.exp(-self.eta*v))
        
        Qm = Nm/nm
        Qv = Nv/nv
        
        std_Qm = np.sqrt(Qm*(1- Qm)/nm)
        std_Qv = np.sqrt(Qv*(1- Qv)/nv)
        
        qber_m = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*m)))/(1-np.exp(-self.eta*m) + self.bg_prob)
        qber_v = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*v)))/(1-np.exp(-self.eta*v) + self.bg_prob)
        

        Qm_u = Qm + a*std_Qm
        Qv_l = Qv - a*std_Qv 
        
        
        Y1_l = m/(m*v-v**2) * (Qv_l*np.exp(v) - Qm_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * self.bg_prob)
        
        e1_u = (Qv_l*qber_v*np.exp(v) - self.bg_prob/2)/(Y1_l*v)

        return Y1_l, e1_u

    def set_initial_conditions(self, x0):
        self.x = x0

    def optimize(self, confidence):
        
        x0 = self.x0
        cons = ({'type': 'ineq', "fun": lambda x: x[0] - x[1] - 0.1},
                {'type': 'ineq', "fun": lambda x: x[0] - 0.0001},
                {'type': 'ineq', "fun": lambda x: x[1] - 0.0001},
                {'type': 'ineq', "fun": lambda x: 1 - x[0]},
                {'type': 'ineq', "fun": lambda x: 1 - x[1]},
                {'type': 'ineq', "fun": lambda x: x[2] - 0.001},
                {'type': 'ineq', "fun": lambda x: 1 - x[2]})

        res = opt.minimize(self.cost_Y1, x0, constraints=cons)

        self.x = res['x']
        print(res)

    

        # m, v, pm = res['x']
        # pv = 1 - pm

        # Qm = (BG_PROB + 1 - np.exp(-ETA*m))
        # Qv = (BG_PROB + 1 - np.exp(-ETA*v))

        # qber_m = (Y0/2 + e_d*(1-np.exp(-eta*m)))/(1-np.exp(-eta*m) + Y0)
        # qber_v = (Y0/2 + e_d*(1-np.exp(-eta*v)))/(1-np.exp(-eta*v) + Y0)

        # Y1 = func_Y1(m, v, Qm, Qv)
        # e1 = func_e1(m, v, Qm, Qv, qber_m, qber_v)

        # print('Q1: ', np.exp(-m)*Y1)
        # print('e1: ', e1)

        # print('QBER m: ', qber_m)
        # print('QBER v: ', qber_v)

        # print('Single photon pulses: ', np.exp(-m)*Y1*m*ra*t)
        # print('Single photon pulses: ', (RA*TIME*pm*np.exp(-m)*m + RA*TIME*pv*np.exp(-v)*v)*Y1)

        