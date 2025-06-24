# -*- coding: utf-8 -*-

import numpy as np 
import scipy as sp


class CostFunc:
    def __init__(self, n, eta, Y0, e_d, confidence, emin=0, ecap=0):
        self.n = n
        self.eta = eta
        self.Y0 = Y0
        self.e_d = e_d
        self.e = (1 - confidence)/6
        self.consts = [None]*3
        
        # Used for second optimization
        self.emin = emin
        self.ecap = ecap
        self.conste = None
        
        
    def lower_bound(self, x, n, e, e1, ind):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        self.consts[ind] = (-(2/e1)**(1/ml) + np.exp(3/(4*np.sqrt(2))**2))
        # if (2/e1)**(1/ml) > np.exp(3/(4*np.sqrt(2))**2):
            # print('Condition not satisfied lower bound : ', x)
            # return 0
        
        D = np.sqrt(2*x*np.log(1/(e1**4/16)))
        return x - D
        
    
    def upper_bound(self, x, n, e, e2, ind):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        self.consts[ind] = (-(1/e2)**(1/ml) + np.exp(1/3))
        # if (1/e2)**(1/ml) > np.exp(1/3):
            # print('Condition not satisfied upper bound : ', x)
            # return 0
        
        D = np.sqrt(2*x*np.log(1/(e2**(3/2))))
        return x + D
    
    
    
    def vals(self, x):
        m, v, pm = x
        pv = 1 - pm
        
        a = sp.stats.norm.ppf(1-self.e)
        
        # Number of sent pulses
        nm = self.n*pm
        nv = self.n*pv

        # Number of received pulses
        Nm = nm*(1 - np.exp(-self.eta*m))
        Nv = nv*(1 - np.exp(-self.eta*v))
        
        Qm = Nm/nm
        Qv = Nv/nv
        
        Qm_u = self.upper_bound(Nm, nm, self.e, self.e, 0)/nm
        Qv_u = self.upper_bound(Nv, nv, self.e, self.e, 1)/nv  
        Qv_l = self.lower_bound(Nv, nv, self.e, self.e, 2)/nv 
        
        # Number of received background
        Nb = self.n*self.Y0
        
        qber_m = (Nb*pm/2 + self.e_d*Nm)/(Nb*pm + Nm)
        qber_v = (Nb*pv/2 + self.e_d*Nv)/(Nb*pv + Nv)
        
        # # Number of coincidences divided by two after sifting
        # qber_m_u = self.upper_bound(Nm*qber_m/2, Nm/2, self.e, self.e, 3)/(Nm/2)
        # qber_v_u = self.upper_bound(Nv*qber_v/2, Nv/2, self.e, self.e, 4)/(Nv/2)
        
        # Number of coincidences divided by two after sifting
        std_qber_m = np.sqrt(qber_m*(1- qber_m)/(Nm/2))
        std_qber_v = np.sqrt(qber_v*(1- qber_v)/(Nv/2))
        
        
        qber_m_u = qber_m + a*std_qber_m
        qber_v_u = qber_v + a*std_qber_v
        
        
        Y1_l = m/(m*v-v**2) * (Qv_l*np.exp(v) - Qm_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * self.Y0)
        ns = self.n*Y1_l*(pm*np.exp(-m)*m + pv*np.exp(-v)*v)
        
        e1_u_1 = (Qv_u*qber_v_u*np.exp(v) - self.Y0/2)/(Y1_l*v)
        e1_u_2 = (Qm_u*qber_m_u*np.exp(m) - self.Y0/2)/(Y1_l*m)
        
        # e1_u = min(e1_u_1, e1_u_2)
        # Decoy pulses give better QBER bounds
        e1_u = e1_u_1

        # print('QBERS:',e1_u_1, e1_u_2)
        # print('Nm: ', Nm)
        
        
        b1 = ((Qm*np.exp(m) - self.Y0)/(2*m) < Y1_l)
        b2 = ((Qv*np.exp(v) - self.Y0)/(2*v) < Y1_l)
        # print(b1, b2)

        sig = Nm*b1 + Nv*b2
        
        # self.consts[5] = (ns - sig/2)
        self.conste = self.emin + self.ecap - e1_u 

        return sig, Y1_l, e1_u
        
        
    
    # def __call__(self, x):
    #     s, n, e = self.vals(x)
    #     # print('--- ', x, s)
    #     return -1000*s
    
    
    def costfunc1(self, x):
        s, n, e = self.vals(x)
        # print('--- ', x, s)
        return -1000*n
        
    
    def costfunc2(self,x):
        sig, n, e = self.vals(x)
        # print('--- ', x, sig, self.consts, self.conste)
        return -100*sig
    
    
    def const_erange(self, x):
        res = self.conste
        if res == None:
            _ = self.vals(x)
            res = self.conste
        
        self.conste = None
        return res
    
        
    def consti(self, x, i):
        res = self.consts[i]
        if res == None:
            _ = self.vals(x)
            res = self.consts[i]
        
        self.consts[i] = None
       
        return res
    
    
    def const1(self, x):
        return self.consti(x, 0)
    
    def const2(self, x):
        return self.consti(x, 1)

    def const3(self, x):
        return self.consti(x, 2)
    
    def const4(self, x):
        return self.consti(x, 3)
    
    def const5(self, x):
        return self.consti(x, 4)
        
    def const6(self, x):
        return self.consti(x, 5)
