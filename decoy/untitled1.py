# -*- coding: utf-8 -*-

import numpy as np 


class CostFunc:
    def __init__(self, n, eta, bg_prob, e_d confidence):
        self.n = n
        self.eta = eta
        self.bg_prob = bg_prob
        self.e_d = e_d
        self.e = (1 - confidence)/3
        
    
    def lower_bound(self, x, n, e, e1):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        if (2/e1)**(1/ml) > np.exp(3/(4*np.sqrt(2))**2):
            print('Condition not satisfied lower bound : ', x)
            # return 0
        
        D = np.sqrt(2*x*np.log(1/(e1**4/16)))
        return x - D
        
    def upper_bound(self, x, n, e, e2):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        if (1/e2)**(1/ml) > np.exp(1/3):
            print('Condition not satisfied upper bound : ', x)
            # return 0
        
        D = np.sqrt(2*x*np.log(1/(e2**(3/2))))
        return x + D
    
    
    
    def vals(self, x):
        m, v, pm = x
        pv = 1 - pm
        
        # Number of sent pulses
        nm = self.n*pm
        nv = self.n*pv

        # Number of received pulses
        Nm = nm*(1 - np.exp(-self.eta*m))
        Nv = nv*(1 - np.exp(-self.eta*v))
        
        Qm_u = self.upper_bound(Nm, nm, e, e)/nm
        Qv_u = self.upper_bound(Nv, nv, e, e)/nv  
        Qv_l = self.lower_bound(Nv, nv, e, e)/nv 
        
        # Number of received background
        Nb = self.n*self.bg_prob
        
        qber_m = (Nb*pm/2 + self.e_d*Nm)/(Nb*pm + Nm)
        qber_v = (Nb*pv/2 + self.e_d*Nv)/(Nb*pv + Nv)
        
        # Number of coincidences divided by two after sifting
        qber_m_u = self.upper_bound(Nm*qber_m/2, Nm/2, e, e)/(Nm/2)
        qber_v_u = self.upper_bound(Nv*qber_v/2, Nv/2, e, e)/(Nv/2)
        
        
        Y1_l = m/(m*v-v**2) * (Qv_l*np.exp(v) - Qm_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * self.bg_prob)
        ns = self.n*Y1_l*(pm*np.exp(-m)*m + pv*np.exp(-v)*v)
        
        e1_u_1 = (Qv_u*qber_v_u*np.exp(v) - self.bg_prob/2)/(Y1_l*v)
        e1_u_2 = (Qm_u*qber_m_u*np.exp(m) - self.bg_prob/2)/(Y1_l*m)
        
        # e1_u = min(e1_u_1, e1_u_2)
        # Decoy pulses give better QBER bounds
        e1_u = e1_u_1

        # print('QBERS:',e1_u_1, e1_u_2)
        # print('Nm: ', Nm)

        sig = Nv + Nm

        return sig, ns, e1_u
        
        
    
    def __call__(self, x):
        s, n, e = self.vals(x)
        return e
        
        

