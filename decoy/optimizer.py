# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize as opt

class Optimizer:
    
    def __init__(self, channel):
        self.set_channel(channel)
        
        
    def set_channel(self, channel):
        # Divide by 2 due to basis sifting
        self.ra = channel['ra']
        
        self.time = channel['time']
        self.e_d = channel['e_d']
        
        self.det_loss = channel['detector_loss']
        self.channel_loss = channel['channel_loss']
        self.time_filter = channel['time_filter']
        self.jitter = channel['jitter']
        
        
        eta_db = channel['channel_loss'] + channel['detector_loss'] # total loss
        eta = 10**(-eta_db/10) # loss
        
        # Entangled source
        eta_db_ES  = channel['channel_loss'] + 2*channel['detector_loss'] # total loss
        eta_ES = 10**(-eta_db_ES/10)
        
        
        # probability of a background being identified as signal - before sifting
        bg_prob = channel['time_filter']*channel['jitter']*channel['bg_rate']
        bg_prob_ES = channel['time_filter']*np.sqrt(2)*channel['jitter']*channel['bg_rate']*10**(-channel['detector_loss']/10) 
                
        
        # time filter to select a coincidence - units of time jitter sigma
        # time filter adds extra loss
        if channel['time_filter'] == 2:
            eta = eta*0.6826
            eta_ES = eta_ES*0.6826
        if channel['time_filter'] == 4:
            eta = eta*0.954        
            eta_ES = eta_ES*0.954


        self.eta = eta
        self.bg_prob = bg_prob
        
        self.eta_ES = eta_ES
        self.bg_prob_ES = bg_prob_ES
        
        if 'channel_fluctuation' in channel:
            self.channel_loss_fluctuation = channel['channel_fluctuation']
        else:
            self.channel_loss_fluctuation = 0 
        
        
    def set_initial_conditions(self, x0):
        self.x = x0

    def set_confidence(self, conf):
        self.conf = conf
        
    
    
    def restriction(self, x):
        s = self.signal(x)
        n1, _ = self.theoretical_vals(x, self.conf)
        
        # print('---> rest', s, n1)
        
        return 2*n1 - s
    
        
    def signal(self, x):
        m, v, pm = x
        pv = 1 - pm
        # s = self.ra*self.time*self.eta*(pm*(1-np.exp(-m)) + pv*(1-np.exp(-v)))
        s = self.ra*self.time*(pm*(self.bg_prob + 1 - np.exp(-self.eta*m)) + pv*(self.bg_prob + 1 - np.exp(-self.eta*v)))
        
        # ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))
        return s
       
        
    def signal_ES(self):
        return self.ra*self.time*(self.eta_ES + self.bg_prob_ES)
    
    
    def QBER_ES(self):
        return (self.eta_ES*self.e_d + self.bg_prob_ES/2) / (self.bg_prob + self.eta_ES)
    
    
    def prob_fail_ES(self):
        N = self.time/(self.time_filter*self.jitter) # Number of time bins 
        
        Ns = self.signal_ES()
        ps = self.ra*self.time*10**(-self.det_loss/10)/N # probability of Alice data point in bin 
        pb = self.bg_rate*self.time/N # probability of background photons per bin
        
        p = ps*pb
        mu_b = N*p
        
        P_fail = N/2*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
        return P_fail
        

    def prob_fail(self, x):
        N = self.time/(self.time_filter*self.jitter) # Number of time bins 
        
        Ns = self.signal(x)
        ps = self.ra*self.time/N # probability of Alice data point in bin 
        pb = self.bg_rate*self.time/N # probability of background photons per bin
        
        p = ps*pb
        mu_b = N*p
        
        P_fail = N/2*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
        return P_fail



    
    

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
    
    
    
    def theoretical_vals_decoy(self, x, confidence):
        m, v, pm = x
        pv = 1 - pm
        
        e = (1 - confidence)/3
        
        
        # Number of sent pulses
        nm = self.ra*self.time*pm
        nv = self.ra*self.time*pv


        # Number of received pulses
        Nm = nm*(self.bg_prob + 1 - np.exp(-self.eta*m))
        Nv = nv*(self.bg_prob + 1 - np.exp(-self.eta*v))
        
        Qm = Nm
        Qv = Nv
        
        
        Qm_u = self.upper_bound(Qm, nm, e, e)/nm
        Qv_u = self.upper_bound(Qv, nv, e, e)/nv  
        Qv_l = self.lower_bound(Qv, nv, e, e)/nv 
        
        
        
        qber_m = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*m)))/(1-np.exp(-self.eta*m) + self.bg_prob)
        qber_v = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*v)))/(1-np.exp(-self.eta*v) + self.bg_prob)
        
        # Number of coincidences divided by two after sifting
        qber_m_u = self.upper_bound(Nm*qber_m/2, Nm/2, e, e)/(Nm/2)
        qber_v_u = self.upper_bound(Nv*qber_v/2, Nv/2, e, e)/(Nv/2)
        
        
        Y1_l = m/(m*v-v**2) * (Qv_l*np.exp(v) - Qm_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * self.bg_prob)
        ns = self.ra*self.time*Y1_l*(pm*np.exp(-m)*m + pv*np.exp(-v)*v)
        
        e1_u_1 = (Qv_u*qber_v_u*np.exp(v) - self.bg_prob/2)/(Y1_l*v)
        e1_u_2 = (Qm_u*qber_m_u*np.exp(m) - self.bg_prob/2)/(Y1_l*m)
        
        # e1_u = min(e1_u_1, e1_u_2)
        # Decoy pulses give better QBER bounds
        e1_u = e1_u_1

        # print('QBERS:',e1_u_1, e1_u_2)
        # print('Nm: ', Nm)

        return ns, e1_u
    
    
    def theoretical_vals_ES(self, confidence):
        e = (1 - confidence)/3
        
        qber = self.QBER_ES()
        n = self.signal_ES()/2
        
        qber_u = self.upper_bound(qber, n, e, e)
        
        return qber_u