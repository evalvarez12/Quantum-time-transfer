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
        
    
    def cost_Y1(self, x):
        n, e1 = self.theoretical_vals(x, self.conf)
        
        prob_ph = 1- np.exp(-x[0])
        sec_lim = (self.ra*self.time*self.eta*x[2]*prob_ph)/2
        
        cost = - (n - sec_lim)
        
        # if e1 > .2:
            # order_n = np.floor(np.log10(n)) + 1
            # cost = e1*10**(order_n) - n 
        
        # cost = e1
        # print(cost)

        return cost
    
    
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


    def solve(self, x):
        m, v, pm = x
    
        f = lambda mv: self.restriction([mv, v, pm,])
        mr = opt.fsolve(f, m)[0]
        return mr
        
    def cost_S(self, x):
        v, pm = x
        
        m = self.solve([self.x[0], v, pm])

        return -self.signal([m, v, pm])
    
    

    def lower_bound(x, n, e, e1):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        if (2/e1)**ml > np.exp(3/(4*np.sqrt(2))**2):
            print('Condition not satisfied lower bound : ', x)
            return 0
        
        D = np.sqrt(2*x*np.log(1/(e1**4/16)))
        return x - D
        
    def upper_bound(x, n, e, e2):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        if (1/e2)**ml > np.exp(1/3):
            print('Condition not satisfied upper bound : ', x)
            return 0
        
        D = np.sqrt(2*x*np.log(1/(e2**(3/2))))
        return x + D
    
    
    
    def theoretical_vals_fancy_bounds(self, x, confidence):
        m, v, pm = x
        pv = 1 - pm
        
        e = (1 - confidence)/3
        
        
        # Number of sent pulses
        nm = self.ra*self.time*pm
        nv = self.ra*self.time*pv


        # Number of received pulses
        Nm = nm*(self.bg_prob + 1 - np.exp(-self.eta*m))
        Nv = nv*(self.bg_prob + 1 - np.exp(-self.eta*v))
        
        Qm = Nm/nm
        Qv = Nv/nv
        
        
        Qm_u = self.upper_bound(Qm, nm, e, e)
        Qv_u = self.upper_bound(Qv, nv, e, e)  
        Qv_l = self.lower_bound(Qv, nv, e, e) 
        
        
        qber_m = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*m)))/(1-np.exp(-self.eta*m) + self.bg_prob)
        qber_v = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*v)))/(1-np.exp(-self.eta*v) + self.bg_prob)
        
        # Number of coincidences divided by two after sifting
        qber_m_u = self.upper_bound(qber_m, Nm/2, e, e)
        qber_v_u = self.upper_bound(qber_v, Nv/2, e, e)
        
        
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
    
    
    def theoretical_vals(self, x, confidence):
        m, v, pm = x
        pv = 1 - pm
        
        a = sp.stats.norm.ppf(confidence)
        
        # Number of sent pulses
        nm = self.ra*self.time*pm
        nv = self.ra*self.time*pv


        # Number of received pulses
        Nm = nm*(self.bg_prob + 1 - np.exp(-self.eta*m))
        Nv = nv*(self.bg_prob + 1 - np.exp(-self.eta*v))
        
        Qm = Nm/nm
        Qv = Nv/nv
        
        if self.channel_loss_fluctuation != 0:
            Qm_u =  Qm + a*self.channel_loss_fluctuation/np.sqrt(nm)
            Qv_l =  Qv - a*self.channel_loss_fluctuation/np.sqrt(nv)
            Qv_u =  Qv + a*self.channel_loss_fluctuation/np.sqrt(nv)
            
            std_Qmu = np.sqrt(Qm_u*(1- Qm_u)/nm)
            std_Qvu = np.sqrt(Qv_u*(1- Qv_u)/nv)
            std_Qvl = np.sqrt(Qv_l*(1- Qv_l)/nv)
            
        else: 
            Qm_u = Qm
            Qv_l = Qv
            Qv_u = Qv

            std_Qmu = np.sqrt(Qm*(1- Qm)/nm)
            std_Qvu = np.sqrt(Qv*(1- Qv)/nv)
            std_Qvl = np.sqrt(Qv*(1- Qv)/nv)
        
        
        qber_m = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*m)))/(1-np.exp(-self.eta*m) + self.bg_prob)
        qber_v = (self.bg_prob/2 + self.e_d*(1-np.exp(-self.eta*v)))/(1-np.exp(-self.eta*v) + self.bg_prob)
        
        # Number of coincidences divided by two after sifting
        std_qber_m = np.sqrt(qber_m*(1- qber_m)/(Nm/2))
        std_qber_v = np.sqrt(qber_v*(1- qber_v)/(Nv/2))
        
        
        qber_m_u = qber_m + a*std_qber_m
        qber_v_u = qber_v + a*std_qber_v
        
        Qm_u = Qm_u + a*std_Qmu
        Qv_l = Qv_l - a*std_Qvl 
        Qv_u = Qv_u + a*std_Qvu 
        
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
    


    def optimize(self):
        
        x0 = self.x
        cons = ({'type': 'ineq', "fun": lambda x: x[0] - x[1] - 0.01},
                {'type': 'ineq', "fun": lambda x: x[0] - 0.0001},
                {'type': 'ineq', "fun": lambda x: x[1] - 0.0001},
                {'type': 'ineq', "fun": lambda x: 1 - x[0]},
                {'type': 'ineq', "fun": lambda x: 1 - x[1]},
                {'type': 'ineq', "fun": lambda x: x[2] - 0.001},
                {'type': 'ineq', "fun": lambda x: 1 - x[2]})

        res = opt.minimize(self.cost_Y1, x0, constraints=cons)

        self.x = res['x']
        print(res)


    def optimize2(self):
        
        # self.m0 = self.x
        
        
        cons = ({'type': 'ineq', "fun": lambda x: x[0] - 0.1},
                {'type': 'ineq', "fun": lambda x: x[1] - 0.2},
                {'type': 'ineq', "fun": lambda x: 0.4 - x[0]},
                {'type': 'ineq', "fun": lambda x: 0.4 - x[1]})

        res = opt.minimize(self.cost_S, self.x[1:], constraints=cons)

        print(res)

        if not res['success']:
            print('-----> Optimization failed')

        v, pm = res['x']
        m = self.solve([self.x[0], v, pm])
        # return [m, v, pm]
        self.x = np.array([m, v, pm])
    
        