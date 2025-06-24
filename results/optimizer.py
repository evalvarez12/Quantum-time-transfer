# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import costfunc

class Optimizer:
    
    def __init__(self, channel):
        self.set_channel(channel)
        
        
    def set_channel(self, channel):
        # Divide by 2 due to basis sifting
        self.ra = channel['ra']
        
        self.time = channel['time']
        self.e_d = channel['e_d']
        
        self.det_loss = channel['detector_loss']
        # Channel must be a function that returns samples of the channel loss
        self.channel_loss = channel['channel_loss']
        
        self.time_filter = channel['time_filter']
        self.jitter = channel['jitter']
        
        
        # eta_db = channel['channel_loss'] + channel['detector_loss'] # total loss
        # eta = 10**(-eta_db/10) # loss
        
        # Entangled source
        # eta_db_ES  = channel['channel_loss'] + 2*channel['detector_loss'] # total loss
        # eta_ES = 10**(-eta_db_ES/10)
        
        
        # probability of a background being identified as signal - before sifting
        Y0 = channel['time_filter']*channel['jitter']*channel['bg_rate']
        Y0_ES = channel['time_filter']*np.sqrt(2)*channel['jitter']*channel['bg_rate']*10**(-channel['detector_loss']/10) 
                
        
        # time filter to select a coincidence - units of time jitter sigma
        # time filter adds extra loss
        # if channel['time_filter'] == 2:
        #     eta = eta*0.6826
        #     eta_ES = eta_ES*0.6826
        # if channel['time_filter'] == 4:
        #     eta = eta*0.954        
        #     eta_ES = eta_ES*0.954


        # self.eta = eta
        self.Y0 = Y0
        
        # self.eta_ES = eta_ES
        self.Y0_ES = Y0_ES
        
        self.confidence =channel['confidence']
        
        
    def channel_mean_decoy(self):
        N = 100000
        etas = self.channel_loss(N) + self.det_loss
        etas = 10**(-etas/10)
        
        eta_mean = np.average(etas)
        
        if self.time_filter == 2:
            eta_mean = eta_mean*0.6826
        if self.time_filter == 4:
            eta_mean = eta_mean*0.954
            
        return eta_mean
    
    def channel_mean_ES(self):
        N = 100000
        etas = self.channel_loss(N) + 2*self.det_loss
        etas = 10**(-etas/10)
        
        eta_mean = np.average(etas)
        
        if self.time_filter == 2:
            eta_mean = eta_mean*0.6826
        if self.time_filter == 4:
            eta_mean = eta_mean*0.954
            
            
        print('eta mean = ', -10*np.log10(eta_mean))
        return eta_mean
        
        
        
    def get_detector_clicks_decoy(self, m, pm):
        # Number of sent pulses
        N = int(self.time*self.ra*pm)
        
        Ni = 0
        div = 1000000
        sum_clicks = 0
        
        use_mean = True
        
        
        if use_mean:
            mean = self.channel_mean_decoy()
            
            sum_clicks = np.floor((1 - np.exp(-m*mean))*N)
        
        else:
            # Do channel sampling in block to dont saturate mem
            while N > 0:
                
                if N > div:
                    Ni = div
                    N = N - div
                else:
                    Ni = N
                    N = 0
                    
                # For each pulse sample the channel
                etas = self.channel_loss(Ni) + self.det_loss
            
                # Transform from dB to transmissivity 
                etas = 10**(-etas/10)
                
                if self.time_filter == 2:
                    etas = etas*0.6826
                if self.time_filter == 4:
                    etas = etas*0.954 
                
                
                # Probability of click for each pulse
                p_click = 1 - np.exp(-m*etas)
                
                roll = np.random.rand(Ni)
                
                clicks = roll < p_click
                sum_clicks += sum(clicks)
            
        return sum_clicks
        
    def get_detector_clicks_ES(self):
        # Number of sent pulses
        N = int(self.time*self.ra)
        
        Ni = 0
        div = 100000
        
        sum_clicks = 0
        
        use_mean = True
        
        
        if use_mean:
            mean = self.channel_mean_ES()
            
            sum_clicks = np.floor(mean*N)
        
        else:
            # Do channel sampling in block to dont saturate mem
            while N > 0:
                
                if N > div:
                    Ni = div
                    N = N - div
                else:
                    Ni = N
                    N = 0
    
                # For each pulse sample the channel
                etas = self.channel_loss(Ni) + 2*self.det_loss
                # Transform from dB to transmissivity 
                etas = 10**(-etas/10)
                
                if self.time_filter == 2:
                    etas = etas*0.6826
                if self.time_filter == 4:
                    etas = etas*0.954 
                
                roll = np.random.rand(Ni)
                
                clicks = roll < etas
                sum_clicks += sum(clicks)
        return sum_clicks
    
        
    def signal(self, x):
        m, v, pm = x
        pv = 1 - pm
        # s = self.ra*self.time*self.eta*(pm*(1-np.exp(-m)) + pv*(1-np.exp(-v)))
        s = self.ra*self.time*(pm*(self.Y0 + 1 - np.exp(-self.eta*m)) + pv*(self.Y0 + 1 - np.exp(-self.eta*v)))
        
        # ra*time*(pmr*(opt.Y0 + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.Y0 + 1 - np.exp(-eta*vr))))
        return s
       
        
    def values_ES(self):
        N = self.get_detector_clicks_ES()
        Nb = self.Y0_ES*self.time*self.ra
        
        q = (N*self.e_d + Nb/2)/(N + Nb)
        
        # N = N + Nb
        
        a = sp.stats.norm.ppf(self.confidence)
        
        std_q = np.sqrt(q*(1-q)/(N/2))
        q_u = q + a*std_q
        print(q, a*std_q, N)
        
        
        return N, q_u
        
        
    def prob_fail_ES(self):
        N = self.time/(self.time_filter*self.jitter) # Number of time bins 
        
        Ns = self.signal_ES()
        ps = self.ra*self.time*10**(-self.det_loss/10)/N # probability of Alice data point in bin 
        pb = self.Y0*self.time/N # probability of background photons per bin
        
        p = ps*pb
        mu_b = N*p
        
        P_fail = N/2*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
        return P_fail
        

    def prob_fail_decoy(self, x):
        m, v, pm = x
        pv = 1 - pm
        
        
        N = self.time/(self.time_filter*self.jitter) # Number of time bins 
        
        Ns = self.signal(x)
        ps = self.ra*self.time/N # probability of Alice data point in bin 
        pb = self.Y0*self.time/N # probability of background photons per bin
        
        p = ps*pb
        mu_b = N*p
        
        P_fail = N/2*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
        return P_fail


    def prob_fail_decoy(self, Qm, Qv, x):
        m, v, pm = x
        pv = 1 - pm  
        
        N = self.time/(self.time_filter*self.jitter)
        Nm = N*pm
        Nv = N*pv
        
        Nsm = self.ra*self.time*pm*Qm
        Nsv = self.ra*self.time*pv*Qv
        Ns = Nsm + Nsv
        
        ps = self.ra*self.time/N # probability of Alice data point in bin
        pb = self.Y0*self.time/N # probability of background photons per bin

        p = ps*pb
        mu_bm = Nm*p
        mu_bv = Nv*p
        mu_b = N*p
        
        P_failm = Nm/2*sp.special.erfc((Nsm)/np.sqrt(2*mu_bm))
        P_failv = Nv/2*sp.special.erfc((Nsv)/np.sqrt(2*mu_bv))
        P_fail = N/2*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
        
        return P_failm, P_failv, P_fail

    def lower_bound(self, x, n, e, e1):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        if (2/e1)**(1/ml) > np.exp(3/(4*np.sqrt(2))**2):
            print('Condition not satisfied lower bound : ', x, n, ml)
            # return 0
        
        D = np.sqrt(2*x*np.log(1/(e1**4/16)))
        return x - D
        
    def upper_bound(self, x, n, e, e2):
        ml = x - np.sqrt(n*np.log(1/e)/2)
        
        if (1/e2)**(1/ml) > np.exp(1/3):
            print('Condition not satisfied upper bound : ', x, n, ml)
            # return 0
        
        D = np.sqrt(2*x*np.log(1/(e2**(3/2))))
        return x + D
    
    
    
    def theoretical_vals_decoy(self, x, finite=True):
        m, v, pm = x
        pv = 1 - pm
        
        e = (1 - self.confidence)/6
        
        a = sp.stats.norm.ppf(1-e)
        
        
        # Number of sent pulses
        nm = self.ra*self.time*pm
        nv = self.ra*self.time*pv

        # Number of received pulses
        Nm = self.get_detector_clicks_decoy(m, pm)
        Nv = self.get_detector_clicks_decoy(v, pv)


        # Number of received pulses
        # Nm = nm*(self.Y0 + 1 - np.exp(-self.eta*m))
        # Nv = nv*(self.Y0 + 1 - np.exp(-self.eta*v))
        
        Qm = Nm/nm
        Qv = Nv/nv
        
        
        Qm_u = self.upper_bound(Nm, nm, e, e)/nm
        Qv_u = self.upper_bound(Nv, nv, e, e)/nv  
        Qv_l = self.lower_bound(Nv, nv, e, e)/nv 
        
        Nb = self.time*self.ra*self.Y0
        # Number of received background counts
        
        Em = (Nb*pm/2 + self.e_d*Nm)/(Nb*pm + Nm)
        Ev = (Nb*pv/2 + self.e_d*Nv)/(Nb*pv + Nv)
        # Ev = (self.Y0/2 + self.e_d*(1-np.exp(-Qv*v)))/(1-np.exp(-self.eta*v) + self.Y0)
        
        # # Number of coincidences divided by two after sifting
        # Em_u = self.upper_bound(Nm*Em/2, Nm/2, e, e)/(Nm/2)
        # Ev_u = self.upper_bound(Nv*Ev/2, Nv/2, e, e)/(Nv/2)
        
        
        # Number of coincidences divided by two after sifting
        std_Em = np.sqrt(Em*(1- Em)/(Nm/2))
        std_Ev = np.sqrt(Ev*(1- Ev)/(Nv/2))
        
        
        Em_u = Em + a*std_Em
        Ev_u = Ev + a*std_Ev
        
        Em_l = Em - a*std_Em
        Ev_l = Ev - a*std_Ev
        
        if finite == False:
            Em_u = Em 
            Ev_u = Ev
            
            Em_l = Em
            Ev_l = Ev
        
            Qm_u = Qm
            Qv_u = Qv  
            Qv_l = Qv 
        
        Y1_l = m/(m*v-v**2) * (Qv_l*np.exp(v) - Qm_u*np.exp(m)*(v/m)**2 - (m**2-v**2)/m**2 * self.Y0)
        ns = self.ra*self.time*Y1_l*(pm*np.exp(-m)*m + pv*np.exp(-v)*v)
        
        e1_u_2 = (Qv_u*Ev_u*np.exp(v) - self.Y0/2)/(Y1_l*v)
        e1_u_1 = (Qm_u*Em_u*np.exp(m) - self.Y0/2)/(Y1_l*m)
        
        e1_u_3 = (Qm_u*Em_u*np.exp(m) - Qv_l*Ev_l*np.exp(v))/(Y1_l*(m-v))
        e1_u_4 = ((v/m)**2*Qm_u*Em_u*np.exp(m) - Qv_l*Ev_l*np.exp(v) + self.Y0/2*(m**2-v**2)/m**2)/(Y1_l*(v**2-m*v)/m)
        
        es = [e1_u_1, e1_u_2, e1_u_3, e1_u_4]
        
        # e1_u = min(e1_u_1, e1_u_2)
        # Decoy pulses give better QBER bounds
        # e1_u = e1_u_2

        # print('QBERS:',e1_u_1, e1_u_2)
        # print('Nm: ', Nm)

        sig = Nv + Nm
        
        
        nsm = self.ra*self.time*Y1_l*(pm*np.exp(-m)*m)
        nsv = self.ra*self.time*Y1_l*(pv*np.exp(-v)*v)
        
        
        # Estimate Key rate for reference
        Q1 = Y1_l*m*np.exp(-m)
        R = .5*(-Qm*self.H2(Em) + Q1*(1-self.H2(e1_u_1)))
        # print('Key rate: ', R)
        
        
        return Qm, Qv, Y1_l, self.Y0, es
    
    
    # def theoretical_vals_ES(self):
    #     e = (1 - self.confidence)/2
        
    #     # Number of detected coincidences
    #     N = self.get_detector_clicks_ES()
        
    #     # Number of detected background
    #     Nb = self.time*self.ra*self.Y0_ES
        
    #     qber = (N*self.e_d + Nb/2) / (N + Nb)
        
    #     Eu = self.upper_bound(N*qber/2, N/2, e, e)/(N/2)
        
    #     return N, Eu
    
    
    def optimize(self, x0):
        # optimization happens on the mean of the transmissivity
        eta = self.channel_mean_decoy()
        n = self.time*self.ra
        Y0 = self.Y0
        
        cf = costfunc.CostFunc(n, eta, self.Y0, self.e_d, self.confidence)
        
        consts = ({'type': 'ineq', "fun": lambda x: x[0] - x[1] - 0.01},
                  {'type': 'ineq', "fun": lambda x: x[0] - 0.099999},
                  {'type': 'ineq', "fun": lambda x: x[1] - 0.099999},
                  {'type': 'ineq', "fun": lambda x: 1 - x[0]},
                  {'type': 'ineq', "fun": lambda x: 1 - x[1]},
                  {'type': 'ineq', "fun": lambda x: x[2] - 0.001},
                  {'type': 'ineq', "fun": lambda x: 1 - x[2]},
                  {'type': 'ineq', "fun": cf.const1},
                  {'type': 'ineq', "fun": cf.const2},
                  {'type': 'ineq', "fun": cf.const3})
        
        res = sp.optimize.minimize(cf.costfunc1, x0, constraints=consts)
        
        return res
    
    
    
    def optimize2(self, x0, emin, ecap):
        # optimization happens on the mean of the transmissivity
        eta = self.channel_mean_decoy()
        n = self.time*self.ra
        Y0 = self.Y0
        
        cf = costfunc.CostFunc(n, eta, self.Y0, self.e_d, self.confidence, emin, ecap)
        
        consts = ({'type': 'ineq', "fun": lambda x: x[0] - x[1] - 0.01},
                  {'type': 'ineq', "fun": lambda x: x[0] - 0.099999},
                  {'type': 'ineq', "fun": lambda x: x[1] - 0.099999},
                  {'type': 'ineq', "fun": lambda x: 1 - x[0]},
                  {'type': 'ineq', "fun": lambda x: 1 - x[1]},
                  {'type': 'ineq', "fun": lambda x: x[2] - 0.001},
                  {'type': 'ineq', "fun": lambda x: 1 - x[2]},
                  {'type': 'ineq', "fun": cf.const1},
                  {'type': 'ineq', "fun": cf.const2},
                  {'type': 'ineq', "fun": cf.const3},
                  {'type': 'ineq', "fun": cf.const_erange})
        
        res = sp.optimize.minimize(cf.costfunc2, x0, method='COBYLA', constraints=consts)
        
        return res
    
    
    
    

        
        
    def H2(self, x):
        return -x*np.log2(x)-(1-x)*np.log2(1-x)

    
    
    