# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import optimizer
import interpolator


class SatelliteLinkDecoy:
    def __init__(self, channel):
        self.set_channel(channel)
        
        
    def set_channel(self, ra):
        self.ra = ra
        self.bg_rate = 400 
        self.jitter = 364e-12
        self.e_d = 0.05
        self.time_filter = 2
        
        self.confidence = 0.99
        
        self.det_loss = 5
        
        self.time = 1
        
        self.atmosphere_coh = 0.1
        
        Y0 = self.time_filter*self.jitter*self.bg_rate
        self.Y0 = Y0
        
        

    

    def run_pnt(self, ch_loss_avg, ch_loss_std):
        self.set_opt_params(ch_loss_avg)
        
        
        Ni = int(self.time/self.atmosphere_coh)
        
        Nm = 0
        Nv = 0
        Nbg = 0
        
        Nb_avg = self.time*self.ra*self.Y0/2
        
        lavg = np.zeros(Ni)
        for i in range(Ni):
            loss_roll = np.random.normal(ch_loss_avg, ch_loss_std)
            loss = loss_roll + self.det_loss
            
            # print(str(loss) + ' ', end='')
            lavg[i] = loss
        
            a, b = self.get_clicks(loss)
            Nm += a
            Nv += b
            Nbg += np.random.poisson(Nb_avg)
        
        # print('avg rolled = ', np.average(lavg))
        
        return self.theoretical_vals_decoy(Nm, Nv, Nbg)
    
    def set_opt_params(self, ch_loss_avg):
        
        ch_loss = interpolator.Interpolator([ch_loss_avg], 1).get_sampler()
        
        channel = {
            'ra': self.ra,
            'channel_loss': ch_loss,
            'detector_loss': self.det_loss,
            'bg_rate': self.bg_rate,
            'jitter': self.jitter,
            'time_filter': self.time_filter,
            'time': self.time, 
            'e_d': self.e_d,
            'confidence': self.confidence}
        

        
        link = optimizer.Optimizer(channel)
        
        m = 0.5 # signal mean photon number
        v = 0.3 # decoy mean photon number  - required v < m
        pm = 0.2
        x0 = [m, v, pm] 
        
        opt_res = link.optimize(x0)
        
        if not opt_res['success']:
            print("Optimization process fail")
        
        self.x = opt_res['x']
        
        
    def get_opt_params(self):
        return self.x
 
    
    def get_clicks(self, loss):
        # print('click')
        m, v, pm = self.x
        
        eta = 10**(-loss/10)
        
        # Number of sent pulses
        Nm = int(self.atmosphere_coh*self.ra*pm)
        Nv = int(self.atmosphere_coh*self.ra*(1-pm))
        
        p_click_m = 1 - np.exp(-m*eta)
        roll = np.random.rand(Nm)
        clicks_m = np.sum(roll < p_click_m)
        
        p_click_v = 1 - np.exp(-v*eta)
        roll = np.random.rand(Nv)
        clicks_v = np.sum(roll < p_click_v)
            
        return clicks_m, clicks_v
    

    
    
    
    def theoretical_vals_decoy(self, Nm, Nv, Nb):
        m, v, pm = self.x
        pv = 1 - pm
        
        e = (1 - self.confidence)/6
        
        a = sp.stats.norm.ppf(1-e)
        
        
        # Number of sent pulses
        nm = self.ra*self.time*pm
        nv = self.ra*self.time*pv

        # # Number of received pulses
        # Nm = self.get_clicks_decoy(m, pm)
        # Nv = self.get_clicks_decoy(v, pv)


        Qm = Nm/nm
        Qv = Nv/nv
        
        
        Qm_u = self.upper_bound(Nm, nm, e, e)/nm
        Qv_u = self.upper_bound(Nv, nv, e, e)/nv  
        Qv_l = self.lower_bound(Nv, nv, e, e)/nv 
        
        # Number of received background counts - afer sifting
        # Nb_avg = self.time*self.ra*self.Y0/2
        # Nb = np.random.poisson(Nb_avg)
        
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
        
        
        return Qm, Qv, Y1_l, self.Y0, es
    
    
    
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
    
    
    
    
    
    
    # def sample_data(self, loss_avg, loss_std, time):
    #     loss = np.random.normal(loss_avg, loss_std)
        
    #     Nsig = int(0.68*self.signal_rate*time*10**(-loss/10))    
        
        
    #     bg_mean = self.background_rate*time*self.signal_rate*2*self.time_jitter
    #     if self.mode == 'ES':
    #         bg_mean = bg_mean*self.detector_loss
        
    #     Nbg = np.random.poisson(bg_mean)
    #     return Nsig, Nbg
    
    
    # def get_clicks_ES(self, loss_avg, loss_std, time):
    #     Ni = int(time/self.atmosphere_coh)
        
    #     Nsig = 0
    #     Nbg = 0
    #     for i in range(Ni):
    #         a, b = self.sample_data(loss_avg, loss_std, time)
    #         Nsig += a
    #         Nbg += b
            
    #     return Nsig, Nbg
    
        
    # def pnt_ES(self, Nsig, Nbg):
    #     prec = 2*self.time_jitter/np.sqrt(Nsig)
        
    #     qber = (Ns*self.ed + Nb/2)/(Ns + Nb)
        
    #     return prec, qber
    
       