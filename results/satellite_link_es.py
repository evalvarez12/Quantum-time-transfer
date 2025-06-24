# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp

class SatelliteLinkES:

    def __init__(self, ra):
        self.ra = ra
        self.bg_rate = 400 
        self.jitter = 364e-12
        self.e_d = 0.05
        self.time_filter = 2
        
        self.confidence = 0.99
        
        self.det_loss = 5
        
        self.time = 1
        
        self.atmosphere_coh = 0.1
        
        Y0 = self.time_filter*np.sqrt(2)*self.jitter*self.bg_rate*10**(-self.det_loss/10)
        self.Y0 = Y0

    

    def run_pnt(self, ch_loss_avg, ch_loss_std):
      
        Ni = int(self.time/self.atmosphere_coh)
        
        Nsig = [0]*Ni
        Nbg = 0
        
        Nb_avg = self.time*self.ra*self.Y0/2
        
        lavg = np.zeros(Ni)
        for i in range(Ni):
            loss_roll = np.random.normal(ch_loss_avg, ch_loss_std)
            loss = loss_roll + 2*self.det_loss
        
            print(str(loss) + ' ', end='')
            lavg[i] = loss
                
            Nsig[i] = self.get_clicks(loss)
            Nbg += np.random.poisson(Nb_avg)
        
            sig_total = sum(Nsig)
        
        # print('avg = ', np.mean(lavg))
        
        return self.values_ES(sig_total,Nbg)

   
    def values_ES(self, Nsig, Nbg):
        q = (Nsig*self.e_d + Nbg/2)/(Nsig + Nbg)
        
        N = Nsig + Nbg 
        
        a = sp.stats.norm.ppf(self.confidence)
        
        std_q = np.sqrt(q*(1-q)/(N/2))
        q_u = q + a*std_q
        
        return N, q_u  
        
 
    
    def get_clicks(self, loss):
        
        eta = 10**(-loss/10)
        
        eta = eta*0.6826
        
        
        
        # Number of sent pulses
        N = int(self.atmosphere_coh*self.ra)
        
        roll = np.random.rand(N)
        clicks = np.sum(roll < eta)
        
        # print('f eta = ', -10*np.log10(eta), clicks, int(N*eta))
        
        return clicks
    

