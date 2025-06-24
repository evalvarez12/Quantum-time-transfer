# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import optimizer

class SatellitePass:
    def __init__(self, signal_rate, background_rate, sample_time, sec_integration_time, time_jitter,detector_loss, ed, mode):
        self.signal_rate = signal_rate
        self.background_rate = background_rate
        self.sample_time = sample_time
        self.sec_integration_time = sec_integration_time
        
        self.atmosphere_coh = .1
        
        self.detector_loss = 10**(-detector_loss/10)
        self.time_jitter = time_jitter
        
        self.mode = mode
        
        self.ed = ed
        
        
        
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
    
    
    def get_opt_params(self, ch_loss_avg, ch_loss_std, time, decoy_params=0):
        
        confidence = 0.99
        
        channel_uplink = {
            'ra': self.signal_rate,
            'channel_loss': loss_avg,
            'detector_loss': self.detector_loss,
            'bg_rate': self.background_rate,
            'jitter': self.time_jitter,
            'time_filter': 2,
            'time': time, 
            'e_d': self.e_d,
            'confidence': self.confidence}
        
        channel_downlink = {
            'ra': self.signal_rate,
            'channel_loss': loss_avg,
            'detector_loss': self.detector_loss,
            'bg_rate': self.background_rate,
            'jitter': self.time_jitter,
            'time_filter': 2,
            'time': time, 
            'e_d': self.e_d,
            'confidence': self.confidence}
        
        self.downlink = optimizer.Optimizer(channel_downlink)
        self.uplink = optimizer.Optimizer(channel_uplink)
        
        
        
        opt = optimizer() 
        
        
        
        if mode == 'ES'
        
        
        
        
    # def pnt_ES(self, Nsig, Nbg):
    #     prec = 2*self.time_jitter/np.sqrt(Nsig)
        
    #     qber = (Ns*self.ed + Nb/2)/(Ns + Nb)
        
    #     return prec, qber
    
    
    def pnt_decoy
    
    
    
    def get_clicks_decoy(self, loss, time, x):
        m, v, pm = x
        
        
        # Number of sent pulses
        Nm = int(time*self.signal_rate*pm)
        Nv = int(time*self.signal_rate*(1-pm))
        
        p_click_m = 1 - np.exp(-m*eta)
        roll = np.random.rand(Nm)
        clicks_m = sum(roll < p_click_m)
        
        p_click_v = 1 - np.exp(-v*eta)
        roll = np.random.rand(Nv)
        clicks_v = sum(roll < p_click_v)
            
        return clicks_m, clicks_v
    