# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import interpolator
import optimizer

class Satellite:
    
    def __init__(self, ra=20e6, time=1):
        self.passes = 5
        self.ra = ra
        # self.bg_rate = 400
        self.bg_rate = 2400
        self.jitter = 364e-12
        self.e_d = 0.05
        self.time_filter = 2
        
        self.confidence = 0.99
        
        self.detector_loss = 5
        
        self.time = time
        
        self.uplink = None
        self.downlink = None
        
        
        self.p1 = {
            '35-45' : 150
            }
        
        self.p2 = {
            '35-45' : 37,
            '45-55' : 69,
            '55-65' : 96,
            '65-75' : 120,
            '75-85' : 142,
            '85-90-85' : 18,
            '85-75' : 167,
            '75-65' : 144,
            '65-55' : 121,
            '55-45' : 94,
            '45-35' : 61,
            '35-30' : 19
            }
        
        self.p3 = {
            '35-45' : 46,
            '45-55' : 83,
            '55-65' : 122,
            '65-70-65' : 329,
            '65-55' : 148,
            '55-45' : 109,
            '45-35' : 72,
            '35-30' : 27
            }
        
        self.p4 = {
            '30-35' : 22,
            '35-45' : 64,
            '45-55' : 97,
            '55-65' : 125,
            '65-75' : 150,
            '75-85-75' : 359,
            '75-65' : 151,
            '65-55' : 125,
            '55-45' : 98,
            '45-35' : 65,
            '35-30' : 23
            }
        
        self.p5 = {
            '30-35-30' : 59
            }


    def extinction(self, z):
        return np.exp(-0.7 * np.sec(np.deg2rad(z)))
    
    
    
    def get_data_uplink(self, z):
        file = f'data/PAPER_UPLINK_z={z}_wl=0.81_w0=0.35_ap=0.2_H=500_perr=1.2e-06'
        data = sp.io.loadmat(file)
        data = data['res'].transpose()[0]
        
        return data
        
    
    def get_data_downlink(self, z):
        file = f'data/PAPER_DOWNLINK_z={z}_wl=0.81_w0=0.2_ap=0.35_H=500_perr=1.2e-06'
        data = sp.io.loadmat(file)
        data = data['res'].transpose()[0]
        
        return data
        
        
    def setup_channels(self, z, verbose=False, ext=[0, 0]):

        extAtmos = np.exp(-.7/np.cos(np.deg2rad(z)))
        extU = 10**(-ext[0]/10)
        extD = 10**(-ext[1]/10)
        data_uplink = self.get_data_uplink(z) * extAtmos * extU 
        data_downlink = self.get_data_downlink(z) * extAtmos * extD
    
        if verbose:
            print('Loss Uplink avg: ', -10*np.log10(np.mean(data_uplink)) + 10)
            # print('Loss Downlink avg: ', -10*np.log10(np.mean(data_downlink)))
            
            print('expected signal  = ', 20e6*(np.mean(data_uplink)*.1))
            
        data_uplink = -10*np.log10(data_uplink)
        data_downlink = -10*np.log10(data_downlink)


        Nbins = 60
        interpol_uplink = interpolator.Interpolator(data_uplink, Nbins)
        interpol_downlink = interpolator.Interpolator(data_downlink, Nbins)

    
        channel_uplink = {
            'ra': self.ra,
            'channel_loss': interpol_uplink.get_sampler(),
            'detector_loss': self.detector_loss,
            'bg_rate': self.bg_rate,
            'jitter': self.jitter,
            'time_filter': self.time_filter,
            'time': self.time, 
            'e_d': self.e_d,
            'confidence': self.confidence}
        
        channel_downlink = {
            'ra': self.ra,
            'channel_loss': interpol_downlink.get_sampler(),
            'detector_loss': self.detector_loss,
            'bg_rate': self.bg_rate,
            'jitter': self.jitter,
            'time_filter': self.time_filter,
            'time': self.time, 
            'e_d': self.e_d,
            'confidence': self.confidence}
        
        self.downlink = optimizer.Optimizer(channel_downlink)
        self.uplink = optimizer.Optimizer(channel_uplink)
    
    
        # self.downlink_channel = channel_downlink
        # self.uplink_channel = channel_uplink
    
    
    def get_vals_ES(self):
        vals_downlink = self.downlink.values_ES()
        vals_uplink = self.uplink.values_ES()
        return vals_uplink, vals_downlink
    
    
    def get_vals_decoy(self, x0_up, x0_down):
        opt_up = self.uplink.optimize(x0_up)
        opt_down = self.downlink.optimize(x0_down)
        
        # print("Opt UP")
        # print(opt_up)
        # x = opt_up['x']
        
        # vals_up = self.uplink.theoretical_vals_decoy(x)
        
        # print("Opt DOWN")
        # print(opt_down)
        # x = opt_down['x']
        
        # vals_down = self.downlink.theoretical_vals_decoy(x)
        
        return opt_up, opt_down