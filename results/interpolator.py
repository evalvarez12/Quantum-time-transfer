# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


class Interpolator:
    
    def __init__(self, data, Nbins, plot=False):
        N = 100000
        
        y, x = np.histogram(data, bins=Nbins, density=1)
        

            
        
        Ne = int(N/len(y))
        
        y2 = np.outer(y, np.ones(Ne)).flatten()
        x2 = np.linspace(np.min(x), np.max(x), len(y2))
        
        y2 = y2/sum(y2)
        
        if plot:
            plt.figure()    
            plt.plot(x[1:], y)
            plt.figure()
            plt.plot(x2, y2, 'r.')
        
        self.samples = lambda Nx : np.random.choice(x2, size=Nx, p=y2)
        
        
    def get_samples_dB(self, N):
        return self.samples(N)
    
    
    def get_samples(self, N):
        return 10**(-self.get_samples_dB(N)/10)
        
    
    def get_sampler(self):
        return self.samples
    
        
        
        