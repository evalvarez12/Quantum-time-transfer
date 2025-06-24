# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-


import numpy as np
import scipy as sp

def Sph_GT(f, k, Cn2, V, L0, L):
    fun = lambda x: integrand(x, f, V, L0)
    I = sp.integrate.quad(fun, 0, np.inf)
    return 5.211*k**2*Cn2/V * I[0]
    
def integrand(k1, f, V, L0):
    a = (2*np.pi*f/V)**2
    return (a + k1**2 + (1/L0)*(a+k1**2)**(1/2))**(-11/6)

def Sph_K(f, k, Cn2, V, L0, L):
    return 0.0326*k**2*Cn2*f**(-8/3)*V**(5/3)

def Sjitter_GT(f, wl, Cn2, V, L0, L):
    return Sph_GT(f,2*np.pi/wl, Cn2, V, L0, L)*(wl/(2*np.pi*sp.constants.c))**2

def Sjitter_K(f, wl, Cn2, V, L0, L):
    return Sph_K(f,2*np.pi/wl, Cn2, V, L0, L)*(wl/(2*np.pi*sp.constants.c))**2

def Sjitter_K2(f, Cn2, V, L):
    return 0.0326*sp.constants.c**(-2)*L*Cn2*f**(-8/3)*V**(5/4)


def HV_model(h):
    v = 21
    A = 3.2e-13
    Cn2 = 0.00594*(v/27)**2*(1e-5*h)**10*np.exp(-h/1000) \
          +  2.7*1e-16*np.exp(-h/1500) + A*np.exp(-h/100)
  
    return Cn2

def Cn2_layer(hl, hu):
    val = sp.integrate.quad(HV_model, hl, hu)[0]
    return val
    

def CV_profile(h):
    L0 = 4/(1+((h-8500)/2500)**2)
    return L0

def outer_scale(hl, hu):
    val = sp.integrate.quad(CV_profile, hl, hu)[0] / (hu - hl)
    return val





import matplotlib.pyplot as plt

# Phase screen positions
# z = 0 - H = 600 km
H = [600000, 20000, 18938.38, 16882.21, 15337.21, 14057.94, 12937.75, 11918.21, 10961.72, 10040.25, 9129.51, 8205.05, 7239.94,
     6210.12, 5134.73, 4132.22, 3288.03, 2582.98, 1980.88, 1455.63, 989.72, 574.45, 262.94, 100.83] 

V = 0.6
wl = 810e-9
N = 500
fi = 0.001
ff =30

fs = np.linspace(fi, ff, N)

df = 1/ff

Nit = 10000
dts_GT = np.zeros(Nit)
dts_K = np.zeros(Nit)

fGT = 0
fK = 0

zs = [0, 10, 20, 30, 40, 50, 60]
jitter = []

for z in zs:

    z = np.deg2rad(zs)
    
    for hi in range(len(H)-1):
        hu = H[hi]
        hl = H[hi+1]
        
        Cn2 = Cn2_layer(hl, hu)
        L0 = outer_scale(hl, hu)
    
    
        L = (hu - hl) * 1/np.cos(z)
    
        psd_K = Sjitter_K(fs, wl, Cn2, V, L0, L)
        # psd_K = Sjitter_K2(fs, Cn2, V, L)
    
    
        psd_GT = []
        for j in fs: 
            psd_GT += [Sjitter_GT(j, wl, Cn2, V, L0, L)]
    
        for i in range(Nit):
    
            # Generate random phase angles
            phase = np.random.rand(len(psd_K)) * 2 * np.pi 
        
            # Construct complex Fourier coefficients
            fft_coeffs_GT = np.sqrt(psd_GT) * np.exp(1j * phase) 
            
            fft_coeffs_K = np.sqrt(psd_K) * np.exp(1j * phase) 
        
            # Perform the inverse Fourier transform
            signal_GT = np.fft.irfft(fft_coeffs_GT)
            
            signal_K = np.fft.irfft(fft_coeffs_K)
            
            # signal_GT = signal_GT + (fGT - signal_GT[0])
            # signal_K = signal_K + (fK - signal_K[0])
            
            dts_GT[i] += np.sum(signal_GT)*1e12
            dts_K[i] += np.sum(signal_K)*1e12
            
            # fGT = signal_GT[-1]
            # fK = signal_K[-1]
    jitter += [np.std(dts_GT)]
    
print('-----------------------')
# print(np.mean(dts))
print('Time jitter Greenwood Tarazano:', np.std(dts_GT))
print('Time jitter Kolmogorov:', np.std(dts_K))

plt.figure()
plt.plot(zs, jitter, linewidth=3)

plt.xlabel('Zenith angle (deg)',fontsize=13)
plt.ylabel('Jitter (ps)', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)