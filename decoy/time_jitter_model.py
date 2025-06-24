# -*- coding: utf-8 -*-


import numpy as np
import scipy as sp

def Sph_GT(f, k, Cn2, V, L0, L):
    fun = lambda x: integrand(x, f, V, L0)
    I = sp.integrate.quad(fun, 0, np.inf)
    return 5.211*k**2*Cn2*L/V * I[0]
    
def integrand(k1, f, V, L0):
    a = (2*np.pi*f/V)**2
    return (a + k1**2 + (1/L0)*(a+k1**2)**(1/2))**(-11/6)

def Sph_K(f, k, Cn2, V, L0, L):
    return 0.0326*k**2*Cn2*L*f**(-8/3)*V**(5/3)

def Sjitter_GT(f, wl, Cn2, V, L0, L):
    return Sph_GT(f,2*np.pi/wl, Cn2, V, L0, L)*(wl/(2*np.pi*sp.constants.c))**2

def Sjitter_K(f, wl, Cn2, V, L0, L):
    return Sph_K(f,2*np.pi/wl, Cn2, V, L0, L)*(wl/(2*np.pi*sp.constants.c))**2

def Sjitter_K2(f, Cn2, V, L):
    return 0.0326*sp.constants.c**(-2)*L*Cn2*f**(-8/3)*V**(5/4)


# L = 137
L = 12000
V = .63
Cn2 = 3.2e-13
f = .01
L0 = 1.6
wl = 1560e-9 

print(Sjitter_GT(f, wl, Cn2, V, L0, L))
print(Sjitter_K(f, wl, Cn2, V, L0, L))
print(Sjitter_K2(f, Cn2, V, L))

# print(np.sqrt(Sjitter_K(f, Cn2, V, L))*1e12)




import matplotlib.pyplot as plt

N = 500
fi = 0.001
ff = 50

fs = np.linspace(fi, ff, N)

df = 1/ff

psd_K = Sjitter_K(fs, wl, Cn2, V, L0, L)
# psd_K = Sjitter_K2(fs, Cn2, V, L)


psd_GT = []
for i in fs: 
    psd_GT += [Sjitter_GT(i, wl, Cn2, V, L0, L)]


# Generate random phase angles
phase = np.random.rand(len(psd_K)) * 2 * np.pi

# Construct complex Fourier coefficients
fft_coeffs = np.sqrt(psd_K) * np.exp(1j * phase) 

# Perform the inverse Fourier transform
signal = np.fft.irfft(fft_coeffs)


# phase = (np.random.rand(N) + 1j*np.random.rand(N))*np.sqrt(psd)

# signal = np.fft.irfft(phase)


# Print or plot the reconstructed signal
# print("Reconstructed Signal:", signal)

print('---------------------------------')

print(np.mean(signal)*1e12)
print(np.std(signal)*1e12)

print(np.sum(signal)*1e12)

plt.close('all')
plt.figure()
plt.plot(signal)

# plt.figure()
# plt.plot(fs, psd_GT)
# plt.plot(fs, psd_K)

# plt.yscale('log')
# plt.xscale('log')

plt.show()



Nit = 10000
dts_GT = np.zeros(Nit)
dts_K = np.zeros(Nit)

for i in range(Nit):
    # Generate random phase angles
    phase = np.random.rand(len(psd_GT)) * 2 * np.pi 

    # Construct complex Fourier coefficients
    fft_coeffs_GT = np.sqrt(psd_GT) * np.exp(1j * phase) 

    fft_coeffs_K = np.sqrt(psd_K) * np.exp(1j * phase) 

    # Perform the inverse Fourier transform
    signal_GT = np.fft.irfft(fft_coeffs_GT)
    
    signal_K = np.fft.irfft(fft_coeffs_K)

    
    dts_GT[i] = np.sum(signal_GT)*1e12
    dts_K[i] = np.sum(signal_K)*1e12
    
print('-----------------------')
# print(np.mean(dts))
print('Time jitter Greenwood Tarazano:', np.std(dts_GT))
# print('Time jitter Kolmogorov:', 4*np.std(dts_K))
print('Time jitter Kolmogorov:', np.std(dts_K))

plt.figure()
plt.hist(dts_K, 100)