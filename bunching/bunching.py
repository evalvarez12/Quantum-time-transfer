# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:04:17 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def metaPoiss(x, mu, alpha):
    return ((alpha*mu)**(alpha*x))/sp.special.factorial(alpha*x) * np.exp(-mu*alpha)





source_rate = 25e6 # 1 kHz

dt = 500e-12 # integration time - 1 photon every 4e-8

dt = 1e-8
T = .0001 # total time
n = int(T / dt) # number of bins


N = int(source_rate * T)
clicks = np.random.rand(N) * T

nt = np.zeros(n)

for i in range(n):
    t = i * dt
    counts = np.sum(np.logical_and(clicks > t, clicks < t+dt))
    
    nt[i] = counts


g2 = np.average(nt**2)/np.average(nt)**2
# g2 = np.average(nt**2)/np.average(nt)**2

print('g2 = ', g2)
print('avg nt: ', np.average(nt))


plt.hist(nt)

# y, x = np.histogram(nt, bins=5)

# x = x[0:-1]

# plt.figure()
# plt.plot(x, y, 'o')

# f = lambda xa,ma: metaPoiss(xa, ma, 1)

# fit = sp.optimize.curve_fit(f, x, y)

# mu = fit[0][0]
# plt.plot(x, metaPoiss(x, mu, 1))
# print('mu: ', mu)
