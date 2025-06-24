# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 11:56:45 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def metaPoiss(x, mu, alpha):
    return ((alpha*mu)**(alpha*x))/sp.special.factorial(alpha*x) * np.exp(-mu*alpha)




x = np.arange(0, 20, 0.001)
alpha = .1
mu = 1


p1 = metaPoiss(x, mu, 1)
p2 = metaPoiss(x, mu, alpha)

p2 = np.abs(p2)

p1 = p1/np.sum(p1)
p2 = p2/np.sum(p2)
    
n = 10000
s1 = np.random.choice(x, n, p=p1)

s2 = np.random.choice(x, n, p=p2)

s3 = np.random.poisson(mu, n)

m1 = np.average(s1)
m2 = np.average(s2)

var1 = np.average(s1**2) - np.average(s1)**2
var2 = np.average(s2**2) - np.average(s2)**2


print(m1)
print(m2)

print(var1)
print(var2)

n1 = np.zeros(20)
n2 = np.zeros(20)
for i in range(20):
    n1[i] = np.sum(s1 == i)
    n2[i] = np.sum(s2 == i)


print(np.average(n1**2)/np.average(n1)**2)
print(np.average(n2**2)/np.average(n2)**2)

plt.close('all')

# plt.figure()
# plt.plot(x, metaPoiss(x, mu, 1))
# plt.plot(x, metaPoiss(x, mu, alpha))


# plt.hist(s1, bins=50, density=True, color='b', alpha=0.5)
# plt.hist(s2, bins= 50, density=True, color='r', alpha=0.5)
# plt.hist(s3, bins= 50, density=True, color='g', alpha=0.5)

print()

##################################
n = 20
dt = 1
# mu = .1

T = dt*n
Ntot = 10

clicks = np.random.rand(int(Ntot))*T
bins = np.arange(0,T, dt)

counts = np.zeros(50)

nt = np.zeros(n)

for i in range(len(bins)):
    t = bins[i]
    ind =np.sum(np.logical_and(clicks > t, clicks < t+dt))
    counts[ind] += 1

    nt[i] = ind

counts = counts/np.sum(counts)

x = np.arange(0,50,1)

plt.figure()
plt.plot(x, counts, 'o')

f = lambda xa,ma: metaPoiss(xa, ma, 1)

fit = sp.optimize.curve_fit(f, x, counts)

mu = fit[0][0]
plt.plot(x, metaPoiss(x, mu, 1))
print('mu: ', mu)


print(np.average(nt**2)/np.average(nt)**2)

plt.figure()
plt.plot(nt, 'k.')
plt.plot(clicks/dt, 0.5*np.ones(len(clicks)), 'r.')