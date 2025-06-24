# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 10:35:25 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt 
import scipy as sp


def Poiss(x, mu):
    return (mu**x)/sp.special.factorial(x) * np.exp(-mu)


def circle_advance(x, n):
    if n == 0:
        return x
    a = np.zeros(len(x))
    a[:n] = x[-n:]
    a[n:] = x[:-n]    
    return a

def metaPoiss(x, mu, alpha):
    return ((alpha*mu)**(alpha*x))/sp.special.factorial(alpha*x) * np.exp(-mu*alpha)

n = 100000 # number of bins
mu = .4

tau =  1
x = np.arange(0, n, 1)
p = np.abs(Poiss(x, mu))

p = np.zeros(n)
p[1] = .1
p[0] = 1 - p[1]
# 
nt = np.random.choice(x, size=n, p=p)

# g2 = np.average(nt**2)/np.average(nt)**2
g2 = np.average(nt*circle_advance(nt, tau))/np.average(nt)**2


print('g2 = ', g2)
print('avg nt: ', np.average(nt))