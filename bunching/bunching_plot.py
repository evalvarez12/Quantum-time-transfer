# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 16:00:57 2023

@author: vil034
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)

def Poiss(x, mu):
    return (mu**x)/sp.special.factorial(x) * np.exp(-mu)


def circle_advance(x, n):
    if n == 0:
        return x
    a = np.zeros(len(x))
    a[:n] = x[-n:]
    a[n:] = x[:-n]    
    return a

def g2(a,b, tau):
    g = np.average(a*circle_advance(b, tau))/(np.average(a)*np.average(b))
    return g

def g22(a,b, tau):
    if tau == 0:
        return  np.average(a*b)/(np.average(a)*np.average(b))
    
    a = a[:-tau]
    b = b[tau:]
    
    return np.average(a*circle_advance(b, tau))/(np.average(a)*np.average(b))


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def generate_clicks(p_doubles):
    source_rate = 1e6 # - 1 photon every 4e-8

    T = .01 # total time

    time_jitter = 300e-12
    
    N = int(source_rate * T)
    clicksA = np.random.rand(N) * T
    
    clicksBC = np.copy(clicksA) + np.random.normal(0, time_jitter, size=len(clicksA))
    
    clicksA =  clicksA + np.random.normal(0, time_jitter, size=len(clicksA))
    
    roll_doubles = np.random.rand(N) < p_doubles
    doubles = clicksBC[roll_doubles]
    
    clicksBorC = clicksBC[np.logical_not(roll_doubles)]
    
    
    clicksB = clicksBC[:int(len(clicksBorC)/2)]
    clicksC = clicksBC[int(len(clicksBorC)/2):]
    
    
    clicksA = list(np.sort(clicksA))
    clicksB = list(np.sort(clicksB))
    clicksC = list(np.sort(clicksC))
    doubles = list(doubles)
    
    clicksB = clicksB + doubles
    clicksC = clicksC + doubles

    return clicksA, clicksB, clicksC


def g2(p_doubles):
    clicksA, clicksB, clicksC =  generate_clicks(p_doubles)
    
    N = len(clicksA)
    
    time_window = 300e-12
    
    nB = np.zeros(N)
    nC = np.zeros(N)
    

    for i in range(N):
        
        eventA = clicksA[i]
        
        idB = find_nearest(clicksB, eventA)
        idC = find_nearest(clicksC, eventA)
        
        eventB = clicksB[idB]
        eventC = clicksC[idC]
    
        if np.abs(eventA - eventB) < time_window:
            nB[i] = 1
            del clicksB[idB]
            
        if np.abs(eventA - eventC) < time_window:
            nC[i] = 1
            del clicksC[idC]
            
    return g22(nB, nC, 0)
    
    
    
ps = np.arange(0, 1, 0.04)
gs = np.zeros(len(ps))

for i in range(len(ps)):
    print(i)
    gs[i] = g2(ps[i])

plt.plot(ps, gs, 'bo')

plt.xlabel(r'$p_{doubles}$')
plt.ylabel(r'$g^2$')

ref_gs = np.array([0.0126, 0.023])
ref_ps = [0.00338, 0.00611]

plt.plot(ref_ps, ref_gs, 'ro', label='Reference values')
plt.legend()
