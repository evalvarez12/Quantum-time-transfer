# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:51:35 2023

@author: vil034
"""

import numpy as np


fov = 4*np.pi*np.sin(25e-6/2)**2
area = np.pi * (.15/2)**2
filt = 1 
lamb = 810

I0 = 0.004

eta = 0.6**2

# I0 = 2e-7
I = I0 * fov * area * filt

hc = 1.98644586e-25

N = I/hc * lamb*1e-9 
N = N*eta

print('Background counts: ', N*1e-3)
print('---------------------------')

######## Paper comparisson

N = 500000
I = N*hc/(lamb*1e-9)

I0 = I / (fov * area * filt)
print(I0/0.6)



fov2 = 4*np.pi*(np.sin(2*283e-6)**2)
I = I0 * fov2 * area * filt

N = I/hc * lamb*1e-9 
print(N*1e-6)


###############################################

filt = 1
fov = 4*np.pi*(np.sin(566e-6)**2)
a = np.pi * (.15/2)**2
lamb = 785

dk = 200*4

Eph = hc/lamb
N = 1415 - dk
E = Eph * N

L = E/(fov*a*filt)

filt_o = 3
fov_o = 4*np.pi*(np.sin(20e-3)**2)
a_o = np.pi * .15**2

E_o = L * filt_o * fov_o * a_o
N_o = E_o/hc * lamb
print(N_o + dk)

