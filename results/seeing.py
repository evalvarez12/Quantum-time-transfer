# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:35:10 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# plt.rcParams['text.usetex'] = True

ts = sp.io.loadmat('data/Z1_T.mat')
cx = sp.io.loadmat('data/Z1_ZX.mat')
cy = sp.io.loadmat('data/Z1_ZY.mat')

ts = ts['res']
cx = cx['csx']
cy = cy['csy']



wl = 810e-9
d = 0.15


mx = np.arctan(wl*cx/(d*np.pi))
my = np.arctan(wl*cy/(d*np.pi))


plt.hist(mx, 100, density=True)

plt.xlabel(r'$\theta_x$', fontsize=13)
plt.ylabel(r'PDF', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)

# plt.hist(my, 200, alpha=0.5)