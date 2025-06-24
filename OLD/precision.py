# -*- coding: utf-8 -*-
"""
Created on Mon May 29 10:30:18 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import gauss_fit as gf


alice_times = np.load('alice_points.npy')
bob_times = np.load('bb_points.npy')

diff = alice_times - bob_times

pfit = gf.fit_w_hist(diff, 200)


print('Offset: ', pfit[1])
print('Precision: ', pfit[2])
