# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import interpolator

plt.close('all')

d = 7000
file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
data = sp.io.loadmat(file_day %d)

scattering_loss = 0.45 # db/km for 810
data = data['res'].transpose()[0]
data = -10*np.log10(data)+scattering_loss*d*1e-3


Nbins = 60
interpol = interpolator.Interpolator(data, Nbins, plot=True)

N = 5000
samps = interpol.get_samples_dB(N)

plt.figure()
plt.hist(samps, 80, alpha=0.5, density=True)
plt.hist(data, 80, alpha=0.5, density=True)
