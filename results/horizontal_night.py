# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import optimizer 
import interpolator


plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"


################# Loss plot - day & night

file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
file_night = 'data/HORIZONTAL_DISTANCE_NIGHT2_trx=15_L=%i_wl=0.81'

# dists = [7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000]
dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]


m_d = []
std_d = []
m_n = []
std_n = []


scattering_loss = 0.45 # db/km for 810


#################### QBER plot - night - decoy & ES
jitter = 364e-12
time = 2
sig_rate = 20e6
det_loss = 5
time_filter = 2
qber_dev = 0.05

bg_day = 483924
bg_night = bg_day/50

channel = {
    'ra': sig_rate,
    'channel_loss': 0,
    'detector_loss': det_loss,
    'bg_rate': 0,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': qber_dev}


# d = 7000

m = 0.5 # signal mean photon number
v = 0.1 # decoy mean photon number  - required v < m
pm = 0.6 # probability of sigal produced

x = [m, v, pm]

confidence = 0.99

sigs = []
ns = []
qs = []

sigs_es = []
qs_es = []

for d in dists:

    data = sp.io.loadmat(file_night %d)
    data = data['res'].transpose()[0]
    data = -10*np.log10(data)+scattering_loss*d*1e-3

    interpol = interpolator.Interpolator(data, 60)
    channel['channel_loss'] = interpol.get_sampler()

    channel['bg_rate'] = bg_night
    opt = optimizer.Optimizer(channel)

    res= opt.optimize(x, confidence)
    
    print(res)
    
    xo = res['x']

    s, n, q = opt.theoretical_vals_decoy(xo, confidence)
    
    sigs += [s]
    ns += [n]
    qs += [q]
    
    s, q = opt.theoretical_vals_ES(confidence)
    
    sigs_es += [s]
    qs_es += [q]


sigs = np.array(sigs)

plt.figure()
plt.plot(dists, sigs)
plt.plot(dists, ns)
plt.plot(dists, sigs/2, 'r--')

plt.plot(dists, sigs_es, 'o')


plt.figure()
plt.plot(dists, qs)
plt.plot(dists, qs_es)
plt.ylim([0, 0.25])






