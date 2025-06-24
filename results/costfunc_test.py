# -*- coding: utf-8 -*-


import numpy as np
import costfunc 
import scipy as sp
import interpolator
import optimizer
import scipy as sp

jitter = 364e-12
time = 2
sig_rate = 20e6
det_loss = 5
time_filter = 2
qber_dev = 0.05

bg_day = 483924
bg_night = bg_day/50

d = 7000
file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
data = sp.io.loadmat(file_day %d)

scattering_loss = 0.45 # db/km for 810
data = data['res'].transpose()[0]
data = -10*np.log10(data)+scattering_loss*d*1e-3

interpol = interpolator.Interpolator(data, 60)

channel = {
    'ra': sig_rate,
    'channel_loss': 0,
    'detector_loss': det_loss,
    'bg_rate': bg_night,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': qber_dev}

channel['channel_loss'] = interpol.get_sampler()
opt = optimizer.Optimizer(channel)

n = time*sig_rate
bg_prob = opt.bg_prob
eta = opt.channel_mean_decoy()

confidence = 0.99

m = 0.5 # signal mean photon number
v = 0.3 # decoy mean photon number  - required v < m
pm = 0.6 # probability of sigal produced

x = [m, v, pm]

cf = costfunc.CostFunc(n, eta, bg_prob, qber_dev, confidence)

print(cf.vals(x))
print(cf.consts)

x0 = x
consts = ({'type': 'ineq', "fun": lambda x: x[0] - x[1] - 0.01},
          {'type': 'ineq', "fun": lambda x: x[0] - 0.0001},
          {'type': 'ineq', "fun": lambda x: x[1] - 0.0001},
          {'type': 'ineq', "fun": lambda x: 1 - x[0]},
          {'type': 'ineq', "fun": lambda x: 1 - x[1]},
          {'type': 'ineq', "fun": lambda x: x[2] - 0.001},
          {'type': 'ineq', "fun": lambda x: 1 - x[2]},
          {'type': 'ineq', "fun": cf.const1},
          {'type': 'ineq', "fun": cf.const2},
          {'type': 'ineq', "fun": cf.const3},
          {'type': 'ineq', "fun": cf.const4},
          {'type': 'ineq', "fun": cf.const5},
          {'type': 'ineq', "fun": cf.const6})

res = sp.optimize.minimize(cf, x0, constraints=consts)

print(res)

xr = res['x']
print(cf.vals(xr))
print(cf.consts)

