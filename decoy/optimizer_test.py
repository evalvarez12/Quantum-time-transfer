# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-

import optimizer
import numpy as np
import matplotlib.pyplot as plt

ra = 20e6 # pulse rate
bg_rate = 500e3 # background rate
jitter = 364e-12 # total time jitter
e_d = 0.05 # detector QBER
time_filter = 2 # 2 sigma time filter on the coincidence window

channel_loss = 12.5 # in dB
detector_loss = 5 
eta_db = channel_loss + detector_loss # total loss
eta = 10**(-eta_db/10) # loss
time = 1 # total time of prototocol



# background rate - night
# bg_rate = 10e3

channel = {
    'ra': ra,
    'channel_loss': channel_loss,
    'detector_loss': detector_loss,
    'bg_rate': bg_rate,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': e_d}
opt = optimizer.Optimizer(channel)

m = 0.5 # signal mean photon number
v = 0.3 # decoy mean photon number  - required v < m


pm = 0.6 # probability of sigal produced

x = [m, v, pm]


opt.set_initial_conditions(x)
opt.set_confidence(0.99)
# opt.set_confidence(0.5)
print(opt.theoretical_vals_decoy(x, 0.99))
