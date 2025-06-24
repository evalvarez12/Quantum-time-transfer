# -*- coding: utf-8 -*-


import numpy as np
import optimizer
import matplotlib.pyplot as plt
import scipy as sp


# plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

def QBER_2sig(qber_dev, sig_rate, back_rate, eta, jitter, time):
    time_filter = 2*jitter
    time_filter_keep = 0.68
    
    Ns = time_filter_keep*eta*sig_rate*time
    Nb = sig_rate*back_rate*time*time_filter
    
    return (Ns*qber_dev + Nb/2)/(Ns + Nb)

##################### Define parameters

jitter = 364e-12
jitter_es = jitter*np.sqrt(2)

# time = 0.1
sig_rate = 20e6
qber_dev = 0.05
background_day = 429518
background_night = 10000


e_d = 0.05 # detector QBER
time_filter = 2

det_loss = 5 
time = 1 # total time of prototocol

plt.rcParams["font.family"] = "Times New Roman"

channel = {
    'ra': sig_rate,
    'channel_loss': 0,
    'detector_loss': det_loss,
    'bg_rate': 0,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': qber_dev}



# qber_day = QBER_2sig(qber_dev, sig_rate, background_day, loss_day, jitter_es, time)
# qber_night = QBER_2sig(qber_dev, sig_rate, background_night, loss_night, jitter_es, time)

# fig = plt.figure()
# plt.plot(dists, qber_day, label='Day', linewidth=3)
# plt.plot(dists, qber_night, label='Night', linewidth=3)
# plt.xlabel('Distance (m)', fontsize=15)
# plt.ylabel('QBER', fontsize=15)
# plt.xticks(fontsize=13, rotation=0)
# plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)
# plt.ylim([0.0, 0.27])

# fig.set_size_inches(18.5*.3, 10.5*.3)
# fig.savefig('qber_main.pdf', dpi=200)






#####################################


m = 0.55 # signal mean photon number
v = 0.1 # decoy mean photon number  - required v < m
pm = 0.5 # probability of sigal produced

x = [m, v, pm]
    
loss_night = np.load('loss_night.npy')
    
## NiGHT
qber_night_decoy = []
ns_night_decoy = []
ph_num_sig = []
ph_num_decoy = []
signal_night_decoy = []


qber_night_ES = []
signal_night_ES = []

channel['bg_rate'] = background_night
for i in loss_night:
    channel['channel_loss'] = i
    
    opt = optimizer.Optimizer(channel)

    opt.set_initial_conditions(x)
    opt.set_confidence(0.99)

    # opt.optimize2()

    # opt.x[0] = opt.x[0] + 0.2

    ns, qber = opt.theoretical_vals(opt.x, opt.conf)
    
    # print(ns, qber)
    # print('--------------------')
    
    qber_night_decoy += [qber]
    ns_night_decoy += [ns]
    ph_num_sig += [opt.x[0]]
    ph_num_decoy += [opt.x[1]]
    decoy_sig_prob += [opt.x[2]]
    


loss_night_decoy = np.array(loss_night_decoy)
eta = 10**(-(loss_night_decoy + det_loss)/10)
eta = eta*0.6826
ph_num_sig = np.array(ph_num_sig)    
ph_num_decoy = np.array(ph_num_decoy)
decoy_sig_prob = np.array(decoy_sig_prob)
signal_night = sig_rate*time*(opt.bg_prob + decoy_sig_prob*(1 - np.exp(-eta*ph_num_sig)) + (1-decoy_sig_prob)*(1 - np.exp(-eta*ph_num_decoy)))    

bg_prob_es = 2*jitter_es*background_night*10**(-det_loss/10)
coincidences_night = time * sig_rate * (loss_night*0.6826 + bg_prob_es)
      
plt.figure()

plt.plot(dists, qber_night, label="ES")
plt.plot(dists, qber_night_decoy, label='Decoy')
plt.legend()

plt.figure()
plt.plot(dists, coincidences_night, label="ES")
plt.plot(dists, ns_night_decoy, label='Decoy - sp bound')
plt.plot(dists, signal_night, label='Decoy - signal')
plt.plot(dists, signal_night/2, 'r--')
plt.legend()





### Day

m = 0.3 # signal mean photon number
v = 0.45 # decoy mean photon number  - required v < m
pm = 0.6 # probability of sigal produced

x = [m, v, pm]

qber_day_decoy = []
ns_day_decoy = []
ph_num_sig = []
ph_num_decoy = []
decoy_sig_prob = []

channel['bg_rate'] = background_day
for i in loss_day_decoy:
    channel['channel_loss'] = i
    
    opt = optimizer.Optimizer(channel)

    opt.set_initial_conditions(x)
    opt.set_confidence(0.99)

    # opt.optimize2()

    ns, qber = opt.theoretical_vals(opt.x, opt.conf)
    
    # print(ns, qber)
    # print('--------------------')
    
    qber_day_decoy += [qber]
    ns_day_decoy += [ns]
    ph_num_sig += [opt.x[0]]
    ph_num_decoy += [opt.x[1]]
    decoy_sig_prob += [opt.x[2]]
    


loss_day_decoy = np.array(loss_day_decoy)
eta = 10**(-(loss_day_decoy + det_loss)/10)
ph_num_sig = np.array(ph_num_sig)    
ph_num_decoy = np.array(ph_num_decoy)
decoy_sig_prob = np.array(decoy_sig_prob)
signal_day = sig_rate*time*(opt.bg_prob + decoy_sig_prob*(1 - np.exp(-eta*ph_num_sig)) + (1-decoy_sig_prob)*(1 - np.exp(-eta*ph_num_decoy)))     

bg_prob_es = 2*jitter_es*background_day*10**(-det_loss/10)
coincidences_day = time * sig_rate * (loss_day*0.6826 + bg_prob_es)

plt.figure()

plt.plot(dists, qber_day, label="ES")
plt.plot(dists, qber_day_decoy, label='Decoy')
plt.legend()

plt.ylim([0, .5] )

plt.figure()
plt.plot(dists, coincidences_day, label="ES")
plt.plot(dists, ns_day_decoy, label='Decoy - sp bound')
plt.plot(dists, signal_day, label='Decoy - signal')
plt.plot(dists, signal_day/2, 'r--')
plt.legend()


# plt.ylim([0, 500000] )