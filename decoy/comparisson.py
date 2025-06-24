# -*- coding: utf-8 -*-


import numpy as np
import optimizer
import matplotlib.pyplot as plt
import scipy as sp


plt.close('all')
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
    'e_d': qber_dev,
    'channel_fluctuation': 0.005}



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



#########################################
file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
file_night = 'data/HORIZONTAL_DISTANCE_NIGHT2_trx=15_L=%i_wl=0.81'

# dists = [7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000]
dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]


m_d = []
m_n = []


scattering_loss = 0.45 # db/km for 810

for d in dists:
    data_d = sp.io.loadmat(file_day %d)
    data_d = data_d['res']
    data_d = -10*np.log10(np.average(data_d))+scattering_loss*d*1e-3
    # data_day = -10*np.log10(data_day) 

    m_d += [data_d]
    
    data_n = sp.io.loadmat(file_night %d)
    data_n = data_n['res']
    data_n = -10*np.log10(np.average(data_n))+scattering_loss*d*1e-3
    # data_n = -10*np.log10(data_n)
    
    m_n += [data_n]
    
    
    
np.save('loss_night.npy', m_n)
np.save('loss_day.npy', m_d)



#####################################

confidence = 0.998

m = 0.4 # signal mean photon number
v = 0.1 # decoy mean photon number  - required v < m
pm = 0.7 # probability of sigal produced

x = [m, v, pm]
    
loss_night = np.load('loss_night.npy')
    
## NiGHT
qber_night_decoy = []
ns_night_decoy = []
signal_night_decoy = []

# ph_num_sig = []
# ph_num_decoy = []


qber_night_ES = []
signal_night_ES = []

channel['bg_rate'] = background_night
for i in loss_night:
    channel['channel_loss'] = i
    
    opt = optimizer.Optimizer(channel)

    opt.set_initial_conditions(x)
    opt.set_confidence(confidence)

    # opt.optimize2()

    ns, qber = opt.theoretical_vals_decoy(opt.x, opt.conf)
    
    # print(ns, qber)
    # print('--------------------')
    
    qber_night_decoy += [qber]
    ns_night_decoy += [ns]
    signal_night_decoy += [opt.signal(x)]
    # ph_num_sig += [opt.x[0]]
    # ph_num_decoy += [opt.x[1]]
    # decoy_sig_prob += [opt.x[2]]
    
    qber_night_ES += [opt.QBER_ES()]
    signal_night_ES += [opt.signal_ES()]
    

      
plt.figure()

dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]


plt.plot(dists, qber_night_ES, label="ES")
plt.plot(dists, qber_night_decoy, label='Decoy')
plt.legend()
plt.xlabel('Distance')
plt.ylabel('QBER')

plt.ylim([0, .25] )


plt.figure()
plt.plot(dists, signal_night_ES, label="ES")
plt.plot(dists, ns_night_decoy, label='Decoy - sp bound')
plt.plot(dists, signal_night_decoy, label='Decoy - signal')
plt.plot(dists, np.array(signal_night_decoy)/2, 'r--')
plt.legend()
plt.xlabel('Distance')
plt.ylabel('Coincidences')




### Day

# m = 0.4 # signal mean photon number
# v = 0.1 # decoy mean photon number  - required v < m
# pm = 0.6 # probability of sigal produced

# x = [m, v, pm]

# loss_day = np.load('loss_day.npy')

# qber_day_decoy = []
# ns_day_decoy = []
# signal_day_decoy = []

# # ph_num_sig = []
# # ph_num_decoy = []
# # decoy_sig_prob = []

# qber_day_ES = []
# signal_day_ES = []

# channel['bg_rate'] = background_day
# for i in loss_day:
#     channel['channel_loss'] = i
    
#     opt = optimizer.Optimizer(channel)

#     opt.set_initial_conditions(x)
#     opt.set_confidence(confidence)

#     # opt.optimize2()

#     ns, qber = opt.theoretical_vals(opt.x, opt.conf)
    
#     # print(ns, qber)
#     # print('--------------------')
    
#     qber_day_decoy += [qber]
#     ns_day_decoy += [ns]
#     signal_day_decoy += [opt.signal(x)]
    
    
#     # ph_num_sig += [opt.x[0]]
#     # ph_num_decoy += [opt.x[1]]
#     # decoy_sig_prob += [opt.x[2]]
    
#     signal_day_ES += [opt.signal_ES()]
#     qber_day_ES += [opt.QBER_ES()]    
    
    
#     # print('Nb: ', opt.bg_prob_ES*channel['ra']*channel['time'])
#     # print(opt.eta_ES)



# plt.figure()

# plt.plot(dists, qber_day_ES, label="ES")
# plt.plot(dists, qber_day_decoy, label='Decoy')
# plt.legend()
# plt.xlabel('Distance')
# plt.ylabel('QBER')

# plt.ylim([0, .25] )

# plt.figure()
# plt.plot(dists, signal_day_ES, label="ES")
# plt.plot(dists, ns_day_decoy, label='Decoy - sp bound')
# plt.plot(dists, signal_day_decoy, label='Decoy - signal')
# plt.plot(dists, np.array(signal_day_decoy)/2, 'r--')
# plt.legend()
# plt.xlabel('Distance')
# plt.ylabel('Coincidences')

################################################

channel2 = {
    'ra': sig_rate,
    'channel_loss': 0,
    'detector_loss': det_loss,
    'bg_rate': 0,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': 5, 
    'e_d': qber_dev,
    'channel_fluctuation': 0.005}


channel3 = {
    'ra': sig_rate,
    'channel_loss': 0,
    'detector_loss': det_loss,
    'bg_rate': 0,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': 10, 
    'e_d': qber_dev,
    'channel_fluctuation': 0.005}


confidence = 0.998

m = 0.4 # signal mean photon number
v = 0.1 # decoy mean photon number  - required v < m
pm = 0.7 # probability of sigal produced

x = [m, v, pm]
    
loss_night = np.load('loss_night.npy')
    
## NiGHT
# qber_night_decoy = []
ns_night_decoy = []
signal_night_decoy = []


qber_night_decoy2 = []
qber_night_decoy3 = []
# ph_num_sig = []
# ph_num_decoy = []



channel['bg_rate'] = background_night
for i in loss_night:
    channel2['channel_loss'] = i
    channel3['channel_loss'] = i
    
    opt2 = optimizer.Optimizer(channel2)
    opt3 = optimizer.Optimizer(channel3)
    
    opt2.set_initial_conditions(x)
    opt2.set_confidence(confidence)
    
    opt3.set_initial_conditions(x)
    opt3.set_confidence(confidence)

    # opt.optimize2()


    ns2, qber3 = opt2.theoretical_vals_decoy(opt2.x, opt2.conf)
    ns3, qber2 = opt3.theoretical_vals_decoy(opt3.x, opt3.conf)

    
    qber_night_decoy3 += [qber2]
    qber_night_decoy2 += [qber3]

    
    

      
plt.figure()

dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]


plt.plot(dists, qber_night_ES, label="ES")
plt.plot(dists, qber_night_decoy, label='Decoy - 1 sec')

plt.plot(dists, qber_night_decoy2, label='Decoy - 5 sec')
plt.plot(dists, qber_night_decoy3, label='Decoy - 10 sec')


plt.legend()
plt.xlabel('Distance')
plt.ylabel('QBER')

plt.ylim([0, .25] )





