# -*- coding: utf-8 -*-


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


def P_fail_ES_2jitter(ra, rb, channel_loss, detector_loss, jitter, time):
    eta = 10**(-(channel_loss + detector_loss)/10)
    ra = ra*10**(-detector_loss/10)
    jitter = np.sqrt(2)*jitter
    return P_fail_2jitter(ra, rb, eta, jitter, time)

def P_fail_decoy_2jitter(ra, rb, channel_loss, detector_loss, jitter, time):
    eta = 10**(-(channel_loss + detector_loss)/10)
    return P_fail_2jitter(ra, rb, eta, jitter, time)


    
def P_fail_2jitter(ra, rb, eta, sig, T):
    dt = 2*sig
    time_filter_keep = 0.68
    # dt = 2*sig
    # time_filter_keep = 0.68
    # number of bins in the coarse guess
    # N0 = 500
    N0 = T/dt
    
    Ns = eta*ra*T*time_filter_keep
    mub = ra*rb*dt*T
    
    
    
    P_failure = N0*.5*sp.special.erfc((Ns)/np.sqrt(2*mub))
    # P_failure[P_failure > 1] = 1 
    if P_failure > 1:
        P_failure = 1
    return  P_failure



##################### Define parameters

jitter = 364e-12
jitter_es = jitter*np.sqrt(2)

# time = 0.1
sig_rate = 100e6
# sig_rate_es = 100e6

qber_dev = 0.05
background = 400



e_d = 0.05 # detector QBER
time_filter = 2

det_loss = 5 
time = 10 # total time of prototocol

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



#########################################
z = [0, 10 , 20, 30, 40, 50, 60]
loss = np.array([23.3694,   24.0985,   25.3455,   26.6543,   29.1963,   31.7907,   35.5599]) + 20
                



#####################################

confidence = 0.99

m = 0.6 # signal mean photon number
v = 0.2 # decoy mean photon number  - required v < m
pm = 0.6 # probability of sigal produced

x = [m, v, pm]
        
## NiGHT
qber_night_decoy = []
ns_night_decoy = []
signal_night_decoy = []

# ph_num_sig = []
# ph_num_decoy = []


qber_night_ES = []
signal_night_ES = []


pfail = []
pfail_es = []

channel['bg_rate'] = background
for i in loss:
    channel['channel_loss'] = i
    
    opt = optimizer.Optimizer(channel)

    opt.set_initial_conditions(x)
    opt.set_confidence(confidence)

    # opt.optimize2()

    ns, qber = opt.theoretical_vals(opt.x, opt.conf)
    
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
    
    pfail += [P_fail_decoy_2jitter(sig_rate, background, i, det_loss, jitter, time/10)]
    pfail_es += [P_fail_ES_2jitter(sig_rate, background, i, det_loss, jitter, time/10)]

 


## Second parameter
sig_rate = 400e6

channel = {
    'ra': sig_rate,
    'channel_loss': 0,
    'detector_loss': det_loss,
    'bg_rate': 0,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': qber_dev}

qber_night_decoy2 = []
ns_night_decoy2 = []
signal_night_decoy2 = []

# ph_num_sig = []
# ph_num_decoy = []


qber_night_ES2 = []
signal_night_ES2 = []


pfail2 = []
pfail_es2 = []

channel['bg_rate'] = background
for i in loss:
    channel['channel_loss'] = i
    
    opt = optimizer.Optimizer(channel)

    opt.set_initial_conditions(x)
    opt.set_confidence(confidence)

    # opt.optimize2()

    ns, qber = opt.theoretical_vals(opt.x, opt.conf)
    
    # print(ns, qber)
    # print('--------------------')
    
    qber_night_decoy2 += [qber]
    ns_night_decoy2 += [ns]
    signal_night_decoy2 += [opt.signal(x)]
    # ph_num_sig += [opt.x[0]]
    # ph_num_decoy += [opt.x[1]]
    # decoy_sig_prob += [opt.x[2]]
    
    qber_night_ES2 += [opt.QBER_ES()]
    signal_night_ES2 += [opt.signal_ES()]
    
    pfail2 += [P_fail_decoy_2jitter(sig_rate, background, i, det_loss, jitter, time/10)]
    pfail_es2 += [P_fail_ES_2jitter(sig_rate, background, i, det_loss, jitter, time/10)]




signal_night_ES = np.array(signal_night_ES)
signal_night_decoy = np.array(signal_night_decoy)
signal_night_ES2 = np.array(signal_night_ES2)
signal_night_decoy2 = np.array(signal_night_decoy2)
     
plt.figure()

dists = z


plt.plot(dists, qber_night_decoy, label='Decoy - 100MHz')
# plt.plot(dists, qber_night_ES2, label="ES")
plt.plot(dists, qber_night_decoy2, label='Decoy - 400MHz')
plt.plot(dists, qber_night_ES, label="ES")

plt.xlabel('zenith angle (deg)')
plt.ylabel('QBER')
plt.legend()
plt.ylim([0, 0.25])

plt.figure()
plt.plot(dists, signal_night_ES, label="ES")
plt.plot(dists, ns_night_decoy, label='Decoy - sp bound')
plt.plot(dists, signal_night_decoy, label='Decoy - signal')
plt.plot(dists, signal_night_decoy/2, 'r--')
plt.legend()


plt.figure()
plt.plot(dists[0:2], 1e12*2*jitter/np.sqrt(signal_night_decoy[0:2]/10), label="Decoy - 100MHz")
plt.plot(dists[0:4], 1e12*2*jitter/np.sqrt(signal_night_decoy2[0:4]/10), label="Decoy - 400MHz")

plt.plot(dists, 1e12*2*jitter/np.sqrt(signal_night_ES/20), label="ES - 100MHz")
plt.plot(dists[0:-1], 1e12*2*jitter/np.sqrt(signal_night_ES2[0:-1]/20), label="ES - 400MHz")


plt.xlabel('zenith angle (deg)')
plt.ylabel('PNT precision (picosec)')
plt.legend()

print(pfail)
print(pfail_es)


print(pfail2)
print(pfail_es2)
