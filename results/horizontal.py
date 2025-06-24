# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import optimizer 
import interpolator

def avg_std(arr):
    return np.average(arr), np.std(arr)


def log_norm(x, m, sig):
    return np.exp(-np.log(x-m)**2/(2*sig**2))/(x*sig*np.sqrt(2*np.pi))

def norm(x, m, sig):
    return np.exp(-(x-m)**2/(2*sig**2))/(sig*np.sqrt(2*np.pi))

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


channels = []

for d in dists:
    data_d = sp.io.loadmat(file_day %d)
    data_d = data_d['res']
    data_d = -10*np.log10(data_d)+scattering_loss*d*1e-3
    # data_day = -10*np.log10(data_day) 

    m_d += [np.average(data_d)]
    std_d += [np.std(data_d)]
    
    data_n = sp.io.loadmat(file_night %d)
    data_n = data_n['res']
    data_n = -10*np.log10(data_n)+scattering_loss*d*1e-3
    # data_n = -10*np.log10(data_n)
    
    m_n += [np.average(data_n)]
    std_n += [np.std(data_n)]





    # plt.figure()
    # y, x = np.histogram(data_d, 100, density=1)
    # plt.plot(x[1:],y)

    

    # p0 = [np.mean(data_d), np.std(data_d)]
    # res = sp.optimize.curve_fit(norm, x[1:], y, p0)
    # res = sp.optimize.curve_fit(log_norm, x[:-1], y, p0)

    # p1, p2 = res[0]

    # plt.plot(x, log_norm(x, p1, p2), '--')
    # plt.plot(x, norm(x, p1, p2), '--')

    # plt.plot(x, log_norm(x, np.mean(data_d), np.std(data_d)), '--')


    # ts = np.linspace(1, 100, 5000)
    # ps = norm(ts, p1,p2)*(ts[1]-ts[0])
    # channel_loss = lambda N : np.random.choice(ts, size=N, p=ps)

    # channels += [channel_loss]





    
# plt.errorbar(distsi, m_d, std_d)
# plt.errorbar(distsi, m_n, std_n)

fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.plot(dists, m_d, linewidth=3, label='Day')
plt.plot(dists, m_n, linewidth=3, label='Night')

plt.xlabel('Distance (m)',fontsize=13)
plt.ylabel('Loss (dB)', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)


plt.legend()
# fig.savefig('loss.pdf', dpi=200)



# data_d = 10**(-data_d/10)




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

sigs_d = []
ns_d = []
qs_d = []


sigs_es_d = []
qs_es_d = []

sigs_n = []
ns_n = []
qs_n = []

sigs_es_n = []
qs_es_n = []

for d in dists:

    data_d = sp.io.loadmat(file_day %d)
    data_d = data_d['res'].transpose()[0]
    data_d = -10*np.log10(data_d)+scattering_loss*d*1e-3

    interpol = interpolator.Interpolator(data_d, 60)
    channel['channel_loss'] = interpol.get_sampler()


    channel['bg_rate'] = bg_day
    opt = optimizer.Optimizer(channel)

    s, n, q = opt.theoretical_vals_decoy(x, confidence)
    
    sigs_d += [s]
    ns_d += [n]
    qs_d += [q]
    
    s, q = opt.theoretical_vals_ES(confidence)
    
    sigs_es_d += [s]
    qs_es_d += [q]
    
    data_n = sp.io.loadmat(file_night %d)
    data_n = data_n['res'].transpose()[0]
    data_n = -10*np.log10(data_n)+scattering_loss*d*1e-3

    interpol = interpolator.Interpolator(data_n, 60)
    channel['channel_loss'] = interpol.get_sampler()

    channel['bg_rate'] = bg_night
    opt = optimizer.Optimizer(channel)

    # s, n, q = opt.theoretical_vals_decoy(x, confidence)
    s, n, q = opt.optimize(x, confidence)
    
    sigs_n += [s]
    ns_n += [n]
    qs_n += [q]
    
    s, q = opt.theoretical_vals_ES(confidence)
    
    sigs_es_n += [s]
    qs_es_n += [q]


sigs_d = np.array(sigs_d)
sigs_n = np.array(sigs_n)

plt.figure()
plt.plot(dists, sigs_d)
plt.plot(dists, ns_d)
plt.plot(dists, sigs_d/2, 'r--')

plt.plot(dists, sigs_es_d, 'o')

plt.figure()
plt.plot(dists, qs_d)
plt.plot(dists, qs_es_d)
plt.ylim([0, 0.25])


plt.figure()
plt.plot(dists, sigs_n)
plt.plot(dists, ns_n)
plt.plot(dists, sigs_n/2, 'r--')

plt.plot(dists, sigs_es_n, 'o')


plt.figure()
plt.plot(dists, qs_n)
plt.plot(dists, qs_es_n)
plt.ylim([0, 0.25])






