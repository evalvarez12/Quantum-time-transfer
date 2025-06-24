# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 13:28:02 2023

@author: vil034
"""

import model as mod
import numpy as np
import distance_analysis as da
import matplotlib.pyplot as plt
import fine_analysis as fa
from scipy import special

plt.close('all')

new_data = True
c = 299792458

if new_data:
    setup = mod.Model()
    
    background_rate = 10000
    qber = 0.05 *(2/3)
    loss_channel = 20
    loss_detector = 5
    signal_rate = 20
    jitter = 500e-6
    
    setup.set_background_rate(background_rate)
    setup.set_qber(qber)
    setup.set_loss(loss_channel, loss_detector)
    setup.set_jitter(jitter)
    setup.set_signal_rate(signal_rate)
    
    
    time_shift = 33.356
    start_time = 0
    total_protocol_time_sec = .1

    
    
    total_protocol_time = total_protocol_time_sec * 1e6
    ta, tb, qa, qb = setup.generate_time_data(time_shift, start_time, total_protocol_time_sec, verbose=True)
    
    
    
    
    n = 4
    # offset_rough = time_shift - 0.05
    offset_rough = 33.011
    dists, idx_a = da.get_distances_B_to_A(ta, tb, qa, qb, n, offset_rough, limit=50) 
    # dists = da.get_distances_A_to_B(ta, tb, qa, qb, n, offset_rough, limit=50) 
    
    
    bins1 = int(np.ceil((max(dists) - min(dists))/(5*jitter)))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = ax.hist(dists, bins=bins1, fc='k', alpha=0.4)


   
    
    
    offset_ind = np.argmax(n)
    
    offset = (bins[offset_ind] + bins[offset_ind+1])/2
    time_res = bins[1] -  bins[0]
    
    plt.plot(np.ones(15)*time_shift, np.linspace(0,0.9*max(n), 15), 'r-', alpha=0.2, linewidth=4)
    
    plt.show()
    
    
    print("Offset: ", offset)
    print('Precision: ', np.abs(time_res))
    success = time_shift < bins[offset_ind+1] and time_shift > bins[offset_ind]
    print('Success: ', success)
    
    
    if success: 
        ta2, tb2, qa2, qb2 = fa.get_overlap(ta[idx_a], tb, qa[idx_a], qb, offset, time_res)
    
        np.save('alice_time_points', ta2)
        np.save('bob_time_points', tb2)
        np.save('alice_q_points', qa2)
        np.save('bob_q_points', qb2)

        

        print('-------------- Fine Analysis')




        alice_times = np.load('alice_time_points.npy')
        bob_times = np.load('bob_time_points.npy')
        alice_q = np.load('alice_q_points.npy')
        bob_q = np.load('bob_q_points.npy')

        print("Coincidences used: ", len(alice_times))

        diff = bob_times - alice_times
        bins2 = int(np.ceil((max(diff) - min(diff))/(jitter/20)))
        
        pfit = fa.fit_w_hist(diff, bins2)
        
        # time filter = 2 jitter 
        idx = np.where(np.abs(diff - pfit[0]) < np.abs(pfit[1]) )
        
        qber_measured = sum(alice_q[idx] != bob_q[idx])/len(alice_q[idx])
        
        
        plt.hist(diff[idx], fc='k', alpha=0.3)
        
        print('Offset: ', pfit[0])
        print('Precision: ', np.abs(pfit[1]))
        print('Measured QBER: ', qber_measured)
        
        print('Pairs used for QBER: ', len(idx[0]))
        
        
        print('----------------- Central limit theorem')
        
        alice_times2 = alice_times[idx]
        bob_times2 = bob_times[idx]
        diff2 = diff[idx]
        np.random.shuffle(diff2)
        
        n_cl = 2
        n_blocks = int(np.floor(len(diff2)/n_cl))
        
        blocks = np.zeros((n_cl, n_blocks))
        for i in range(n_blocks):
            blocks[:,i] = diff2[i*n_cl:(i+1)*n_cl]
            
        avg1 = np.average(blocks, 0)
        avg2 = np.average(avg1)
        
        print('Offset: ', avg2)
        print('Precision: ', np.std(avg1))
        print('Expected precision: ', np.abs(pfit[1])/np.sqrt(n_cl))
            
        print('Precision meters: ', np.std(avg1)*1e-6*c)
        
        print('REAL OFFSET: ', time_shift)
        print('WANTED PRECISION: 0.1 m',)

print('-------------------- Theory')
time_jitter = jitter
N = total_protocol_time/(time_jitter)
dt = total_protocol_time/N
# each bin has the width of time jitter - time_jitter is the time unit
Ns = signal_rate*total_protocol_time*(10**(-(loss_channel + 2*loss_detector)/10))
ra = signal_rate*total_protocol_time*(10**(-loss_detector/10))/N
rb = background_rate*total_protocol_time_sec/N
p = ra*rb
mu_b = N*p
# sig_b2 = N*p*(1-p)
# print('Avg noise: ', mu_b)
# print('Std noise: ', sig_b2)

P_failure = .5*special.erfc((Ns)/np.sqrt(2*mu_b))
print('Probability of large bin: ', P_failure)
print('Probability of failure: ', P_failure*N)    


n_background_in_peak = ra*N*background_rate*1e-6*2*time_jitter 
qber_est = (Ns*qber + n_background_in_peak/2)/(Ns + n_background_in_peak)
print('QBER: ', qber_est)


####### CC analysis again
# time2 = max(tb2 - ta2) - min(tb2 - ta2)

# offset2, time_res2 = cc.analysis(ta2, tb2, n, time2, True, True)


# print('----------- Fine analysis 1')

# ta2, tb2, qa2, qb2 = fa.get_overlap(tb, ta, qb, qa, offset, time_res)
# print('Data points: ', len(ta2))

# np.save('alice_time_points', ta2)
# np.save('bob_time_points', tb2)
# np.save('alice_q_points', qa2)
# np.save('bob_q_points', qb2)



# alice_times = np.load('alice_time_points.npy')
# bob_times = np.load('bob_time_points.npy')
# alice_q = np.load('alice_q_points.npy')
# bob_q = np.load('bob_q_points.npy')

# diff = alice_times - bob_times

# pfit = fa.fit_w_hist(diff, 200)

# qber = sum(alice_q != bob_q)/len(alice_q)


# print('Offset: ', pfit[1])
# print('Precision: ', pfit[2])
# print('Floor: ', pfit[3])
# print('Measured QBER: ', qber)

# print('----------- Fine analysis 2')


# # m1 = (diff - pfit[1]) >  - 2*pfit[2]
# # m2  = (diff - pfit[1]) < 2*pfit[2]

# # m = np.logical_and(m1, m2)

# m = np.where(np.abs(diff - pfit[1]) < pfit[2])

# diff2 = diff[m]
# qa3 = alice_q[m]
# qb3 = bob_q[m]

# # if len()

# print('Data points: ', len(diff2))



# pfit2 = fa.fit_w_hist(diff2, 200)
# qber = sum(qa3 != qb3)/len(qa3)

# print('Offset: ', pfit2[1])
# print('Precision: ', pfit2[2])
# print('Floor: ', pfit2[3])
# print('Measured QBER: ', qber)

# print('Mean: ', np.average(diff2))

# print('Std: ', np.std(diff2))

