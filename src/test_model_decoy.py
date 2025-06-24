# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 15:23:03 2023

@author: vil034
"""


import model_decoy as mod
import numpy as np
import distance_analysis as da
import matplotlib.pyplot as plt
import fine_analysis as fa

plt.close('all')

new_data = True
c = 299792458

if new_data:
    mu1 = 0.4
    mu2 = 0.5
    
    setup = mod.ModelDecoy(mu1, mu2)
    
    setup.set_background_rate(1e4)
    setup.set_qber(.01)
    setup.set_loss(30, 0)
    # setup.set_jitter(0)
    setup.set_pulse_rate(20)
    
    
    time_shift = 33.356
    start_time = 0
    total_protocol_time_sec = 0.1

    
    
    total_protocol_time = total_protocol_time_sec * 1e6
    ta, tb, qa, qb = setup.generate_time_data(time_shift, start_time, total_protocol_time_sec, verbose=True)
    
    
    
    
    n = 2
    offset_rough = time_shift - 0.05
    # offset_rough = 33.011
    dists, idx_a = da.get_distances_B_to_A(ta, tb, qa, qb, n, offset_rough, limit=50) 
    # dists = da.get_distances_A_to_B(ta, tb, qa, qb, n, offset_rough, limit=50) 
    
    
    bins1 = 500
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = ax.hist(dists, bins=bins1, fc='k', alpha=0.3)


    plt.plot(np.ones(15)*time_shift, np.arange(0,15)*10, 'r--')
    
    plt.show()
    
    # need to repeat this for other peaks if first fails
    offset_ind = np.argmax(n)
    
    offset = (bins[offset_ind] + bins[offset_ind+1])/2
    time_res = bins[1] -  bins[0]
    
    print("Offset: ", offset)
    print('Precision: ', np.abs(time_res))
    success = time_shift < bins[offset_ind+1] and time_shift > bins[offset_ind]
    print('Success: ', success)
    
    ta2, tb2, qa2, qb2 = fa.get_overlap(ta[idx_a], tb, qa[idx_a], qb, offset, time_res)
    
    np.save('alice_time_points', ta2)
    np.save('bob_time_points', tb2)
    np.save('alice_q_points', qa2)
    np.save('bob_q_points', qb2)


print('-------------- Fine Analysis')


# sorted(iterable, key=key, reverse=True)[:n]

alice_times = np.load('alice_time_points.npy')
bob_times = np.load('bob_time_points.npy')
alice_q = np.load('alice_q_points.npy')
bob_q = np.load('bob_q_points.npy')

diff = bob_times - alice_times
bins2 = 600

pfit = fa.fit_w_hist(diff, bins2)


idx = np.where(np.abs(diff - pfit[0]) < 2*np.abs(pfit[1]) )

qber = sum(alice_q[idx] != bob_q[idx])/sum(alice_q[idx] != 0)

# need tom compute qber from signal pulses only

plt.hist(diff[idx], fc='k', alpha=0.3)

print('Offset: ', pfit[0])
print('Precision: ', np.abs(pfit[1]))
print('Measured QBER: ', qber)

print('Pairs detected: ', len(idx[0]))



# sorted_indices = np.argsort(arr)[::-1]
# largest_indices = sorted_indices[:3]
# print(largest_indices)



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

