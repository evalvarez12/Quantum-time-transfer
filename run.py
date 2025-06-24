# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:20:06 2023

"""
import sys
sys.path.append('./src')

import model as mod
import numpy as np
import distance_analysis as da
import matplotlib.pyplot as plt
import fine_analysis as fa
import scipy as sp

plt.close('all')

tags = True
setup = mod.Model(tags=tags)

# NOTE: Time units are in microsecod
background_rate_sec = 37000 # background counts per second 
qber = 0.05 # device qber not considering backgroud radiation
loss_channel = 20 # loss of the channel - only afects Bob
loss_detector = 5 # loss of the detection system - affects Alice and Bob
signal_rate = .05 # average entangled pair generation - 20 MHz
jitter = 364e-6
# jitter = 10e-6 

background_rate = background_rate_sec*1e-6
setup.set_background_rate(background_rate)
setup.set_qber(qber)
setup.set_loss(loss_channel, loss_detector)
setup.set_jitter(jitter)
setup.set_signal_rate(signal_rate)

setup.set_detector_dead_time(0)

time_shift = 33.356 # clock offset to be found in microseconds
start_time = 0 # zero position of the Alice clock
total_protocol_time_sec = .1 # total duration of the protocol in seconds
total_protocol_time = total_protocol_time_sec * 1e6 # protocol time microsec

# STEP 1: Generate data using a monte carlo simulation - for simplicity only 1 basis is used on the polarization
# ta = timing data Alice
# tb = timing data Bob - includes background radiation
# qa = quantum states prepared by Alice
# qb = quantum states measured by Bob - includes background radiation
ta, tb, qa, qb = setup.generate_time_data(time_shift, start_time, total_protocol_time_sec, verbose=True)


# STEP 2: Find the peak using the rough estimate from the GPS
offset_rough = 33.011 # rough offset given by the GPS - approx 90 meters off 
n = 8 # this parameter sets the number of timestamps to each side of the rough_offset to be used to find the peak
dists, idx_a = da.get_distances_B_to_A(ta, tb, qa, qb, n, offset_rough, limit=50) # find the distances and the indexes of the timestamps involved in Alice data



# Plot a histogram to find the peak
# at this point the peak width is smaller that the bin width - we a doing a coarse search
bins1 = int(np.ceil((max(dists) - min(dists))/(5*jitter)))

fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(dists, bins=bins1, fc='k', alpha=0.4)

offset_ind = np.argmax(n)

offset = (bins[offset_ind] + bins[offset_ind+1])/2
time_res = bins[1] -  bins[0]

# Highlight the largest bin where the peak will be
plt.plot(np.ones(15)*time_shift, np.linspace(0,max(n), 15), 'r-', alpha=0.2, linewidth=4)
plt.title('Rough exploration to find the peak')
plt.show()


print("Offset: ", offset)
print('Rough precision: ', np.abs(time_res))
success = time_shift < bins[offset_ind+1] and time_shift > bins[offset_ind]
print('Success: ', success)


# If the biggest bin does not contain the true offset abort 
if not success: 
    sys.exit("Incorrect time bin selected")

# Step 3: Identify all the coincidences that are in the peak
print('-------------- Fine Analysis')

# Find all the data inside the biggest bin
alice_times, bob_times, alice_q, bob_q = fa.get_overlap(ta[idx_a], tb, qa[idx_a], qb, offset, time_res)

if tags:
    tt = [np.where(tb == x)[0][0] for x in bob_times]
    tags1 = setup.tags_data[tt]


# Second histogram to find peak
diff = bob_times - alice_times

print('Num coincidences inside big bin: ', len(diff))
print('Avg of coincidence diff: ', np.mean(diff))
print('Std of diff: ', np.std(diff))

bins2 = int(np.ceil((max(diff) - min(diff))/(jitter/20)))

pfit = fa.fit_w_hist(diff, bins2)

# With the peak find use only the data points that are less that one std away from the mean - this time filter is good for QBER 
idx = np.where(np.abs(diff - pfit[0]) < np.abs(pfit[1]) )


# plot shaded are of the only coincidences used
# plt.hist(diff[idx], fc='k', alpha=0.3)

print("Coincidences in the peak: ", len(idx[0]))
print('Offset: ', pfit[0])
print('Precision raw: ', np.abs(pfit[1]))

print('Gauss amplitude: ', pfit[2])
print('Noise floor: ', pfit[3])

# print('Pairs used for QBER: ', len(idx[0]))

print('----------------- Security')
qber_limit = 0.25
qber_data = alice_q[idx] != bob_q[idx]
qber_measured = np.average(qber_data)
qber_std = np.std(qber_data)/np.sqrt(len(idx[0])/2) # number of data points divided over 2 to account for basis mismatch

confidence = 1-sp.special.erfc((qber_limit-qber_measured)/(np.sqrt(2)*qber_std))/2 
print('Measured QBER: ', qber_measured)
print('Sigma: ', qber_std)
print('Sigmas distance: ', (qber_limit - qber_measured)/qber_std)
print('Confidence: ', confidence)


print('----------------- Central limit theorem')

# Extract the data points of the coincidences used
alice_times2 = alice_times[idx]
bob_times2 = bob_times[idx]
diff2 = diff[idx]
np.random.shuffle(diff2)

if tags:
    tags2 = tags1[idx]

# Number of blocks used with the central limit theorem
n_cl = 2
n_blocks = int(np.floor(len(diff2)/n_cl))

blocks = np.zeros((n_cl, n_blocks))
for i in range(n_blocks):
    blocks[:,i] = diff2[i*n_cl:(i+1)*n_cl]
    
avg1 = np.average(blocks, 0)
avg2 = np.average(avg1)

print('Offset: ', avg2)
print('REAL OFFSET: ', time_shift)

print('Improved precision: ', np.std(avg1))
print('Expected precision: ', np.abs(pfit[1])/np.sqrt(n_cl))
    
print('Precision meters: ', np.std(avg1)*1e-6*sp.constants.c)

print('WANTED PRECISION: 0.1 m',)

print('----------------- Central limit theorem - simple')
print('Offset: ', np.average(diff2))

print('Improved precision: ', 2*np.std(diff2)/len(diff2))
    
print('Precision meters: ', 2*np.std(diff2)/len(diff2)*1e-6*sp.constants.c)



print('-------------------- Theory')
# Final step compare with theory to se if results are similar to what is expected
# Results might not be perfect as theory does not accounts for small effects such as detector dead time

# each bin has the width of time jitter - jitter is the time unit
total_jitter = np.sqrt(2)*jitter
N = total_protocol_time/(total_jitter)
dt = total_protocol_time/N

Ns = signal_rate*total_protocol_time*(10**(-(loss_channel + 2*loss_detector)/10)) # number of coincidences
ra = signal_rate*total_protocol_time*(10**(-loss_detector/10))/N # number of Alice data points
rb = background_rate*total_protocol_time/N # Number of background photons for Bob
p = ra*rb
mu_b = N*p
# sig_b2 = N*p*(1-p)
# print('Avg noise: ', mu_b)
# print('Std noise: ', sig_b2)

# Probability of misidentifiying the peak due to low SNR
P_failure = .5*sp.special.erfc((Ns)/np.sqrt(2*mu_b))
print('Probability of large bin: ', P_failure)
print('Probability of failure: ', P_failure*N)    

# QBER estimation based on statistics
print('---> 4 sig time filter')
print('Num of coincidences: ', Ns)
Nb = ra*N*background_rate*4*total_jitter 
print('Num background photons misidentified as signal: ', Nb)
qber_est = (Ns*qber + Nb/2)/(Ns + Nb)
print('QBER: ', qber_est)

print('---> 2 sig time filter')
Ns2  = 0.68*Ns
print('Num of coincidences: ', Ns2)
Nb2 = Nb/2
print('Num background photons misidentified as signal: ', Nb2)
print('QBER - 2 sig: ', (Ns2*qber + Nb2/2)/(Ns2 + Nb2))
# Estimated QBER always a bit lower than measured one for some reason - needs a closer look


if tags:
    # Tags are used on each timestamp tp tagg each one as signal or background
    # this allows for a more close look into how many background get into the Gauss peak
    print('-------------------- Tags')
    print('Total tags: ', len(tags2))
    print('Background tagged: ', np.sum(tags2))
