# -*- coding: utf-8 -*-

import sys
sys.path.append('./src')
import numpy as np
import TimeTagger as tt
import distance_analysis as da
import matplotlib.pyplot as plt
import fine_analysis as fa
import scipy as sp
plt.close('all')








# data_1 = r'Y:\work\E91 Sources\Correlated Source\10-07-2024\TT_data.1.ttbin'
# data_2 = r'Y:\work\E91 Sources\Correlated Source\10-07-2024\TT_data.2.ttbin'
# data_3 = r'Y:\work\E91 Sources\Correlated Source\10-07-2024\TT_data.3.ttbin'

# data_1 = r'Y:\work\E91 Sources\Correlated Source\15-07-2024\uncoated fiber_30.5C ppktp\TT_chan_1_2.1.ttbin'


######## 24/07 data
# data_1 = r'Y:\work\E91 Sources\Correlated Source\24-07-2024\10nm_BPF_780nm_arm\coincidence_win_scan\100_ps\TT_data_100_10_0_0_2024-07-24 16-01-59.ttbin'


data_1 = r'Y:\work\E91 Sources\Correlated Source\24-07-2024\10nm_BPF_780nm_arm\coincidence_win_scan\100_ps\TT_data_100_10_0_-500_2024-07-24 16-01-32.ttbin'


TT_data_100_10_0_0_2024-07-24 16-01-59

TT_data_100_10_0_0_2024-07-24 16-02-05



TT_data_100_10_0_-100_2024-07-24 16-01-53




# Read the time 
fr1 = tt.FileReader(data_1)
# fr2 = tt.FileReader(data_2)


# n_samples = 10000000
# ta = frA.getData(n_samples)
# tb = frB.getData(n_samples)

n_samples = 500000
n_it = 10
dists_total = np.array([])

a_range = []
b_range = []
a_len = []
b_len = []
# stream1 = fr1.getData(n_samples*150)


if True:
    i = 0
    while i < n_it and fr1.hasData():
        stream1 = fr1.getData(n_samples)
        # tb = frB.getData(n_samples)
        
        
        
        
        t1 = stream1.getTimestamps()
        ch1 = stream1.getChannels()
        
        ta = t1[ch1 == 1]
        tb = t1[ch1 == 2]
        
        print('ta size: ', len(ta))
        a_len += [len(ta)]
        a_range += [[min(ta), max(ta)]]
        print(min(ta), max(ta))
        
        print('tb size: ', len(tb))
        b_len += [len(tb)]
        b_range += [[min(tb), max(tb)]]
        print(min(tb), max(tb))
        
        print('Min max diffs: ', min(ta) - min(tb), max(ta) - max(tb))
        
        
        # Cross correlation analysis
        
        # cc = sp.signal.correlate(ta, tb, 'same')
        
        # plt.hist(cc, 100, ec='k', alpha=0.5)
        
        
        
        
        
        
        # STEP 1: Find the peak using the rough estimate from the GPS
        offset_rough = 0 # rough offset given by the GPS - approx 90 meters off 
        n = 100 # this parameter sets the number of timestamps to each side of the rough_offset to be used to find the peak
        window = 50000  # 2 nanoseconds
        
        # dists = da.get_distances_B_to_A_simple(ta, tb, window) # find the distances and the indexes of the timestamps involved in Alice data
        dists, _ = da.get_distances_B_to_A_noq(ta, tb, n, offset_rough, limit=window) # find the distances and the indexes of the timestamps involved in Alice data
        
        print('------> ', i)
        # print('Dists : ', dists)
        
        dists_total = np.append(dists_total, dists)
        
        i += 1
        

    
    np.save('dists.npy', dists_total)


dist_data = np.load('dists.npy')


Nbins = 500
fig = plt.figure()
(n, bins, patches) = plt.hist(dist_data, Nbins, ec='k', alpha=0.5)
plt.xlabel('time difference (ps)')

bin_width = bins[1] - bins[0]

print('Bin width (ps): ', bin_width)

print('Data points: ', len(dist_data))
print('Total time (ms):')
print(a_range[0][0]*1e-9, a_range[-1][-1]*1e-9, (a_range[-1][-1] - a_range[0][0])*1e-9)
print(b_range[0][0]*1e-9, b_range[-1][-1]*1e-9, (b_range[-1][-1] - b_range[0][0])*1e-9)


print('a data points: ', sum(a_len))
a_rate = 1e12*sum(a_len)/(a_range[-1][-1] - a_range[0][0])
print('a rate (sec): ', a_rate)


print('b data points: ', sum(b_len))
b_rate = 1e12*sum(b_len)/(b_range[-1][-1] - b_range[0][0])
print('b rate (sec): ', b_rate)

print('acc rate per bin (sec): ', a_rate*b_rate*bin_width*1e-12)
print('# acc coincidences in bin: ', a_rate*b_rate*bin_width*1e-12*(a_range[-1][-1] - a_range[0][0])*1e-12)
# Plot a histogram to find the peak
# at this point the peak width is smaller that the bin width - we a doing a coarse search


# 120 ps hist bin






# fig = plt.figure()
# plt.hist(dist_data, 600, ec='k', alpha=0.5)


# ax = fig.add_subplot(111)

# offset_ind = np.argmax(n)

# offset = (bins[offset_ind] + bins[offset_ind+1])/2
# time_res = bins[1] -  bins[0]

# # Highlight the largest bin where the peak will be
# plt.plot(np.ones(15)*0, np.linspace(0,max(n), 15), 'r-', alpha=0.2, linewidth=4)
# plt.title('Rough exploration to find the peak')
# plt.show()