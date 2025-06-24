# -*- coding: utf-8 -*-

import sys
sys.path.append('./src')
import numpy as np
import TimeTagger as tt
import distance_analysis as da
import matplotlib.pyplot as plt
import fine_analysis as fa
import scipy as sp
import os

import glob


plt.close('all')

# # All files and directories ending with .txt and that don't begin with a dot:
# print(glob.glob("/home/adam/*.txt")) 
# # All files and directories ending with .txt with depth of 2 folders, ignoring names beginning with a dot:
# print(glob.glob("/home/adam/*/*.txt")) 


# os.chdir(r'Y:\work\E91 Sources\Correlated Source\24-07-2024\10nm_BPF_780nm_arm\coincidence_win_scan\100_ps')
os.chdir(r'Y:\work\E91 Sources\Correlated Source\24-07-2024\10nm_BPF_780nm_arm\coincidence_win_scan\1000_ps')


print(os.getcwd())


files = glob.glob('*.ttbin')
files1 = glob.glob('*.1.ttbin')

files = [i for i in files if i not in files1]

a_range = []
b_range = []
a_len = []
b_len = []

dists_total = np.array([])

windows = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
coincidences = np.array([0]*10) 
tagged = np.array([0]*10)

noise_perc = 0 # noise percentage to be added to tb


if True: 
    for fi in files:
        
    # fi = files[-20]
        print('------> ', fi)
        fr = tt.FileReader(fi)
        
        delay_ch2 = int(fi.split('_')[5]) 
        print('Delay: ', delay_ch2)
            
        stream = fr.getData(1000000)        
        
        ts = stream.getTimestamps()
        ch = stream.getChannels()
        
        ta = np.array(ts[ch == 1])
        tb = np.array(ts[ch == 2])
        
        tb = tb - delay_ch2
        
        print('ta size: ', len(ta))
        a_len += [len(ta)]
        min_ta = min(ta) 
        max_ta = max(ta)
        a_range += [[min_ta, max_ta]]
        print(min_ta, max_ta)
        
        print('tb size: ', len(tb))
        b_len += [len(tb)]
        min_tb = min(tb) 
        max_tb = max(tb)
        b_range += [[min_tb, max_tb]]
        print(min_tb, max_tb)
        
        print('Min max diffs: ', min(ta) - min(tb), max(ta) - max(tb))
        
        # Add noise to tb
        noise_len = int(noise_perc*len(tb))
        noise = np.random.randint(min_tb, max_tb, noise_len, dtype='int64')
        tags_b = np.array([0]*len(tb))
        tags_noise = np.array([1]*noise_len)
        
        tb = np.concatenate((tb, noise))
        sort_ind = np.argsort(tb)
        tags = np.concatenate((tags_b, tags_noise))
        
        tb = tb[sort_ind]
        tags = tags[sort_ind]
    
    
    
        offset_rough = 0 # rough offset given by the GPS - approx 90 meters off 
        n = 100 # this parameter sets the number of timestamps to each side of the rough_offset to be used to find the peak
        window = 5000  # in ps
        
        # dists = da.get_distances_B_to_A_simple(ta, tb, window) # find the distances and the indexes of the timestamps involved in Alice data
        dists, _ = da.get_distances_B_to_A_noq(ta, tb, n, offset_rough, limit=window) # find the distances and the indexes of the timestamps involved in Alice data
        
        
        # print('Dists : ', dists)
        
        dists_total = np.append(dists_total, dists)
        
        
        for j in range(len(windows)):
            c,tg = da.get_coincidences_A_to_B_wtags(ta, tb, tags, offset_rough, windows[j])
            coincidences[j] += c
            tagged[j] += tg
        
    print('-----------------------------------------------')

        
    print('Data points: ', len(dists_total))
    print('Total time (sec):')
    time = (a_range[-1][-1] - a_range[0][0])*1e-12
    print(a_range[0][0]*1e-12, a_range[-1][-1]*1e-12, time)
    print(b_range[0][0]*1e-12, b_range[-1][-1]*1e-12, (b_range[-1][-1] - b_range[0][0])*1e-12)




    print('a data points: ', sum(a_len))
    a_rate = sum(a_len)/(time)
    print('a rate (sec): ', a_rate)


    print('b data points: ', sum(b_len))
    b_rate = sum(b_len)/(time)
    b_rate2 = 1e12*sum(b_len)/(b_range[-1][-1] - b_range[0][0])
    
    print('b rate (sec): ', b_rate)
    print('b rate 2 (sec): ', b_rate2)


    coincidences_rate = coincidences/time 
    print(coincidences_rate)



    
    np.save('dists.npy', dists_total)

dist_data = np.load('dists.npy')


# dist_data = dist_data[abs(dist_data) < 1500] 

Nbins = 100
fig = plt.figure()
fig.set_size_inches(18.5*.6, 10.5*.6)
(n, bins, patches) = plt.hist(dist_data, Nbins, ec='k', alpha=0.7)
plt.xlabel('time difference (ps)')

bin_width = bins[1] - bins[0]

print('Bin width (ps): ', bin_width)
# print('acc rate per bin (sec): ', a_rate*b_rate*bin_width*1e-12)
# print('# acc coincidences in bin: ', a_rate*b_rate*bin_width*1e-12*(a_range[-1][-1] - a_range[0][0])*1e-12)

# Gaussian function for fitting 
def Gaussian(t, amp, mean, std, floor):
    return amp*np.exp(-1/2*((t-mean)/std)**2) + floor


my = np.max(n)
param_init = [my, 500, 500, 0]
# print(param_init)
    
popt, pcov = sp.optimize.curve_fit(Gaussian, bins[1:], n, p0=param_init, maxfev=2000)

        # plt.figure()
        # plt.plot(x,y,color='black',label='Data',alpha=0.6)
plt.plot(bins[1:], Gaussian(bins[1:], popt[0], popt[1], popt[2], popt[3]) , label='Fitted Gaussian',color='red',linestyle='dotted', linewidth=3)
plt.xlabel(r"$\Delta t$", fontsize=13)
plt.ylabel(r"Coincidences",  fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)


print('---------- Gaussian fit')
print('Amplitude : ', popt[0])
print('Mean : ', popt[1])
print('Std : ', popt[2])
print('Floor : ', popt[3])


qber = tagged/(2*coincidences)
print('coincidences :', coincidences)
print('tagged :', tagged)


fig = plt.figure()
fig.set_size_inches(18.5*.6, 10.5*.6)
plt.plot(windows, coincidences, 'o')
plt.xlabel(r"Coincidence window", fontsize=14)
plt.ylabel(r"Coincidences",  fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)


fig = plt.figure()
fig.set_size_inches(18.5*.6, 10.5*.6)
plt.plot(windows, qber, 'o')
plt.xlabel(r"Coincidence window", fontsize=13)
plt.ylabel(r"QBER",  fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)


print('-----')
accidentals = a_rate*b_rate*bin_width*1e-12*time

background_rate = np.sqrt(popt[3]/(bin_width*1e-12*time)) 

print(f'Accidentals = {accidentals}')
print(f'Background_rate = {background_rate}')

loss_a = coincidences_rate/a_rate

loss_a =-10*np.log10(loss_a)

loss_b = coincidences_rate/b_rate

loss_b = -10*np.log10(loss_b)

pair_generation_rate = a_rate*b_rate / coincidences_rate


fig = plt.figure()
fig.set_size_inches(18.5*.6, 10.5*.6)
plt.plot(windows, loss_a, 'o', label='A')
plt.plot(windows, loss_b, 'o', label='B')
plt.xlabel(r"Coincidence window", fontsize=13)
plt.ylabel(r"Loss (dB)",  fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend()

fig = plt.figure()
fig.set_size_inches(18.5*.6, 10.5*.6)
plt.plot(windows, pair_generation_rate, 'o')
plt.xlabel(r"Coincidences window", fontsize=13)
plt.ylabel(r"Pair generation rate (Hz)",  fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)



