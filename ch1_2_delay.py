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
from tqdm import tqdm

import glob


plt.close('all')

# # All files and directories ending with .txt and that don't begin with a dot:
# print(glob.glob("/home/adam/*.txt")) 
# # All files and directories ending with .txt with depth of 2 folders, ignoring names beginning with a dot:
# print(glob.glob("/home/adam/*/*.txt")) 


os.chdir(r'C:\Users\vil034\OneDrive - CSIRO\Documents\laser data\aft_Ch3_delay_adjust\500_ps')


print(os.getcwd())


files = glob.glob('*.ttbin')
files1 = glob.glob('*.1.ttbin')

files = [i for i in files if i not in files1]

a_range = []
b_range = []
t1_len = []
b_len = []

dists_total = np.array([])

coincidence_window = 500


files = [files[202]]
# files = files[200:205]
# print(files)

streamSize= 500000

offset_rough = 0 # rough offset given by the GPS - approx 90 meters off 
n = 50 # this parameter sets the number of timestamps to each side of the rough_offset to be used to find the peak
window = 2000  # in ps

if True:    
    
    for fi in files:
        while True:
            print('File: ', fi)
            fr = tt.FileReader(fi)
            
            # delay_ch2 = int(fi.split('_')[5]) 
            # print('Delay: ', delay_ch2)
                
            stream = fr.getData(streamSize)        
            
            ts = stream.getTimestamps()
            ch = stream.getChannels()
            
            
            if not len(ts):
                break
            
            t1 = np.array(ts[ch == 1])
            tb = np.array(ts[ch == 2])
            # tb = np.array(ts[ch == 3])
            
            
            print('t1 size: ', len(t1))
            t1_len += [len(t1)]
            min_t1 = min(t1) 
            max_t1 = max(t1)
            a_range += [[min_t1, max_t1]]
            print('time interval : ', 1e-12*(max_t1 - min_t1))
            
            
            print('tb size: ', len(tb))
            b_len += [len(tb)]
            min_tb = min(tb) 
            max_tb = max(tb)
            b_range += [[min_tb, max_tb]]
            print('time interval : ', 1e-12*(max_tb - min_tb))
            
            
            # dists = da.get_distances_B_to_A_simple(ta, tb, window) # find the distances and the indexes of the timestamps involved in Alice data
            dists, _ = da.get_distances_B_to_A_noq(t1, tb, n, offset_rough, limit=window) # find the distances and the indexes of the timestamps involved in Alice data
            
            
            # print('Dists : ', dists)
            
            dists_total = np.append(dists_total, dists)
            break



    print('Data points: ', len(dists_total))
    print('Total time (sec):')
    time = (a_range[-1][-1] - a_range[0][0])*1e-12
    print(a_range[0][0]*1e-12, a_range[-1][-1]*1e-12, time)
    print(b_range[0][0]*1e-12, b_range[-1][-1]*1e-12, (b_range[-1][-1] - b_range[0][0])*1e-12)
    
    
    print('a data points: ', sum(t1_len))
    t1_rate = sum(t1_len)/(time)
    print('a rate (sec): ', t1_rate)
    
    
    print('b data points: ', sum(b_len))
    b_rate = sum(b_len)/(time)
    b_rate2 = 1e12*sum(b_len)/(b_range[-1][-1] - b_range[0][0])
    
    print('b rate (sec): ', b_rate)
    print('b rate 2 (sec): ', b_rate2)        
    
    np.save('dists1.npy', dists_total)


# print('g2 = ', (np.average(n1*n2))/(np.average(n1)*np.average(n2)))
print('............. Distances 1-2')
dist_data = np.load('dists1.npy')

     
Nbins = 200
fig = plt.figure()
fig.set_size_inches(18.5*.6, 10.5*.6)
(coinc, bins, patches) = plt.hist(dist_data, Nbins, ec='k', alpha=0.7)
plt.xlabel('time difference (ps)')

bin_width = bins[1] - bins[0]

print('Bin width (ps): ', bin_width)
# print('acc rate per bin (sec): ', a_rate*b_rate*bin_width*1e-12)
# print('# acc coincidences in bin: ', a_rate*b_rate*bin_width*1e-12*(a_range[-1][-1] - a_range[0][0])*1e-12)

# Gaussian function for fitting 
def Gaussian(t, amp, mean, std, floor):
    return amp*np.exp(-1/2*((t-mean)/std)**2) + floor


my = np.max(coinc)
param_init = [my, -1000, 500, 0]
# print(param_init)
    
popt, pcov = sp.optimize.curve_fit(Gaussian, bins[1:], coinc, p0=param_init, maxfev=2000)

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
 


# Second delay

dists_total = np.array([])
a_range = []
b_range = []
t1_len = []
b_len = []


if True:    
    
    for fi in files:
        while True:
            print('File: ', fi)
            fr = tt.FileReader(fi)
            
            # delay_ch2 = int(fi.split('_')[5]) 
            # print('Delay: ', delay_ch2)
                
            stream = fr.getData(streamSize)        
            
            ts = stream.getTimestamps()
            ch = stream.getChannels()
            
            
            if not len(ts):
                break
            
            t1 = np.array(ts[ch == 1])
            # tb = np.array(ts[ch == 2])
            tb = np.array(ts[ch == 3])
            
            
            print('t1 size: ', len(t1))
            t1_len += [len(t1)]
            min_t1 = min(t1) 
            max_t1 = max(t1)
            a_range += [[min_t1, max_t1]]
            print('time interval : ', 1e-12*(max_t1 - min_t1))
            
            
            print('tb size: ', len(tb))
            b_len += [len(tb)]
            min_tb = min(tb) 
            max_tb = max(tb)
            b_range += [[min_tb, max_tb]]
            print('time interval : ', 1e-12*(max_tb - min_tb))
            
            
            # dists = da.get_distances_B_to_A_simple(ta, tb, window) # find the distances and the indexes of the timestamps involved in Alice data
            dists, _ = da.get_distances_B_to_A_noq(t1, tb, n, offset_rough, limit=window) # find the distances and the indexes of the timestamps involved in Alice data
            
            
            # print('Dists : ', dists)
            
            dists_total = np.append(dists_total, dists)
            break



    print('Data points: ', len(dists_total))
    print('Total time (sec):')
    time = (a_range[-1][-1] - a_range[0][0])*1e-12
    print(a_range[0][0]*1e-12, a_range[-1][-1]*1e-12, time)
    print(b_range[0][0]*1e-12, b_range[-1][-1]*1e-12, (b_range[-1][-1] - b_range[0][0])*1e-12)
    
    
    print('a data points: ', sum(t1_len))
    t1_rate = sum(t1_len)/(time)
    print('a rate (sec): ', t1_rate)
    
    
    print('b data points: ', sum(b_len))
    b_rate = sum(b_len)/(time)
    b_rate2 = 1e12*sum(b_len)/(b_range[-1][-1] - b_range[0][0])
    
    print('b rate (sec): ', b_rate)
    print('b rate 2 (sec): ', b_rate2)        
    
    np.save('dists2.npy', dists_total)



print('............. Distances 1-3')
dist_data = np.load('dists2.npy')

     
Nbins = 200
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
param_init = [my, -1000, 500, 0]
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