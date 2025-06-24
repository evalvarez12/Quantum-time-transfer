# -*- coding: utf-8 -*-

import sys
sys.path.append('./src')
import numpy as np
import TimeTagger as tt
import distance_analysis as da
import cross_corr as cc
import matplotlib.pyplot as plt
import fine_analysis as fa
import scipy as sp
import os
from tqdm import tqdm

import glob



def rotate_array(a, x):
    return np.concatenate((a[-x:], a[:-x]))


def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

plt.close('all')

# # All files and directories ending with .txt and that don't begin with a dot:
# print(glob.glob("/home/adam/*.txt")) 
# # All files and directories ending with .txt with depth of 2 folders, ignoring names beginning with a dot:
# print(glob.glob("/home/adam/*/*.txt")) 


os.chdir(r'C:\Users\vil034\OneDrive - CSIRO\Documents\laser data\aft_Ch3_delay_adjust\500_ps')


print(os.getcwd())


all_files = glob.glob('*.ttbin')
files1 = glob.glob('*.1.ttbin')

all_files = [i for i in all_files if i not in files1]

a_range = []
b_range = []
a_len = []
b_len = []

dists_total = np.array([])

# windows = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
# coincidences = np.array([0]*10) 
# tagged = np.array([0]*10)


delay_arm1 = -1307
delay_arm2 = -1264


# Analog array parameters
dt = 500 # units picosec
N = int(.250/(dt*1e-12))

# taus = np.arange(-2000, 2000, 800)
# taus = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]
taus = [0]
coincidence_window = dt
g2s = []


files = [all_files[202]]
# files = files[400:450]
print(files)

Nfiles = len(files)

nc12_sum = 0
nc1_sum = 0 
nc2_sum = 0
     

# for fi in tqdm(files):
for tau in taus:

    print('tau = ', tau)
    for fi in files:
        print()
    # fi = files[-20]
        # print('File: ', fi)
        fr = tt.FileReader(fi)
        
        # delay_ch2 = int(fi.split('_')[5]) 
        # print('Delay: ', delay_ch2)
            
        stream = fr.getData(1000000)        
        
        ts = stream.getTimestamps()
        ch = stream.getChannels()
        
        t1 = np.array(ts[ch == 1])
        t2 = np.array(ts[ch == 2])
        t3 = np.array(ts[ch == 3])
        
        
        l = len(t1)
        print('t1 size: ', l)
        print('t2 size: ', len(t2))
        print('t3 size: ', len(t3))
        
        t2 = t2 - delay_arm1 + tau
        t3 = t3 - delay_arm2
        
        
        # nt1 = cc.to_analog(t1, N, dt)
        # nt2 = cc.to_analog(t2, N, dt)
        # nt3 = cc.to_analog(t3, N, dt)
        
        
        # nc1 = nt1*nt2
        # nc2 = nt1*nt3
        
        
        # # nc1 = rotate_array(nc1, 51)        
        
        
        # g2 = np.average(nc1*nc2)/(np.average(nc1)*np.average(nc2))
        
        # print()
        # print('g2 : ', g2)
        
        # print('nc1 : ', np.sum(nc1))
        # print('nc2 : ', np.sum(nc2))
        # print('nc1 nc2 : ', np.sum(nc1*nc2))
        
        # nc12_sum += np.sum(nc1*nc2)
        # nc1_sum += np.sum(nc1)
        # nc2_sum += np.sum(nc2)
    
        
        # print('Alternative g2 calculation')
        n1 = np.zeros(l, dtype=bool)
        n2 = np.zeros(l, dtype=bool)
        
        for i in range(l):
            ti = t1[i]
            i2 = find_nearest_idx(t2, ti)
            i3 = find_nearest_idx(t3, ti)
            
            n1[i] = np.abs(ti - t2[i2]) < coincidence_window
            n2[i] = np.abs(ti - t3[i3]) < coincidence_window
            
            
            
        print()    
        print('n1 : ', sum(n1))
        print('n2 : ', sum(n2))
        print('n1 n2 : ', sum(n1*n2))
        
        print(np.average(n1*n2)/(np.average(n1)* np.average(n2)))
        g2i = (np.average(n1*n2))/(np.average(n1)*np.average(n2))
        # g2 += g2i
        print('g2 = ', g2i)
            
        
        g2s += [g2i]
            

    
    
# a = nc12_sum/N/Nfiles
# b = nc1_sum/N/Nfiles
# c = nc2_sum/N/Nfiles
# print('Averaged g2')
# print('g2 : ', a/b/c)
    
    # g2s += [g2/len(files)]
    
    # print(g2s)
        
# print('g2 = ', (np.average(n1*n2))/(np.average(n1)*np.average(n2))

# taus = range(0,1800, 100)
# g2s = np.load('g2_data.npy')

# plt.plot(taus, g2s, 'o-')
# plt.xlabel(r'$\tau$ (picosec)', fontsize=16)
# plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=16)
# plt.yticks(fontsize=16)
# plt.xticks(fontsize=16)
# plt.tight_layout()


