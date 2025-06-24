# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'

# dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]

dist = 10000



scattering_loss = 0.45 # db/km for 810

data_day = sp.io.loadmat(file_day %dist)
data_day = data_day['res']*10**(-scattering_loss*dist*1e-3/10)
data_day = np.transpose(data_day)[0]


plt.show()


ts = data_day
loss = -10*np.log10(data_day)


fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.hist(loss, 100)

plt.xlabel('Loss (dB)', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)



N = int(20e3)
ps = []

ts_sample = np.random.choice(ts, N)
    
    
    # print('Mean loss: ', np.mean(loss))
    # print('Var loss: ', np.std(loss)**2)
    
    
roll = np.random.rand(N)
clicks = roll < ts_sample
    
p = sum(clicks)/N


# print('Total clicks: ', sum(clicks))
# print('Probability: ', p)
    
# print('Distance from pop mean: ', p - mean)
# print('Distance variance expected: ', np.sqrt(p*(1-p)/N))
    

print('-------------- In dB') 
print('Loss mean: ', np.mean(loss))
print('Loss std: ', np.std(loss))
print('-----------')    
print('Loss mean: ', np.mean(ts))
print('Loss std: ', np.std(ts))
      
print(' Channel sample ----------------')
print('Loss sample mean: ', np.mean(ts_sample))
print('Loss sample std: ', np.std(ts_sample))

print('Estimation pres:', np.abs(np.mean(ts) - np.mean(ts_sample)))
print ('Expected pres: ', np.sqrt(np.std(ts_sample)/N))


print('Clicks ---------------------')

print('Total clicks: ', np.sum(clicks))

print('Mean clicks: ', np.mean(clicks))
print('Std clicks: ', np.std(clicks))

print('Std clicks theory: ', np.sqrt(p*(1-p)))
print('Estimation pres: ', np.abs(p - np.mean(ts)))
print('Expected pres: ', np.sqrt(p*(1-p)/N))