# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:20:06 2023

"""
import sys
sys.path.append('.././src')

import model_moving as mod
import numpy as np
import distance_analysis as da
import matplotlib.pyplot as plt
import fine_analysis as fa
import scipy as sp
from tqdm import tqdm


plt.close('all')

c = sp.constants.c

tags = True
setup = mod.ModelMoving(tags=tags)

# NOTE: Time units are in microsecod
background_rate_sec = 37000 # background counts per second 
qber = 0.05 # device qber not considering backgroud radiation
loss_channel = 10 # loss of the channel - only afects Bob
loss_detector = 5 # loss of the detection system - affects Alice and Bob
signal_rate_sec = 500000 # average entangled pair generation - 20 MHz
signal_rate = signal_rate_sec*1e-12
jitter = 364
velocity = 100 # m/s
# jitter = 10e-6 

background_rate = background_rate_sec*1e-12
setup.set_background_rate(background_rate)
setup.set_qber(qber)
setup.set_loss(loss_channel, loss_detector)
setup.set_jitter(jitter)
setup.set_signal_rate(signal_rate)
setup.set_velocity(velocity)

setup.set_detector_dead_time(0)

time_shift = 33356000 # clock offset to be found in microseconds
start_time = 0 # zero position of the Alice clock
total_protocol_time_sec = .1 # total duration of the protocol in seconds
total_protocol_time = total_protocol_time_sec * 1e12 # protocol time picosec

# STEP 1: Generate data using a monte carlo simulation - for simplicity only 1 basis is used on the polarization
# ta = timing data Alice
# tb = timing data Bob - includes background radiation
# qa = quantum states prepared by Alice
# qb = quantum states measured by Bob - includes background radiation
ta, tb, qa, qb = setup.generate_time_data(time_shift, start_time, total_protocol_time_sec, verbose=True)


# STEP 2: Find the peak using the rough estimate from the GPS
offset_rough = 33011000 # rough offset given by the GPS - approx 90 meters off 
velocity_rough = 100 # units of m/s 
n = 4 # this parameter sets the number of timestamps to each side of the rough_offset to be used to find the peak
limit = 1*1e6

print('Approximated positon :', 1e-12*offset_rough*c,  '+-', 1e-12*limit*c)


dists, idx_b = da.get_distances_moving2(ta, tb, qa, qb, n, offset_rough, velocity_rough, limit) # find the distances and the indexes of the timestamps involved in Alice data




minb = np.min(tb)
maxb = np.max(tb)

res_t =2*jitter
# imgy = np.arange(minb, maxb, step=res_t)
# imgx = np.arange(offset_rough- limit, offset_rough+limit, step=res_t)

print('Res img : ', res_t*1e-12*c)
print('Dists : ', len(dists))

y = []
x = []

for i in tqdm(range(len(dists))):
    di = dists[i]
    
    points = dists[i]
    # print(np.abs(di - offset_rough))
    # print(points)
    point_y = tb[idx_b[i]]
    # y_ind = np.argmin(np.abs(point_y - imgy))
    y_ind = int(point_y/res_t)
        
    for pi in points:
        # x_ind = np.argmin(np.abs(pi - imgx))
        x_ind = int(pi/res_t)
        
        y += [y_ind]
        x += [x_ind]
        


np.save('x.npy', x)
np.save('y.npy', y)


x = np.load('x.npy')
y = np.load('y.npy')


plt.figure()
plt.scatter(x, y, marker='.')
plt.xlabel('Time difference (time jitter units)', fontsize=14)
plt.ylabel('Time Bob (time jitter units)', fontsize=14)

yref = np.linspace(np.min(y), np.max(y))
plt.plot(int(offset_rough/res_t)*np.ones(50), yref, 'g--')

print('Img points : ', len(x))



#########
#########
###### fix zero point in the data



def hough_transform0(x_points, y_points, rho_res):
    
    x_min = np.min(x_points)
    x_max = np.max(x_points)
    
    N_rho = int((x_max - x_min)/rho_res) + 1
    accumulator = np.zeros(N_rho)
    
    print('N_rho : ', N_rho,)
    accumulator = np.zeros(N_rho, dtype=int)

    # Hough Transform: vote in the accumulator for each point
    for x, y in zip(x_points, y_points):
        rho_index = int((x - x_min)/rho_res)
        accumulator[rho_index] += 1
    return accumulator


def hough_transform(x_points, y_points, theta_res=1, rho_res=1):
    # Convert angles to radians
    theta_range = np.deg2rad(np.arange(0, 180, theta_res))
    N_theta = len(theta_range)
    # Create an empty accumulator for rho and theta values
    # print(np.max(x_points), np.max(y_points))
    
    
    max_rho = int(np.sqrt(np.int64(np.max(x_points))**2 + np.int64(np.max(y_points))**2))
    min_rho = int(np.sqrt(np.int64(np.min(x_points))**2 + np.int64(np.min(y_points))**2))
    # rhos = np.arange(-max_rho, max_rho, rho_res)
    N_rhos = int((max_rho - min_rho)/rho_res) + 1
    
    print('N_rhos : ', N_rhos, '  N theta : ', N_theta)
    accumulator = np.zeros((N_rhos, N_theta), dtype=int)

    # Hough Transform: vote in the accumulator for each point
    for x, y in zip(x_points, y_points):
        for theta_index, theta in enumerate(theta_range):
            rho = y * np.cos(theta) - x * np.sin(theta)
            # rho_index = np.argmin(np.abs(rhos - rho))  # Find the closest rho bin
            rho_index = int((rho-min_rho)/rho_res)
            accumulator[rho_index, theta_index] += 1   # Increment the accumulator

    return accumulator, theta_range, max_rho





def line_space(xs, ys ):
    n = len(xs)
    nn = int(((n-1)**2 + n - 1)/2)
    
    print('n: ', n)
    print('nn : ', nn)
    
    ys = ys*1e-7

    ms = np.zeros(nn)
    bs = np.zeros(nn)
    
    k = 0
    for i in range(n):
        x1 = xs[i]
        y1 = ys[i]
        for j in range(i+1, n):
            x2 = xs[j]
            y2 = ys[j]
            
            yd = y2 - y1
            xd = x2 - x1
            if yd == 0 or xd ==0:
                ms[k] = -1
                bs[k] = -1
            else:
                ms[k] = (x2 - x1)/yd
                bs[k] = x1 - ms[k]*y1
            k += 1
                
    ms = ms[ms != -1]
    bs = bs[bs != -1]
    print('k: ', k)
    
    inds = np.where(np.logical_and(ms < 100, ms > -100))[0]
    

    return ms[inds], bs[inds]
    # return ms, bs


def find_hough_peak(accumulator):
    x, y = np.where(accumulator == np.max(accumulator))
    return (x[0], y[0])



# Run the Hough Transform
rho_res = 200
accumulator = hough_transform0(x, y, rho_res)


peak = np.where(accumulator == np.max(accumulator))[0]

# plt.figure()
# plt.scatter(x, y, marker='.')
# plt.xlabel('Time difference (picosec)', fontsize=14)
# plt.ylabel('Time Alice', fontsize=14)

yref = np.linspace(np.min(y), np.max(y))
# plt.plot(int(offset_rough/res_t)*np.ones(50), yref, 'g--')

xa = np.min(x)+(peak)*rho_res
xb = np.min(x)+(peak+1)*rho_res
plt.plot(xa*np.ones(50), yref, 'r--')
plt.plot(xb*np.ones(50), yref, 'r--')


x_ind = np.where(np.logical_and(x > xa, x < xb))
x_red = x[x_ind]
y_red = y[x_ind]

plt.figure()
plt.scatter(x_red, y_red*1e-7, marker='.')
plt.xlabel('Time difference (time jitter units 10-7)', fontsize=14)
plt.ylabel('Time Alice (time jitter units)', fontsize=14)




ms, bs = line_space(x_red, y_red)



# plt.figure()
# plt.scatter(ms, bs, marker='.')

# plt.figure()
# plt.hist(ms, 1000)
# plt.figure()
# plt.hist(bs, 1000)


ms_hist, ms_bins = np.histogram(ms, 1000)
bs_hist, bs_bins = np.histogram(bs, 1000)

# plt.figure()
# plt.plot(ms_bins[:-1], ms_hist, '-')
# plt.figure()
# plt.plot(bs_bins[:-1], bs_hist, '-')


mmax = np.max(ms_hist)
mind = np.where(ms_hist == mmax)[0][0]
ma = ms_bins[mind]
mb = ms_bins[mind+1]

bmax = np.max(bs_hist)
bind = np.where(bs_hist == bmax)[0][0]
ba = bs_bins[bind]
bb = bs_bins[bind+1]


print('m in ' , ma*1e5, mb*1e5)
print('b in ' , ba*res_t, bb*res_t)

print('m precision: ', (mb - ma)*1e5)
print('b precision: ', (bb - ba)*res_t)  


m = np.average(ms)
b = np.average(bs)*res_t
print('Averages')
print('m = ', m*1e5)
print('b = ', b)


print('Time shift: ', time_shift)
print('Velocity: ', 100/c*1e12)




y_line = np.linspace(np.min(y_red*1e-7), np.max(y_red*1e-7))
plt.plot(b/res_t+m*y_line, y_line, 'r--')


# rho_res = 5
# accumulator, theta_range, max_rho = hough_transform(x, y, 5,5)


# # Display the Hough accumulator
# plt.figure()
# plt.imshow(accumulator, extent=[np.rad2deg(theta_range.min()), np.rad2deg(theta_range.max()), -max_rho, max_rho], aspect='auto')
# plt.xlabel('Theta (degrees)')
# plt.ylabel('Rho')
# plt.title('Hough Space')
# plt.colorbar(label='Votes')



# peak = find_hough_peak(accumulator)

# print('Peaks: ', peak)

# rhos = np.arange(-max_rho, max_rho, rho_res)
# # Visualize the detected peaks in Hough space
# rho_idx, theta_idx = peak
# plt.scatter(np.rad2deg(theta_range[theta_idx]), -rhos[rho_idx], color='red')

# plt.show()











# plt.figure()
# plt.imshow(accumulator, extent=[np.rad2deg(theta_range.min()), np.rad2deg(theta_range.max()), rhos.min(), rhos.max()], aspect='auto')
# plt.xlabel('Theta (degrees)')
# plt.ylabel('Rho')
# plt.title('Hough Space with Peaks')
# plt.colorbar(label='Votes')
# plt.show()


# plt.figure()
# plt.scatter(x, y, marker='.')
# plt.xlabel('Time difference (picosec)', fontsize=14)
# plt.ylabel('Time Alice', fontsize=14)


# xline = np.linspace(min(x), max(x))
# yline = (-xline*np.cos(theta_range[theta_idx]) + rhos[rho_idx])/np.sin(theta_range[theta_idx])

# yline = np.linspace(min(y), max(y))
# xline = (rhos[rho_idx] - yline*np.sin(theta_range[theta_idx]))/np.cos(theta_range[theta_idx])

# plt.plot(xline, yline, 'r--')