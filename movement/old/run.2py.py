# -*- coding: utf-8 -*-


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
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import LinearRegression

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

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
velocity = 300 # m/s
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


dds = []
dvs = []

err_count = 0
plot = True

for ii in range(1):

    
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
    
    res_t = jitter/4
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
    
    if plot:
        plt.figure()
        plt.title('Data points')
        plt.scatter(x, y, marker='.')
        plt.xlabel('Time difference (time jitter units)', fontsize=14)
        plt.ylabel('Time Bob (time jitter units)', fontsize=14)
        
        yref = np.linspace(np.min(y), np.max(y))
        plt.plot(int(offset_rough/res_t)*np.ones(50), yref, 'g--')
    
    print('Img points : ', len(x))
    
    
    
    
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
    
    
    
    # Run the Hough Transform
    rho_res = 2000
    accumulator = hough_transform0(x, y, rho_res)
    
    
    peak = np.where(accumulator == np.max(accumulator))[0]
    
    if len(peak) > 1:
        err_count += 1
        continue
    
    if plot:
        fig = plt.figure()
        plt.title('Data points 2')
        fig.set_size_inches(18.5*.3, 18.5*.335)
        plt.scatter(x, y, marker='.')
        plt.xlabel('Time difference (time jitter units)', fontsize=14)
        plt.ylabel('Time Bob (time jitter units)', fontsize=14)
    
    yref = np.linspace(np.min(y), np.max(y))
    # plt.plot(int(offset_rough/res_t)*np.ones(50), yref, 'g--')
    
    xa = np.min(x)+(peak)*rho_res
    xb = np.min(x)+(peak+1)*rho_res
    if plot:
        plt.plot(xa*np.ones(50), yref, 'r--')
        plt.plot(xb*np.ones(50), yref, 'r--')
        plt.tight_layout()
    
    
    x_ind = np.where(np.logical_and(x > xa, x < xb))
    # x_red = x[x_ind]*1e-3
    # y_red = y[x_ind]*1e-7
    x_red = x[x_ind]
    y_red = y[x_ind]
    
    
    if plot:
        plt.figure()
        plt.title('Reduced data points')
        plt.scatter(x_red, y_red, marker='.')
        plt.xlabel('Time difference (time jitter units)', fontsize=14)
        plt.ylabel('Time Bob (time jitter units)', fontsize=14)
    
    
    
    
    
    y = y_red
    X = x_red.reshape(-1, 1)
    
    # X = y_red.reshape(-1, 1)
    # y = x_red
    
    # Step 2: Fit the model using RANSAC
    ransac = RANSACRegressor(LinearRegression(), residual_threshold=290556242.5/8)
    ransac.fit(X, y)
    inlier_mask = ransac.inlier_mask_
    outlier_mask = ~inlier_mask
    
    # Step 3: Plot the data
    # line_X = np.arange(X.min(), X.max())[:, np.newaxis]
    line_X = np.linspace(np.min(X), np.max(X))
    line_X = line_X.reshape(-1, 1)
    line_y_ransac = ransac.predict(line_X)
    
    
    if plot:
        fig = plt.figure()
        fig.set_size_inches(18.5*.3, 10.5*.3)
        plt.scatter(X[inlier_mask], y[inlier_mask], color="yellowgreen", marker="o", label="Inliers")
        plt.scatter(X[outlier_mask], y[outlier_mask], color="red", marker="x", label="Outliers")
        plt.plot(line_X, line_y_ransac, color="cornflowerblue", linewidth=2, label="fit")
        plt.legend(loc="upper left")
        plt.xlabel('Time difference (time jitter units)', fontsize=14)
        plt.ylabel('Time Bob (time jitter units)', fontsize=14)
        # plt.title("RANSAC Line Fitting")
        plt.ylim([np.min(y), np.max(y)])
        plt.tight_layout()
        plt.show()
    
    
    
    
    
    
    b = ransac.estimator_.intercept_
    m = ransac.estimator_.coef_
    
    # m  = m*1e-3
    print('m = ', m)
    print('b = ', b*res_t)
    
    print('x intercept = ', -res_t*b/(m))
    print('velocity m = ', c/m)
    
    print('Delta diff: ', np.abs(time_shift+res_t*b/(m)))
    print('Delta vel: ', np.abs(300 - c/m))
    
    print('--- Ref vals ')
    print('Time shift: ', time_shift)
    print('Velocity: ', 300)
    

    dds += [np.abs(time_shift+res_t*b/(m))]
    dvs += [np.abs(300 - c/m)]


print('total loss: ', setup.get_loss())
print('Error diff: ', np.average(dds), np.std(dds))
print('Error vels: ', np.average(dvs), np.std(dvs))
print('Err count: ', err_count)

sss= '''
total loss:  15
Error diff:  101.68230206441135 28.396622426748348
Error vels:  0.15888312564611395 0.11841693254406396
Err count:  0

total loss:  20
Error diff:  107.33200778983533 59.09182694690947
Error vels:  0.2918876511773601 0.20476116943432737
Err count:  0

total loss:  25
Error diff:  139.22031704321503 90.76985826399805
Error vels:  0.5617551453127027 0.46085025043166755
Err count:  0

total loss:  30
Error diff:  228.57382622763816 193.2234814101607
Error vels:  1.1561640104359012 0.9338750023485088
Err count:  2
'''


ll = [15, 20, 25, 30]
ed = [28.3966, 59.091, 90.769, 193.223]
ev = np.array([ 0.118, 0.204, 0.460, 0.933])
ev = ev/c * 1e9

fig, ax1 = plt.subplots(figsize=(8, 6))

fig.set_size_inches(18.5*.3, 10.5*.3)


COLOR_TIME = "#325fdb"
COLOR_VELOCITY = "#c4661d"

# Instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()  

ax1.plot(ll, ed, '-o', color=COLOR_TIME)
ax2.plot(ll, ev, '-*', color=COLOR_VELOCITY)

ax1.set_xlabel('Total loss (dB)', fontsize=14)
ax1.set_ylabel('Time precision (ps)', color=COLOR_TIME, fontsize=14)
ax2.set_ylabel('Drift  precision (ns/s)', color=COLOR_VELOCITY, fontsize=14)
ax1.tick_params(axis="y", labelcolor=COLOR_TIME)
ax2.tick_params(axis="y", labelcolor=COLOR_VELOCITY)

ax2.set_ylim(0, 3.5)
ax1.set_ylim(0, 200)

plt.tight_layout()
plt.show()




# COLOR_TEMPERATURE = "#69b3a2"
# COLOR_PRICE = "#3399e6"

# fig, ax1 = plt.subplots(figsize=(8, 8))
# ax2 = ax1.twinx()

# ax1.plot(date, temperature, color=COLOR_TEMPERATURE, lw=3)
# ax2.plot(date, price, color=COLOR_PRICE, lw=4)

# ax1.set_xlabel("Date")
# ax1.set_ylabel("Temperature (Celsius °)", color=COLOR_TEMPERATURE, fontsize=14)
# ax1.tick_params(axis="y", labelcolor=COLOR_TEMPERATURE)

# ax2.set_ylabel("Price ($)", color=COLOR_PRICE, fontsize=14)
# ax2.tick_params(axis="y", labelcolor=COLOR_PRICE)

# fig.suptitle("Temperature down, price up", fontsize=20)
# fig.autofmt_xdate()

# plt.show()

