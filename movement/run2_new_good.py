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



plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

c = sp.constants.c

tags = True
setup = mod.ModelMoving(tags=tags)

# NOTE: Time units are in microsecod
background_rate_sec = 454590 # background counts per second 136377 454590
qber = 0.0 # device qber not considering backgroud radiation
loss_channel = 13 # loss of the channel - only afects Bob
loss_detector = 6 # loss of the detection system - affects Alice and Bob
signal_rate_sec = 500000 # average entangled pair generation - 20 MHz
signal_rate = signal_rate_sec*1e-12
jitter = 305
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

qbers = []

err_count = 0
plot = 1

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
    
    
    dists, idx_b, qa_red = da.get_distances_moving2(ta, tb, qa, qb, n, offset_rough, velocity_rough, limit) # find the distances and the indexes of the timestamps involved in Alice data
    
    
    
    
    minb = np.min(tb)
    maxb = np.max(tb)
    
    # res_t = jitter/4
    res_t = 1
    # imgy = np.arange(minb, maxb, step=res_t)
    # imgx = np.arange(offset_rough- limit, offset_rough+limit, step=res_t)
    
    print('Res img : ', res_t*1e-12*c)
    print('Dists : ', len(dists))
    
    y = []
    x = []
    
    qy = []
    qx = []
    
    for i in tqdm(range(len(dists))):
        points = dists[i]
        
        # print(np.abs(di - offset_rough))
        # print(points)
        point_y = tb[idx_b[i]]
        
        point_qy = qb[idx_b[i]]
        qa_points = qa_red[i]
        # y_ind = np.argmin(np.abs(point_y - imgy))
        y_ind = int(point_y/res_t)
            
        for pi in range(len(points)):
            # x_ind = np.argmin(np.abs(pi - imgx))
            x_ind = int(points[pi]/res_t)
            
            y += [y_ind]
            x += [x_ind]
            
            qy += [point_qy]
            qx += [qa_points[pi]]
    
    x = np.array(x)
    y = np.array(y)
    qx = np.array(qx)
    qy = np.array(qy)
    
    # np.save('x.npy', x)
    # np.save('y.npy', y)
    
    
    # x = np.load('x.npy')
    # y = np.load('y.npy')
    
    # if plot:
        # plt.figure()
        # plt.title('Data points')
        # plt.scatter(x*1e-12, y*1e-12, marker='.')
        # plt.xlabel('Time difference (s)', fontsize=14)
        # plt.ylabel('Time Bob (s)', fontsize=14)
        
        # yref = np.linspace(np.min(y), np.max(y))
        # plt.plot(int(offset_rough/res_t)*np.ones(50), yref, 'g--')
    
    print('Img points : ', len(x))
    
    
    
    

    
    
    
    # Run the Hough Transform
    rho_res = 300000
    accumulator = hough_transform0(x, y, rho_res)
    
    
    peak = np.where(accumulator == np.max(accumulator))[0]
    
    if len(peak) > 1:
        # err_count += 1
        print("__________________ ERROR")
        continue
    
    if plot:
        fig = plt.figure()
        plt.title('Data points 2')
        fig.set_size_inches(18.5*.3, 18.5*.335)
        plt.scatter(x, y, marker='.')
        plt.xlabel('Time difference (ps)', fontsize=14)
        plt.ylabel('Time Bob (ps)', fontsize=14)
        plt.ylim([np.min(y), np.max(y)])
        plt.xlim([np.min(x), np.max(x)])
    
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
    
    qx = qx[x_ind]
    qy = qy[x_ind]
    
    
    if plot:
        plt.figure()
        plt.title('Reduced data points')
        plt.scatter(x_red, y_red, marker='.')
        plt.xlabel('Time difference (time jitter units)', fontsize=14)
        plt.ylabel('Time Bob (time jitter units)', fontsize=14)
    
    
    
    
    
    y = x_red
    X = y_red.reshape(-1, 1)
    
    # y = y_red
    # X = x_red.reshape(-1, 1)

    # Step 2: Fit the model using RANSAC
    ransac = RANSACRegressor(LinearRegression(), residual_threshold=43699/16)
    # ransac = RANSACRegressor(LinearRegression())

    ransac.fit(X, y)
    inlier_mask = ransac.inlier_mask_
    outlier_mask = ~inlier_mask
    
    # Step 3: Plot the data
    # line_X = np.arange(X.min(), X.max())[:, np.newaxis]
    line_X = np.linspace(np.min(X), np.max(X))
    line_X = line_X.reshape(-1, 1)
    line_y_ransac = ransac.predict(line_X)
    
    
    
    # Computing QBER here
    if 1:
        qx = qx[inlier_mask]
        qy = qy[inlier_mask]
       
        qber = sum(qx != qy)/len(qx)
    
    if plot:
        fig = plt.figure()
        fig.set_size_inches(18.5*.3, 10.5*.3)
        plt.scatter(y[inlier_mask].reshape(-1,1), X[inlier_mask].reshape(1,-1)[0], color="yellowgreen", marker="o", label="Inliers")
        plt.scatter(y[outlier_mask].reshape(-1,1), X[outlier_mask].reshape(1,-1)[0], color="red", marker="x", label="Outliers")
        plt.plot(line_y_ransac, line_X, color="cornflowerblue", linewidth=2, label="fit")
        plt.legend(loc="upper left", fontsize=14)
        plt.xlabel('Time difference (ps)', fontsize=14)
        plt.ylabel('Time Bob (ps)', fontsize=14)
        # plt.title("RANSAC Line Fitting")
        plt.ylim([0, 1e11])
        plt.xlim([3.325e7, 3.350e7])
        plt.tight_layout()
        plt.show()
    
    
    
    
    
    
    b = ransac.estimator_.intercept_
    m = ransac.estimator_.coef_
    
    # m  = m*1e-3
    print('m = ', m)
    print('b = ', b*res_t)
    
    # print('x intercept = ', -res_t*b)
    # print('velocity m = ', m)
    
    print('Delta diff: ', np.abs(time_shift-res_t*b))
    print('Delta vel: ', np.abs(300/c - m))
    
    print('--- Ref vals ')
    print('Time shift: ', time_shift)
    print('Velocity: ', 300/c)
    

    dds += [np.abs(time_shift-res_t*b)]
    dvs += [np.abs(300/c - m)]


    print("QBER = ", qber)
    qbers += [qber]

dds = np.array(dds)
dvs = np.array(dvs)
qbers = np.array(qbers)

errors = dds > 1000
dds = dds[np.bitwise_not(errors)]
dvs = dvs[np.bitwise_not(errors)]
err_count += sum(errors)

qbers = qbers[~errors]

print('------------- o -------------')
print('total loss: ', setup.get_loss())
print('Error diff: ', np.average(dds), np.std(dds))
print('Error vels: ', np.average(dvs), np.std(dvs))
print('Err count: ', err_count)

print('QBER: ', np.average(qbers), np.std(qbers))

sss= '''
total loss:  15
Error diff:  40.820970111228526 28.717149530143885
Error vels:  4.2112580162928574e-10 3.631042014850005e-10
Err count:  0

total loss:  20
Error diff:  48.584917260706426 44.36871016956952
Error vels:  7.673983046326043e-10 5.723820864914807e-10
Err count:  0

total loss:  25
Error diff:  80.36225747670979 62.2007251709788
Error vels:  1.2498071652930094e-09 1.0474605586832894e-09
Err count:  0

total loss:  30
Error diff:  181.9473782433297 135.27061144275962
Error vels:  3.2689495175947456e-09 2.728498965785199e-09
Err count:  7
'''

'''
Field trial calculation
background_rate_sec = 136377
integration = 0.1

total loss:  15
Error diff:  31.92344809278846 19.828175218752314
Error vels:  3.1130199810633156e-10 2.0841026068651124e-10
Err count:  0
QBER:  0.00415939968505794 0.0014940990559353378

total loss:  20
Error diff:  51.355698643177746 39.67255147998592
Error vels:  7.03665415335368e-10 5.54448175227094e-10
Err count:  0
QBER:  0.011313337360835443 0.004389284591218541

total loss:  25
Error diff:  90.96978184863924 81.31993724344188
Error vels:  1.4622598392817626e-09 1.2185100628862228e-09
Err count:  0
QBER:  0.032676481171940434 0.014602986147228877

total loss:  30
Error diff:  64.4291210360825 0.0
Error vels:  1.2546212681612296e-09 0.0
Err count:  0
QBER:  0.08695652173913043 0.054602986147228877



background_rate_sec = 454590
integration = 0.1
total loss:  15
Error diff:  36.574893046095966 25.776354806653355
Error vels:  3.654807880195206e-10 2.802758251617822e-10
Err count:  0
QBER:  0.011154727882613043 0.002721929792795183

total loss:  20
Error diff:  75.45171325311065 118.51948513002833
Error vels:  1.5767393672420756e-09 2.64222141900338e-09
Err count:  0
QBER:  0.031646912622270665 0.0076079605794522635

total loss:  25
Error diff:  220.88224605224647 232.77947074594593
Error vels:  9.697633729683261e-09 3.0256607898839145e-08
Err count:  41
QBER:  0.108823529411764706 0.06785534721

total loss:  30
Error diff:  284.6084078512258 224.5720284898256
Error vels:  6.778590641175835e-09 7.304467008601469e-09
Err count:  53
QBER:  0.3478630821891544 0.1084694974500704

'''

ll = [15, 20, 25, 30]




qber3 = [0.004, 0.011, 0.032, 0.086]
std3 = [0.0014, 0.0043, 0.0146, 0.0246]

qber10 = [0.011, 0.031, 0.108, 0.202]
std10 = [0.0027, 0.0076, 0.0278, 0.0346]



fig, ax1 = plt.subplots(figsize=(8, 6))

fig.set_size_inches(18.5*.3, 10.5*.3)

ax1.set_xlabel('Total loss (dB)', fontsize=14)
ax1.set_ylabel('QBER', fontsize=14)
ax1.set_ylim(0, 0.25)

ax1.errorbar(ll, qber3, std3, label='3nm filter')
ax1.errorbar(ll, qber10, std10, label='10nm filter')
ax1.legend()

# COLOR_TIME = "#325fdb"
# COLOR_VELOCITY = "#c4661d"

# # Instantiate a second axes that shares the same x-axis
# ax2 = ax1.twinx()  

# ax1.plot(ll, ed, '-o', color=COLOR_TIME)
# ax2.plot(ll, ev, '-*', color=COLOR_VELOCITY)

# ax1.set_xlabel('Total loss (dB)', fontsize=14)
# ax1.set_ylabel('Timing precision (ps)', color=COLOR_TIME, fontsize=14)
# ax2.set_ylabel('Drift  precision (ns/s)', color=COLOR_VELOCITY, fontsize=14)
# ax1.tick_params(axis="y", labelcolor=COLOR_TIME)
# ax2.tick_params(axis="y", labelcolor=COLOR_VELOCITY)

# # ax2.set_yscale('log')
# # ax1.set_yscale('log')
# ax2.set_ylim(.3, 3)
# ax1.set_ylim(20, 200)


# plt.tight_layout()
# plt.show()




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

