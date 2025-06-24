
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
from scipy.stats import median_abs_deviation




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
background_rate_sec = 5330000/50 # background counts per second 
qber = 0.05 # device qber not considering backgroud radiation
loss_channel = 20 # loss of the channel - only afects Bob
loss_detector = 5 # loss of the detection system - affects Alice and Bob
signal_rate_sec = 10000000 # average entangled pair generation - 20 MHz
signal_rate = signal_rate_sec*1e-12
jitter = 302
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
plot = 1

for ii in range(0):
    print("-------------------------------->", ii)

    
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
        err_count += 1
        continue
    
    if plot:
        ########################## ATTENUATION
        mask = np.random.choice(len(x), size=int(0.01*len(x)), replace=False)
        xp = 1e-12*x[mask]
        yp = 1e-12*y[mask]
        ############################
        
        
        fig = plt.figure()
        plt.title('Data points 2')
        fig.set_size_inches(18.5*.3, 18.5*.335)
        plt.scatter(xp, yp, marker='.')
        plt.xlabel('Time difference (ps)', fontsize=14)
        plt.ylabel('Time Bob (ps)', fontsize=14)
        plt.ylim([np.min(yp), np.max(yp)])
        plt.xlim([np.min(xp), np.max(xp)])
    
    yref = np.linspace(np.min(y), np.max(y))
    # plt.plot(int(offset_rough/res_t)*np.ones(50), yref, 'g--')
    
    xa = np.min(x)+(peak)*rho_res
    xb = np.min(x)+(peak+1)*rho_res
    if plot:
        plt.plot(1e-12*xa*np.ones(50), 1e-12*yref, 'r--')
        plt.plot(1e-12*xb*np.ones(50), 1e-12*yref, 'r--')
        plt.tight_layout()
    
    
    x_ind = np.where(np.logical_and(x > xa, x < xb))
    # x_red = x[x_ind]*1e-3
    # y_red = y[x_ind]*1e-7
    x_red = x[x_ind]
    y_red = y[x_ind]
    
    qx = qx[x_ind]
    qy = qy[x_ind]
    
    
    if 0:
        plt.figure()
        plt.title('Reduced data points')
        plt.scatter(x_red, y_red, marker='.')
        plt.xlabel('Time difference (time jitter units)', fontsize=14)
        plt.ylabel('Time Bob (time jitter units)', fontsize=14)
    
    
    y = x_red
    X = y_red.reshape(-1, 1)
    
    # Step 2: Fit the model using RANSAC
    # res_threshold = median_abs_deviation(y)/(16)
    res_threshold = 1000
    print('Res thres = ', res_threshold)
    # res_threshold = 500
    ransac = RANSACRegressor(LinearRegression(), residual_threshold=res_threshold, max_trials=400)
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
       
        qber = sum(qx == qy)/len(qx)
    
    if plot:
        xi = X[inlier_mask]
        yi = y[inlier_mask]
        
        xo = X[outlier_mask]
        yo = y[outlier_mask]
        
        maski = np.random.choice(len(xi), size=int(0.01*len(xi)), replace=False)
        masko = np.random.choice(len(xo), size=int(0.01*len(xo)), replace=False)
        
        xi = xi[maski]
        yi = yi[maski]
        
        xo = xo[masko]
        yo = yo[masko]
        
        
        fig = plt.figure()
        fig.set_size_inches(18.5*.3, 10.5*.3)
        plt.scatter(yi.reshape(-1,1)*1e-6, xi.reshape(1,-1)[0]*1e-12, color="yellowgreen", marker="o", label="Inliers")
        plt.scatter(yo.reshape(-1,1)*1e-6, xo.reshape(1,-1)[0]*1e-12, color="red", marker="x", label="Outliers")
        plt.plot(line_y_ransac*1e-6, line_X*1e-12, color="cornflowerblue", linewidth=2, label="fit")
        plt.legend(loc="upper left", fontsize=14)
        plt.xlabel(r'Time difference ($\mu$s)', fontsize=14)
        plt.ylabel('Time Bob (s)', fontsize=14)
        # plt.title("RANSAC Line Fitting")
        plt.ylim([0, .1])
        plt.xlim([33.25, 33.50])
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


dds_all = np.array(dds)
dvs_all = np.array(dvs)


errors = np.array(dds) > 2000 #########################
dds = dds_all[np.bitwise_not(errors)]
dvs = dvs_all[np.bitwise_not(errors)]
err_count += sum(errors)

print('------------- o -------------')
print('total loss: ', setup.get_loss())
print('Error diff: ', np.average(dds), np.std(dds))
print('Error vels: ', np.average(dvs), np.std(dvs))
print('Err count: ', err_count)

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


NEW DATA
total loss:  15
Error diff:  71.59543872497976 69.43533832345268
Error vels:  1.0348693215643574e-09 1.1681709536410141e-09
Err count:  0

total loss:  20
Error diff:  85.41459809165448 81.03803403539081
Error vels:  1.4761841818392408e-09 1.6406780965425883e-09
Err count:  0

total loss:  25
Error diff:  86.9264834022522 71.93654850899173
Error vels:  1.3446890834059514e-09 1.9265391614368513e-09
Err count:  0

total loss:  30
Error diff:  214.15342382738987 221.332321124156
Error vels:  4.133052902666961e-09 5.597230814015387e-09
Err count:  10

total loss:  35
Error diff:  625.7618935500582 552.688390515286
Error vels:  2.7859963869797134e-08 3.8050363880291894e-08
Err count:  55


'''


ll = [15, 20, 25, 30]
ed = [69.43533832345268, 101.03803403539081, 151.93654850899173, 221.332321124156]
ev = np.array([ 1.1681709536410141e-09, 2.0406780965425883e-09, 3.2265391614368513e-09, 5.597230814015387e-09])
ev = ev*1e9


fig, ax1 = plt.subplots(figsize=(8, 6))

fig.set_size_inches(18.5*.3, 10.5*.3)


COLOR_TIME = "#FF2B06"
COLOR_VELOCITY = "#06DAFF"

# Instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()  

ax1.plot(ll, ed, '-o', color=COLOR_TIME)
ax2.plot(ll, ev, '-*', color=COLOR_VELOCITY)

ax1.set_xlabel('Total loss (dB)', fontsize=14)
ax1.set_ylabel('Timing precision (ps)', color=COLOR_TIME, fontsize=14)
ax2.set_ylabel('Drift  precision (ns/s)', color=COLOR_VELOCITY, fontsize=14)
ax1.tick_params(axis="y", labelcolor=COLOR_TIME)
ax2.tick_params(axis="y", labelcolor=COLOR_VELOCITY)

# ax2.set_yscale('log')
# ax1.set_yscale('log')
ax2.set_ylim(.5,6)
ax1.set_ylim(20, 250)


plt.tight_layout()
plt.show()