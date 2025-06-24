# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-

import optimizer
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import interpolator
import scipy.constants as spc
import argparse


import warnings
warnings.filterwarnings('ignore')

plt.close('all')

############################################
# CALCULATION OF THE BACKGROUND COUNTS

# cn2_day = 2.47e-13
cn2_day = 1.52e-14
cn2_night = 1.25e-15
wl = 810e-9

det_loss = 5
fov_angle_AO = 20e-6
# fov_angle_nAO = 100e-6 # 182 - old value
fov = 4*np.pi*np.sin(fov_angle_AO/2)**2

spectral_width = 2.3
# spectral_width = 0.1
sky_radiance = 0.004
a = 0.075 # aperture radious 


background_day = spectral_width*sky_radiance*fov*np.pi*(a**2)*wl/(spc.c* spc.h)*10**(-det_loss/10)
background_night = spectral_width*(sky_radiance/50)*fov*np.pi*(a**2)*wl/(spc.c* spc.h)*10**(-det_loss/10)

# print('Background day: ', background_day)
# print('Background night: ', background_night)
###################################################################

calc = 0


if calc:
    parser = argparse.ArgumentParser(
            description="Parameters to compute the protocol"
        )
    
    parser.add_argument("--mode", required=True, type=str)
    parser.add_argument("--mode2", required=True, type=str)
    parser.add_argument("--ra", required=True, type=int)
    args = parser.parse_args()
    
    ras = args.ra
    mode = args.mode
    mode2 = args.mode2

    print('Running --', mode, mode2, ras)


    ra = ras*10e6 # pulse rate
    jitter = 345e-12 # total time jitter
    e_d = 0.05 # detector QBER
    time_filter = 2 # 2 sigma time filter on the coincidence window
    
    detector_loss = det_loss
    time = 1 # total time of prototocol
    
    file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
    file_night = 'data/HORIZONTAL_DISTANCE_NIGHT2_trx=15_L=%i_wl=0.81'




    

    if mode2 == 'decoyAO':
        background_day = 114538
        ext_loss = 0
    elif mode2 == 'esAO':
        background_day = 263437
        ext_loss = 0
    elif mode2 == 'decoy':
        background_day = 2319399
        ext_loss = 0
    elif mode2 == 'es' :
        background_day = 5334618
        ext_loss = 0
    elif mode2 == 'decoyET':
        background_day = 146344
        ext_loss = 2

    # Day or night settings
    if mode == 'night':
        bg_rate = background_day/50 # background rate
        file_data = file_night
    elif mode == 'day':
        bg_rate = background_day # background rate
        file_data = file_day
    
    
    
    
    confidence = 0.99
    
    channel = {
        'ra': ra,
        'channel_loss': 0,
        'detector_loss': detector_loss + ext_loss,
        'bg_rate': bg_rate,
        'jitter': jitter,
        'time_filter': time_filter,
        'time': time, 
        'e_d': e_d,
        'confidence': confidence}
    
    
    
    qs1 = []
    ss1 = []
    
    qs_ES = []
    ss_ES = []

ds = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]
# ds = [26000]
if calc:
    
    print(' ----- Background: ', bg_rate)
    for di in ds:
        # d = 10000
        d = di
        
        data = sp.io.loadmat(file_data %d)
        
        scattering_loss = 0.45 # db/km for 810
        data = data['res'].transpose()[0]
        data = -10*np.log10(data)+scattering_loss*d*1e-3
        
        # print()
        # print()
        print('----------------> d = ', d)
        # print('channel loss mean: ', np.mean(data))
        
        
        
        
        Nbins = 60
        interpol = interpolator.Interpolator(data, Nbins, plot=False)
        
        channel['channel_loss'] = interpol.get_sampler()
        
        opt = optimizer.Optimizer(channel)
        
        
        m = 0.5 # signal mean photon number
        v = 0.1 # decoy mean photon number  - required v < m
        pm = 0.1
        x = [m, v, pm] 
        
        finite = True
        
        res= opt.optimize(x)
        
        # print(res)
        
        xo = res['x']
        print('optimized values: ', xo)
        
        Qm, Qv, Y1, Y0, q_ests = opt.theoretical_vals_decoy(xo, finite)
        
        
        Nm = ra*time*Qm*xo[2]
        Nv = ra*time*Qv*(1-xo[2])
    
        b1 = ((Qm*np.exp(xo[0])-Y0)/(2*xo[0]) < Y1)
        b2 = ((Qv*np.exp(xo[1])-Y0)/(2*xo[1]) < Y1)
        
        
        s2 = b1 * Nm + b2 * Nv 
        q2 = q_ests[1]
        
        ss1 += [s2]
        qs1 += [q2]
        
        
        # print('------ ES')
        
        s_ES, q_ES = opt.values_ES()
        
        ss_ES += [s_ES]
        qs_ES += [q_ES]
        


    if mode == 'day':
        if mode2 == 'decoy':
            np.save('horizontal_s' + str(int(ras)) + '_day_dec.npy', ss1)
            np.save('horizontal_q' + str(int(ras)) + '_day_dec.npy', qs1)
            print(qs1)
            print(ss1)
        if mode2 == 'decoyAO':
            np.save('horizontal_s' + str(int(ras)) + '_day_decAO.npy', ss1)
            np.save('horizontal_q' + str(int(ras)) + '_day_decAO.npy', qs1)
            print(qs1)
            print(ss1)
        if mode2 == 'decoyET':
            np.save('horizontal_s' + str(int(ras)) + '_day_decET.npy', ss1)
            np.save('horizontal_q' + str(int(ras)) + '_day_decET.npy', qs1)
            print(qs1)
            print(ss1)
        elif mode2 == 'es':
            np.save('horizontal_s' + str(int(ras)) + '_day_es.npy', ss_ES)
            np.save('horizontal_q' + str(int(ras)) + '_day_es.npy', qs_ES)
            print(ss_ES)
            print(qs_ES)
        elif mode2 == 'esAO':
            np.save('horizontal_s' + str(int(ras)) + '_day_esAO.npy', ss_ES)
            np.save('horizontal_q' + str(int(ras)) + '_day_esAO.npy', qs_ES)
            print(ss_ES)
            print(qs_ES)
    elif mode == 'night':
        if mode2 == 'decoy':
            np.save('horizontal_s' + str(int(ras)) + '_night_dec.npy', ss1)
            np.save('horizontal_q' + str(int(ras)) + '_night_dec.npy', qs1)
            print(qs1)
            print(ss1)
        if mode2 == 'decoyAO':
            np.save('horizontal_s' + str(int(ras)) + '_night_decAO.npy', ss1)
            np.save('horizontal_q' + str(int(ras)) + '_night_decAO.npy', qs1)
            print(qs1)
            print(ss1)
        if mode2 == 'decoyET':
            np.save('horizontal_s' + str(int(ras)) + '_night_decET.npy', ss1)
            np.save('horizontal_q' + str(int(ras)) + '_night_decET.npy', qs1)
            print(qs1)
            print(ss1)
        elif mode2 == 'es':
            np.save('horizontal_s' + str(int(ras)) + '_night_es.npy', ss_ES)
            np.save('horizontal_q' + str(int(ras)) + '_night_es.npy', qs_ES)
            print(ss_ES)
            print(qs_ES)
        elif mode2 == 'esAO':
            np.save('horizontal_s' + str(int(ras)) + '_night_esAO.npy', ss_ES)
            np.save('horizontal_q' + str(int(ras)) + '_night_esAO.npy', qs_ES)
            print(ss_ES)
            print(qs_ES)



s20_day = np.load('horizontal_s2_day_dec.npy')
q20_day = np.load('horizontal_q2_day_dec.npy')
s20_night = np.load('horizontal_s2_night_dec.npy')
q20_night = np.load('horizontal_q2_night_dec.npy')


s20_dayAO = np.load('horizontal_s2_day_decAO.npy')
q20_dayAO = np.load('horizontal_q2_day_decAO.npy')
s20_nightAO = np.load('horizontal_s2_night_decAO.npy')
q20_nightAO = np.load('horizontal_q2_night_decAO.npy')

s20_dayET = np.load('horizontal_s2_day_decET.npy')
q20_dayET = np.load('horizontal_q2_day_decET.npy')
s20_nightET = np.load('horizontal_s2_night_decET.npy')
q20_nightET = np.load('horizontal_q2_night_decET.npy')


s100_day = np.load('horizontal_s10_day_dec.npy')
q100_day = np.load('horizontal_q10_day_dec.npy')
s100_night = np.load('horizontal_s10_night_dec.npy')
q100_night = np.load('horizontal_q10_night_dec.npy')

s100_dayAO = np.load('horizontal_s10_day_decAO.npy')
q100_dayAO = np.load('horizontal_q10_day_decAO.npy')
s100_nightAO = np.load('horizontal_s10_night_decAO.npy')
q100_nightAO = np.load('horizontal_q10_night_decAO.npy')


s100_dayET = np.load('horizontal_s10_day_decET.npy')
q100_dayET = np.load('horizontal_q10_day_decET.npy')
s100_nightET = np.load('horizontal_s10_night_decET.npy')
q100_nightET = np.load('horizontal_q10_night_decET.npy')

s20_day_es = np.load('horizontal_s1_day_es.npy')
q20_day_es = np.load('horizontal_q1_day_es.npy')
s20_night_es = np.load('horizontal_s1_night_es.npy')
q20_night_es = np.load('horizontal_q1_night_es.npy')

s20_day_esAO = np.load('horizontal_s1_day_esAO.npy')
q20_day_esAO = np.load('horizontal_q1_day_esAO.npy')
s20_night_esAO = np.load('horizontal_s1_night_esAO.npy')
q20_night_esAO = np.load('horizontal_q1_night_esAO.npy')

s40_day_es = np.load('horizontal_s4_day_es.npy')
q40_day_es = np.load('horizontal_q4_day_es.npy')
s40_night_es = np.load('horizontal_s4_night_es.npy')
q40_night_es = np.load('horizontal_q4_night_es.npy')


s40_day_esAO = np.load('horizontal_s4_day_esAO.npy')
q40_day_esAO = np.load('horizontal_q4_day_esAO.npy')
s40_night_esAO = np.load('horizontal_s4_night_esAO.npy')
q40_night_esAO = np.load('horizontal_q4_night_esAO.npy')



plt.rcParams["font.family"] = "Times New Roman"

cs = plt.rcParams['axes.prop_cycle'].by_key()['color']







#########################################

linetypes = ['o-', 's-', 'v-', '^-']


jitter = 345e-12 # total time jitter


s20_day = np.array(s20_day)/10
s100_day = np.array(s100_day)/10
s100_night = np.array(s100_night)/10
s20_night = np.array(s20_night)/10

s20_dayAO = np.array(s20_dayAO)/10
s100_dayAO = np.array(s100_dayAO)/10
s100_nightAO = np.array(s100_nightAO)/10
s20_nightAO = np.array(s20_nightAO)/10

s20_dayET = np.array(s20_dayET)/10
s100_dayET = np.array(s100_dayET)/10
s100_nightET = np.array(s100_nightET)/10
s20_nightET = np.array(s20_nightET)/10

t20_night = 2*jitter/np.sqrt(s20_night)*1e12
t20_day = 2*jitter/np.sqrt(s20_day)*1e12
t100_day = 2*jitter/np.sqrt(s100_day)*1e12
t100_night = 2*jitter/np.sqrt(s100_night)*1e12

t20_nightAO = 2*jitter/np.sqrt(s20_nightAO)*1e12
t20_dayAO = 2*jitter/np.sqrt(s20_dayAO)*1e12
t100_dayAO = 2*jitter/np.sqrt(s100_dayAO)*1e12
t100_nightAO = 2*jitter/np.sqrt(s100_nightAO)*1e12

t20_nightET = 2*jitter/np.sqrt(s20_nightET)*1e12
t20_dayET = 2*jitter/np.sqrt(s20_dayET)*1e12
t100_dayET = 2*jitter/np.sqrt(s100_dayET)*1e12
t100_nightET = 2*jitter/np.sqrt(s100_nightET)*1e12

t20_day = t20_day/np.sqrt(2)
t20_night = t20_night/np.sqrt(2)
t100_day = t100_day/np.sqrt(2)
t100_night = t100_night/np.sqrt(2)

t20_dayAO = t20_dayAO/np.sqrt(2)
t20_nightAO = t20_nightAO/np.sqrt(2)
t100_dayAO = t100_dayAO/np.sqrt(2)
t100_nightAO = t100_nightAO/np.sqrt(2)

t20_dayET = t20_dayET/np.sqrt(2)
t20_nightET = t20_nightET/np.sqrt(2)
t100_dayET = t100_dayET/np.sqrt(2)
t100_nightET = t100_nightET/np.sqrt(2)

s20_day_es = np.array(s20_day_es)/10
s20_night_es = np.array(s20_night_es)/10
s40_day_es = np.array(s40_day_es)/10
s40_night_es = np.array(s40_night_es)/10

s20_day_esAO = np.array(s20_day_esAO)/10
s20_night_esAO = np.array(s20_night_esAO)/10
s40_day_esAO = np.array(s40_day_esAO)/10
s40_night_esAO = np.array(s40_night_esAO)/10

t20_day_es = 2*jitter/np.sqrt(s20_day_es)*1e12
t20_night_es = 2*jitter/np.sqrt(s20_night_es)*1e12
t40_day_es = 2*jitter/np.sqrt(s40_day_es)*1e12
t40_night_es = 2*jitter/np.sqrt(s40_night_es)*1e12

t20_day_esAO = 2*jitter/np.sqrt(s20_day_esAO)*1e12
t20_night_esAO = 2*jitter/np.sqrt(s20_night_esAO)*1e12
t40_day_esAO = 2*jitter/np.sqrt(s40_day_esAO)*1e12
t40_night_esAO = 2*jitter/np.sqrt(s40_night_esAO)*1e12

t20_day_es = t20_day_es/np.sqrt(2)
t20_night_es = t20_night_es/np.sqrt(2)
t40_day_es = t40_day_es/np.sqrt(2)
t40_night_es = t40_night_es/np.sqrt(2)

t20_day_esAO = t20_day_esAO/np.sqrt(2)
t20_night_esAO = t20_night_esAO/np.sqrt(2)
t40_day_esAO = t40_day_esAO/np.sqrt(2)
t40_night_esAO = t40_night_esAO/np.sqrt(2)


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.2)
ax.set_title('dec prec')
ax.set_ylabel('Precision (ps)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
ax.plot(ds[:len(t20_night)], t20_night, linetypes[0], color=cs[0], label='20 MHz - night', linewidth=3)
ax.plot(ds[:len(t100_night)], t100_night, linetypes[1], color=cs[2], label='100 MHz - night', linewidth=3)
ax.plot(ds[:len(t20_day)], t20_day, linetypes[2], color=cs[1], label='20 MHz - day', linewidth=3)
ax.plot(ds[:len(t100_day)], t100_day, linetypes[3], color=cs[3], label='100 MHz - day', linewidth=3)

ax.set_xlim([6500, 30000])
ax.set_ylim([0, 50])

plt.tight_layout()
fig.savefig('time_precision_decoy1.pdf', dpi=200)

fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.2)
ax.set_title('dec prec AO')
ax.set_ylabel('Precision (ps)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
ax.plot(ds[:len(t20_night)], t20_nightAO, linetypes[0], color=cs[0], label='20 MHz AO - night', linewidth=3)
ax.plot(ds[:len(t100_night)], t100_nightAO, linetypes[1], color=cs[2], label='100 MHz AO - night', linewidth=3)
ax.plot(ds[:len(t20_day)], t20_dayAO, linetypes[2], color=cs[1], label='20 MHz AO - day', linewidth=3)
ax.plot(ds[:len(t100_day)], t100_dayAO, linetypes[3], color=cs[3], label='100 MHz AO - day', linewidth=3)

ax.set_xlim([6500, 30000])
ax.set_ylim([0, 50])

plt.tight_layout()
fig.savefig('time_precision_decoy2.pdf', dpi=200)

fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.2)
ax.set_title('dec prec ET')
ax.set_ylabel('Precision (ps)', fontsize=15)
# ax.set_xlabel('Distance (m)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
ax.plot(ds[:len(t20_night)], t20_nightET, linetypes[0], color=cs[0], label='Night - 20 MHz', linewidth=3)
ax.plot(ds[:len(t100_night)], t100_nightET, linetypes[1], color=cs[2], label='Night - 100 MHz', linewidth=3)
ax.plot(ds[:len(t20_day)], t20_dayET, linetypes[2], color=cs[1], label='Day - 20 MHz', linewidth=3)
ax.plot(ds[:len(t100_day)], t100_dayET, linetypes[3], color=cs[3], label='Day - 100 MHz', linewidth=3)
plt.legend(fontsize=15)

ax.set_xlim([6500, 30000])
ax.set_ylim([0, 50])
plt.tight_layout()
fig.savefig('time_precision_decoy3.pdf', dpi=200)
# ax.set_xlabel('Distance (m)', fontsize=15)
# ax.set_ylabel('Precision (ps)', fontsize=15)

# plt.xticks(fontsize=13, rotation=0)
# plt.yticks(fontsize=13, rotation=0)

# ax.set_xlim([6500, 30000])
# plt.tight_layout()

# fig.savefig('time_precision_decoy.pdf', dpi=200)



fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.2)
ax.set_title('dec qber')
ax.set_ylabel('QBER', fontsize=15)
ax.set_ylim([0.00, 0.25])
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)
ax.set_ylim([0.00, 0.25])
ax.set_xlim([6500, 30000])
ax.plot(ds[:len(q20_night)], q20_night, linetypes[0], color=cs[0], label='Night - 20 MHz', linewidth=3)
ax.plot(ds[:len(q100_night)], q100_night, linetypes[1], color=cs[2], label='Night - 100 MHz', linewidth=3)
ax.plot(ds[:len(q20_day)], q20_day, linetypes[2], color=cs[1], label='Day - 20 MHz', linewidth=3)
ax.plot(ds[:len(q100_night)], q100_day, linetypes[3], color=cs[3], label='Day - 100 MHz', linewidth=3)
plt.tight_layout()
fig.savefig('qbers_decoy1.pdf', dpi=200)


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.2)
ax.set_title('dec qber AO')
ax.set_ylabel('QBER', fontsize=15)
ax.set_ylim([0.00, 0.25])
ax.set_xlim([6500, 30000])

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
ax.plot(ds[:len(q20_night)], q20_nightAO, linetypes[0], color=cs[0], label='Night AO - 20 MHz', linewidth=3)
ax.plot(ds[:len(q100_night)], q100_nightAO, linetypes[1], color=cs[2], label='Night AO - 100 MHz', linewidth=3)
ax.plot(ds[:len(q20_day)], q20_dayAO, linetypes[2], color=cs[1], label='Day AO - 20 MHz', linewidth=3)
ax.plot(ds[:len(q100_night)], q100_dayAO, linetypes[3], color=cs[3], label='Day AO - 100 MHz', linewidth=3)
plt.tight_layout()
fig.savefig('qbers_decoy2.pdf', dpi=200)

fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.2)
ax.set_title('dec qber ET')
ax.set_ylabel('QBER', fontsize=15)
# ax.set_xlabel('Distance (m)', fontsize=15)
ax.set_ylim([0.00, 0.25])
ax.set_xlim([6500, 30000])

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
ax.plot(ds[:len(q20_night)], q20_nightET, linetypes[0], color=cs[0], label='Night - 20 MHz', linewidth=3)
ax.plot(ds[:len(q100_night)], q100_nightET, linetypes[1], color=cs[2], label='Night - 100 MHz', linewidth=3)
ax.plot(ds[:len(q20_day)], q20_dayET, linetypes[2], color=cs[1], label='Day - 20 MHz', linewidth=3)
ax.plot(ds[:len(q100_night)], q100_dayET, linetypes[3], color=cs[3], label='Day - 100 MHz', linewidth=3)
plt.tight_layout()
fig.savefig('qbers_decoy3.pdf', dpi=200)

# ax.set_xlabel('Distance (m)', fontsize=15)
# ax.set_ylabel('QBER', fontsize=15)

# plt.xticks(fontsize=13, rotation=0)
# plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)

# ax.set_ylim([0.00, 0.25])
# ax.set_xlim([6500, 30000])

# plt.tight_layout()

# fig.savefig('qbers_decoy.pdf', dpi=200)






fig = plt.figure()
plt.title('ES QBER')
plt.ylim([0.0, 0.25])
fig.set_size_inches(18.5*.3, 10.5*.2)
plt.xlim([6500, 30000])
plt.plot(ds, q20_night_es, linetypes[0], color=cs[0], label='Night - 10 MHz', linewidth=3)
plt.plot(ds, q40_night_es, linetypes[1], color=cs[2], label='Night - 40 MHz', linewidth=3)
plt.plot(ds, q20_day_es, linetypes[2], color=cs[1], label='Day - 10 MHz', linewidth=3)
plt.plot(ds, q40_day_es, linetypes[3], color=cs[3], label='Day - 40 MHz', linewidth=3)
# plt.xlabel('Distance (m)', fontsize=15)
plt.ylabel('QBER', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

plt.tight_layout()
fig.savefig('qber_main1.pdf', dpi=200)


fig = plt.figure()
plt.title('ES AO QBER')
plt.ylim([0.0, 0.25])
fig.set_size_inches(18.5*.3, 10.5*.2)
plt.xlim([6500, 30000])
plt.plot(ds, q20_night_esAO, linetypes[0], color=cs[0], label='Night - 40 MHz', linewidth=3)
plt.plot(ds, q40_night_esAO, linetypes[1], color=cs[2], label='Night - 1 MHz', linewidth=3)
plt.plot(ds, q20_day_esAO, linetypes[2], color=cs[1], label='Day - 40 MHz', linewidth=3)
plt.plot(ds, q40_day_esAO, linetypes[3], color=cs[3], label='Day - 20 MHz', linewidth=3)
# plt.xlabel('Distance (m)', fontsize=15)
plt.ylabel('QBER', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.tight_layout()
fig.savefig('qber_main2.pdf', dpi=200)


# plt.xlabel('Distance (m)', fontsize=15)
# plt.ylabel('QBER', fontsize=15)
# plt.xticks(fontsize=13, rotation=0)
# plt.yticks(fontsize=13, rotation=0)
# # plt.legend(fontsize=15)
# plt.ylim([0.0, 0.25])

# fig.set_size_inches(18.5*.3, 10.5*.2)
# plt.xlim([6500, 30000])

# plt.tight_layout()

# fig.savefig('qber_main.pdf', dpi=200)



fig = plt.figure()
plt.title('ES PREC')
# plt.xlabel('Distance (m)', fontsize=15)
plt.ylabel('Precision (ps)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)
fig.set_size_inches(18.5*.3, 10.5*.2)
plt.ylim([0.0, 70])
plt.xlim([6500, 30000])
plt.plot(ds, t20_night_es, linetypes[0], color=cs[0], label='Night - 20 MHz', linewidth=3)
plt.plot(ds, t40_night_es, linetypes[1], color=cs[2], label='Night - 10 MHz', linewidth=3)
plt.plot(ds, t20_day_es, linetypes[2], color=cs[1], label='Day - 20 MHz', linewidth=3)
plt.plot(ds, t40_day_es, linetypes[3], color=cs[3], label='Day - 40 MHz', linewidth=3)
plt.tight_layout()
fig.savefig('time_precision1.pdf', dpi=200)


fig = plt.figure()
plt.title('ES PREC AO')
# plt.xlabel('Distance (m)', fontsize=15)
plt.ylabel('Precision (ps)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)
fig.set_size_inches(18.5*.3, 10.5*.2)
plt.ylim([0.0, 70])
plt.xlim([6500, 30000])
plt.plot(ds, t20_night_esAO, linetypes[0], color=cs[0], label='Night - 20 MHz', linewidth=3)
plt.plot(ds, t40_night_esAO, linetypes[1], color=cs[2], label='Night - 10 MHz', linewidth=3)
plt.plot(ds, t20_day_esAO, linetypes[2], color=cs[1], label='Day - 20 MHz', linewidth=3)
plt.plot(ds, t40_day_esAO, linetypes[3], color=cs[3], label='Day - 40 MHz', linewidth=3)
plt.tight_layout()
fig.savefig('time_precision2.pdf', dpi=200)

# plt.xlabel('Distance (m)', fontsize=15)
# plt.ylabel('Precision (ps)', fontsize=15)
# plt.xticks(fontsize=13, rotation=0)
# plt.yticks(fontsize=13, rotation=0)
# # plt.legend(fontsize=15)

# fig.set_size_inches(18.5*.3, 10.5*.2)
# plt.ylim([0.0, 100])
# plt.xlim([6500, 30000])

# plt.tight_layout()

# fig.savefig('time_precision.pdf', dpi=200)