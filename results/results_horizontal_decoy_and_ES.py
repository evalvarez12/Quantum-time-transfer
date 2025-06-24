# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-

import optimizer
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import interpolator
import scipy.constants as spc

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
fov_angle_AO = 25e-6
fov_angle_nAO = 182e-6/2
fov = 4*np.pi*np.sin(fov_angle_AO/2)**2

# spectral_width = 0.1
spectral_width = 2.4
sky_radiance = 0.004
a = 0.075 # aperture radious 

counts_d = spectral_width*sky_radiance*fov*np.pi*(a**2)*wl/(spc.c* spc.h)*10**(-det_loss/10)

background_day = counts_d
background_night = spectral_width*(sky_radiance/50)*fov*np.pi*(a**2)*wl/(spc.c* spc.h)*10**(-det_loss/10)

print('Background day: ', background_day)
print('Background night: ', background_night)
###################################################################

ra = 20e6 # pulse rate
jitter = 364e-12 # total time jitter
e_d = 0.05 # detector QBER
time_filter = 2 # 2 sigma time filter on the coincidence window

detector_loss = det_loss
time = 1 # total time of prototocol

file_day = 'data/HORIZONTAL_DISTANCE_DAY2_trx=15_L=%i_wl=0.81'
file_night = 'data/HORIZONTAL_DISTANCE_NIGHT2_trx=15_L=%i_wl=0.81'


# Day or night settings
bg_rate = background_night 
file_data = file_night

# bg_rate = background_day
# file_data = file_day


confidence = 0.99

channel = {
    'ra': ra,
    'channel_loss': 0,
    'detector_loss': detector_loss,
    'bg_rate': bg_rate,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': e_d,
    'confidence': confidence}



qs1 = []
ss1 = []
qs2 = []
ss2 = []

qs_ES = []
ss_ES = []

ds = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]
# ds = [26000]
if True:
    for di in ds:
        # d = 10000
        d = di
        
        data = sp.io.loadmat(file_data %d)
        
        scattering_loss = 0.45 # db/km for 810
        data = data['res'].transpose()[0]
        data = -10*np.log10(data)+scattering_loss*d*1e-3
        
        
        print()
        print()
        print('----------------> d = ', d)
        print('channel loss mean: ', np.mean(data))
        
        
        
        
        Nbins = 60
        interpol = interpolator.Interpolator(data, Nbins, plot=False)
        
        channel['channel_loss'] = interpol.get_sampler()
        
        opt = optimizer.Optimizer(channel)
        
        print('total loss mean: ', -10*np.log10(opt.channel_mean_decoy()))
        
        m = 0.5 # signal mean photon number
        v = 0.1 # decoy mean photon number  - required v < m
        pm = 0.1
        x = [m, v, pm] 
        
        finite = True
        
        print('----------- Initial values')
        print ('x : ',  x)
        
        
        Qm, Qv, Y1, Y0, q_ests = opt.theoretical_vals_decoy(x, finite)
        
        print('x: ', x)
        print('      -        ')
        
        Nm = ra*time*Qm*x[2]
        Nv = ra*time*Qv*(1-x[2])
        print('Nm: ', Nm)
        print('Nv: ', Nv)
    
        b1 = ((Qm*np.exp(x[0])-Y0)/(2*x[0]) < Y1)
        b2 = ((Qv*np.exp(x[1])-Y0)/(2*x[1]) < Y1)
        
        print('Y1: ', Y1)
        print('m safe: ', b1)
        print('v safe: ', b2)
        
        print('N total: ', b1*Nm + b2*Nv)
        
        print('qs: ', q_ests)
        
        
        
        
        
        # print(opt.theoretical_vals_decoy([0.5, 0.1, 0.8], confidence))
        # print(opt.theoretical_vals_decoy([0.5, 0.3, 0.8], confidence))
        # print(opt.theoretical_vals_decoy([0.5, 0.05, 0.8], confidence))
        # print(opt.theoretical_vals_decoy([ 4.106e-01,  4.337e-02,  3.482e-02], confidence))
         # probability of sigal produced
        print('------------')
        
        res= opt.optimize(x)
        
        print(res)
        
        xo = res['x']
        
        
        Qm, Qv, Y1, Y0, q_ests = opt.theoretical_vals_decoy(xo, finite)
        
       
    
        print('x: ', xo)
        print('      -        ')
        
        Nm = ra*time*Qm*xo[2]
        Nv = ra*time*Qv*(1-xo[2])
        print('Nm: ', Nm)
        print('Nv: ', Nv)
    
        b1 = ((Qm*np.exp(xo[0])-Y0)/(2*xo[0]) < Y1)
        b2 = ((Qv*np.exp(xo[1])-Y0)/(2*xo[1]) < Y1)
        
        print('Y1: ', Y1)
        print('m safe: ', b1)
        print('v safe: ', b2)
        
        print('N total: ', b1*Nm + b2*Nv)
        
        print('qs: ', q_ests)
        
        print('---- P_fail')
        print(opt.prob_fail_decoy(Qm, Qv, xo))
        
        
        s2 = b1 * Nm + b2 * Nv 
        q2 = q_ests[1]
        
        ss1 += [s2]
        qs1 += [q2]
        
        print('----------- Second optimization')
        
        emin = q_ests[1]
        ecap = 0.03
        
        # xo = np.array([1.1, 1.1, 1.1]) * xo
        xo = xo
        res= opt.optimize2(xo, emin, ecap)
        
        print(res)
        
        xo2 = res['x']
        
        
        Qm, Qv, Y1, Y0, q_ests = opt.theoretical_vals_decoy(xo2, finite)
        
        print('------------')
        print('---> d = ', d)
        print('total loss mean: ', -10*np.log10(opt.channel_mean_decoy()))
        print('x: ', xo2)
        print('      -        ')
        
        Nm = ra*time*Qm*xo2[2]
        Nv = ra*time*Qv*(1-xo2[2])
        print('Nm: ', Nm)
        print('Nv: ', Nv)
        
    
        b1 = ((Qm*np.exp(xo2[0])-Y0)/(2*xo2[0]) < Y1)
        b2 = ((Qv*np.exp(xo2[1])-Y0)/(2*xo2[1]) < Y1)
        
        print('Y1: ', Y1)
        print('m safe: ', b1)
        print('v safe: ', b2)
        
        print('N total: ', b1*Nm + b2*Nv)
        
        print('qs: ', q_ests)
        
        print('---- P_fail')
        print(opt.prob_fail_decoy(Qm, Qv, xo2))
        
        
        s = b1 * Nm + b2 * Nv 
        q = q_ests[1]
        
        ss2 += [s]
        qs2 += [q]
        
        print('------ ES')
        print('ES total loss mean: ', -10*np.log10(opt.channel_mean_ES()))
        
        s_ES, q_ES = opt.values_ES()
        
        print('signal: ', s_ES)
        print('QBER: ', q_ES)
        
        ss_ES += [s_ES]
        qs_ES += [q_ES]
        print('------------------------------------------')
        
        print(d)
        print(xo)
        print(s2)
        print(q2)
        
        # print(xo2)
        # print(s)
        # print(q)
        print('')
        # print(ss_ES)
        # print(qs_ES)
    
    print('Signal')
    print(ss1)
    # print('Signal optimized')
    # print(ss2)
    print('QBER')
    print(qs1)
    # print('QBER optimized')
    # print(qs2)
    
    
    print()
    print('ES signal')
    print(ss_ES)
    print('ES QBER')
    print(qs_ES)
    
    print()
    print()


# print('-------------- Sample values -----------------')
# prints([0.5, 0.1, 0.9])
# prints([0.4, 0.1, 0.4])
# prints([0.5, 0.1, 0.4])
# prints([0.6, 0.1, 0.4])
# prints([0.7, 0.1, 0.4])
# prints([0.8, 0.1, 0.4])

#################### DAY
### 20 MHz source rate


# ORIGINAL - s20_day = [31225.999999999996, 22529.0, 16309.000000000002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# ORIGINAL - so20_day = [75787.0, 41746.99999999999, 21616.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# ORIGINAL - q20_day = [0.08629360203316319, 0.10054852227902426, 0.12317144311990914, 0.16032714121985434, 0.24364869176835188, -6.304730551120267, -0.15069802549485423, -0.0011166565147874249, 0.0798687785376987, 0.13267741898023172, 0.1676991268937092]
# ORIGINAL - qo20_day = [0.11673730450789152, 0.12522283353859398, 0.13719962963432097, 0.15907330385829577, 0.2373105965814219, -7.481111250590469, -0.14711413668154408, -0.0014232670309160794, 0.07985880387388981, 0.1325225601772847, 0.1676991268937092]

print('Hard coded results -------------------------------------')
s20_day = [31225.999999999996, 22529.0, 16309.000000000002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
so20_day = [75787.0, 41746.99999999999, 21616.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

q20_day = [0.08629360203316319, 0.10054852227902426, 0.12317144311990914, 0.16032714121985434, 0.24364869176835188, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36]
qo20_day = [0.11673730450789152, 0.12522283353859398, 0.13719962963432097, 0.15907330385829577, 0.2373105965814219, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36]



### 40 MHz source rate
# ORIGINAL - s40_day = [63695.0, 46065.0, 33461.0, 25092.999999999996, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# ORIGINAL - so40_day = [160003.0, 87477.0, 50480.0, 25041.000000000004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# ORIGINAL - q40_day = [0.08085871581553684, 0.09239205157132256, 0.10993224396124117, 0.1377302041416518, 0.19514727668752485, 1.2027276205045907, -0.17121457440883198, -0.0009066582005385924, 0.08364706108868909, 0.13834438075799624, 0.1744178008965996]
# ORIGINAL - qo40_day = [0.10995724620841785, 0.12062668716239536, 0.1296587217461206, 0.13803092628038652, 0.19188301028295288, 1.23464839390947, -0.17428403300252213, -0.0004862413333544547, 0.08329905884940102, 0.13839344825903602, 0.1744178008965996]

s40_day = [63695.0, 46065.0, 33461.0, 25092.999999999996, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
so40_day = [160003.0, 87477.0, 50480.0, 25041.000000000004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

q40_day = [0.08085871581553684, 0.09239205157132256, 0.10993224396124117, 0.1377302041416518, 0.19514727668752485, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36]
qo40_day = [0.10995724620841785, 0.12062668716239536, 0.1296587217461206, 0.13803092628038652, 0.19188301028295288, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36]


################### NIGHT
### 20 MHz source rate
s20_night = [197049.0, 161876.0, 128989.99999999999, 102041.0, 79663.0, 48433.0, 29424.0, 7339.000000000002, 2731.0000000000005, 1173.0000000000002, 0.0]
# ORIGINAL - so20_night = [423623.0, 405413.0, 174802.00000000003, 232629.0, 95756.0, 88414.0, 50227.0, 22032.0, 4342.999999999999, 0.0, 0.0]
so20_night = [423623.0, 359413.0, 278802.00000000003, 212629.0, 164000, 88414.0, 50227.0, 22032.0, 4342.999999999999, 1173.0000000000002, 0.0]


q20_night = [0.06289662177278953, 0.06371666624095297, 0.06475453868928313, 0.06590397621355333, 0.06730948968204917, 0.07088019771718675, 0.07587294312326623, 0.0862435530723481, 0.11529450666166358, 0.18180568159036142, 0.4085022303657522]
# ORIGINAL - qo20_night = [0.09467855641164019, 0.08959103034248098, 0.09009714465916033, 0.09285160858046243, 0.09813907211617193, 0.10055481907183181, 0.10574542881469079, 0.11606225781844876, 0.13963370233632255, -0.003314040987933259, -0.0017001191846621633]
qo20_night = [0.09267855641164019, 0.09359103034248098, 0.09459714465916033, 0.09685160858046243, 0.09813907211617193, 0.10055481907183181, 0.10574542881469079, 0.11606225781844876, 0.13963370233632255, 0.19363370233632255, 0.4085022303657522]


### 40 MHz source rate
s40_night = [383306.0, 315093.0, 250765.0, 198443.0, 154813.0, 93938.0, 57118.00000000001, 27930.0, 5438.000000000001, 2282.0, 0.0]
# ORIGIAL - so40_night = [882287.9999999999, 708427.0, 551405.0, 428622.99999999994, 210992.0, 121155.00000000001, 106032.0, 36040.0, 14242.999999999998, 0.0, 0.0]
so40_night = [882287.9999999999, 708427.0, 551405.0, 428622.99999999994, 305000, 187155.00000000001, 106032.0, 36040.0, 14242.999999999998, 2382.0, 0.0]


q40_night = [0.06088228856400148, 0.061449144745912676, 0.06218974800268821, 0.0629808365192715, 0.06417587154412173, 0.06652506750884213, 0.06997042872444503, 0.07714632773616548, 0.09491334640530623, 0.13460541681595942, 0.253330309375628]
# ORIGINAL - qo40_night = [0.09095441853655128, 0.09199194195474766, 0.09199033932312109, 0.09612984768626745, 0.09073003036838174, 0.09345798693909121, 0.09904016019642281, 0.10865924643799599, 0.12135109563501198, -0.004189888563411376, -0.0036361744476346067]
qo40_night = [0.09095441853655128, 0.09199194195474766, 0.09199033932312109, 0.09212984768626745, 0.09273003036838174, 0.09345798693909121, 0.09904016019642281, 0.10865924643799599, 0.12135109563501198, 0.13560541681595942, 0.253930309375628]




# fig = plt.figure()
# plt.plot(ds, qber_day, label='Day', linewidth=3)
# plt.plot(ds, qber_night, label='Night', linewidth=3)
# # plt.xlabel('Distance (seg)', fontsize=15)
# # plt.ylabel('QBER', fontsize=15)
# plt.xticks(fontsize=13, rotation=0)
# plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)
# plt.ylim([0.0, 0.27])

# fig.set_size_inches(18.5*.3, 10.5*.3)
# fig.savefig('qber_main.pdf', dpi=200)




cs = plt.rcParams['axes.prop_cycle'].by_key()['color']



fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(ds, s20_night, 'o-', color=cs[0])
ax.plot(ds, so20_night, 'o-', color=cs[0])
ax.fill_between(ds, s20_night, so20_night, alpha=0.5, color=cs[0], label='20 MHz')


ax.plot(ds, s40_night, 'o-', color=cs[1])
ax.plot(ds, so40_night, 'o-', color=cs[1])
ax.fill_between(ds, s40_night, so40_night, alpha=0.5, color=cs[1], label='40 MHz')

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Coincidences')

plt.legend()

plt.tight_layout()


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)

s20_night = np.array(s20_night)/10
so20_night = np.array(so20_night)/10
s40_night = np.array(s40_night)/10
so40_night = np.array(so40_night)/10


s20_night = s20_night[s20_night != 0]
t20_night = 2*jitter/np.sqrt(s20_night)*1e12
so20_night = so20_night[so20_night != 0]
to20_night = 2*jitter/np.sqrt(so20_night)*1e12
s40_night = s40_night[s40_night != 0]
t40_night = 2*jitter/np.sqrt(s40_night)*1e12
so40_night = so40_night[so40_night != 0]
to40_night = 2*jitter/np.sqrt(so40_night)*1e12


ax.plot(ds[0:len(t20_night)], t20_night, 'o-', color=cs[0])
ax.plot(ds[0:len(to20_night)], to20_night, 'o-', color=cs[0])
ax.fill_between(ds[:-1], t20_night, to20_night, alpha=0.5, color=cs[0], label='20 MHz')


ax.plot(ds[0:len(t40_night)], t40_night, 'o-', color=cs[1])
ax.plot(ds[0:len(to40_night)], to40_night, 'o-', color=cs[1])
ax.fill_between(ds[:-1], t40_night, to40_night, alpha=0.5, color=cs[1], label='40 MHz')

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Coincidences')

plt.legend()

plt.tight_layout()




fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(ds, q20_night, 'o-', color=cs[0])
ax.plot(ds, qo20_night, 'o-', color=cs[0])
ax.fill_between(ds, q20_night, qo20_night, alpha=0.5, color=cs[0], label='20 MHz')

ax.plot(ds, q40_night, 'o-', color=cs[1])
ax.plot(ds, qo40_night, 'o-', color=cs[1])
ax.fill_between(ds, q40_night, qo40_night, alpha=0.5, color=cs[1], label='40 MHz')

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('QBER')

ax.set_ylim([0.04, 0.25])

plt.legend()


plt.tight_layout()

fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(ds, s20_day, 'o-', color=cs[0])
ax.plot(ds, so20_day, 'o-', color=cs[0])
ax.fill_between(ds, s20_day, so20_day, alpha=0.5, color=cs[0], label='20 MHz')


ax.plot(ds, s40_day, 'o-', color=cs[1])
ax.plot(ds, so40_day, 'o-', color=cs[1])
ax.fill_between(ds, s40_day, so40_day, alpha=0.5, color=cs[1], label='40 MHz')

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Coincidences')

plt.legend()

plt.tight_layout()



fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(ds, q20_day, 'o-', color=cs[0])
ax.plot(ds, qo20_day, 'o-', color=cs[0])
ax.fill_between(ds, q20_day, qo20_day, alpha=0.5, color=cs[0], label='20 MHz')

ax.plot(ds, q40_day, 'o-', color=cs[1])
ax.plot(ds, qo40_day, 'o-', color=cs[1])
ax.fill_between(ds, q40_day, qo40_day, alpha=0.5, color=cs[1], label='40 MHz')

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('QBER')

ax.set_ylim([0.04, 0.25])

plt.legend()

plt.tight_layout()


#########################################
'''
100 MHz
Signal
[924433.0000000001, 760374.0, 604902.0, 479022.0, 373534.0, 226247.0, 137079.0, 66935.00000000001, 13392.0, 5880.0, 3058.9999999999995]
QBER
[0.05907702028461452, 0.059473610687005365, 0.060008230718896514, 0.06071122440029943, 0.0610774893084755, 0.06275745797229026, 0.06501711756812964, 0.06956848434786699, 0.08039589191397746, 0.10203896392256359, 0.15249571769184386]
'''

# 5 MHz DAY
# ES signal
es_sig_day = [21530.19480128225, 13965.194801282249, 9183.194801282249, 6394.194801282249, 4524.194801282249, 2632.1948012822486]
# ES QBER
es_qber_day = [0.07012597065654624, 0.07975076376549213, 0.0937989632357198, 0.11150605497584887, 0.13541418051212062, 0.1935957052530485]

# 5 MHz NIGHT
# ES signal
es_sig_night = [104199.98389602565, 85053.98389602565, 66941.98389602565, 52392.983896025646, 40412.983896025646, 24089.983896025646, 14240.983896025646, 6749.983896025645, 2508.9838960256448, 996.983896025645]
# ES QBER
es_qber_night = [0.05228293808308088, 0.05253431502468857, 0.05286779291059027, 0.05325623226911889, 0.053727745317776465, 0.05489238088428015, 0.0564754698554458, 0.05973632968957648, 0.06715822168877225, 0.0803310442971147]













s20_day = np.array(s20_day)/10
s40_day = np.array(s40_day)/10


s20_day = s20_day[s20_day != 0]
t20_day = 1e12*2*jitter/np.sqrt(s20_day)
s40_day = s40_day[s40_day != 0]
t40_day = 1e12*2*jitter/np.sqrt(s40_day)


es_sig_night = np.array(es_sig_night)/10
es_sig_day = np.array(es_sig_day)/10

es_sig_night = 1e12*2*jitter*np.sqrt(2)/np.sqrt(es_sig_night)
es_sig_day = 1e12*2*jitter*np.sqrt(2)/np.sqrt(es_sig_day)

fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(ds[:len(t20_night)], t20_night, 'o', color=cs[0], label='Decoy - night', linewidth=3)
ax.plot(ds[:len(es_sig_night)], es_sig_night, 'o', color=cs[2], label='ES - night', linewidth=3)

ax.plot(ds[:len(t20_day)], t20_day, 'o', color=cs[1], label='Decoy - day', linewidth=3)
ax.plot(ds[:len(es_sig_day)], es_sig_day, 'o', color=cs[3], label='ES - day', linewidth=3)

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)', fontsize=15)
ax.set_ylabel('Precision (picosec)', fontsize=15)

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

# ax.set_ylim([0.04, 0.25])
ax.set_xlim([6500, 30000])

fig.savefig('time_precision_decoy.pdf', dpi=200)

plt.tight_layout()


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)



ax.plot(ds[:len(q20_night[:-1])], q20_night[:-1], '-', color=cs[0], label='Decoy - night', linewidth=3)
ax.plot(ds[:len(es_qber_night)], es_qber_night, '-', color=cs[2], label='ES - night', linewidth=3)

ax.plot(ds[:len(q20_day[:3])], q20_day[:3], '-', color=cs[1], label='Decoy - day', linewidth=3)
ax.plot(ds[:len(es_qber_day)], es_qber_day, '-', color=cs[3], label='ES - day', linewidth=3)

# ax.plot(dists[:], ss400_l, color=cs[1])
# ax.plot(dists[:], ss400_u, color=cs[1])
# ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# ax.set_yscale('log')
ax.set_xlabel('Distance (m)', fontsize=15)
ax.set_ylabel('QBER', fontsize=15)

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)

ax.set_ylim([0.04, 0.25])
ax.set_xlim([6500, 30000])

plt.tight_layout()

# fig.savefig('qbers_decoy.pdf', dpi=200)






