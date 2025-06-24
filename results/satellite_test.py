
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import satellite

import warnings
warnings.filterwarnings("ignore")


ra = 100e6
sat_link = satellite.Satellite(ra=ra, time=1)

mode = 'up'
mode2 = 'decoy'

zs = [0, 10, 20, 30, 40, 50, 60, 70, 80]

up_sig_ES = []
up_qber_ES = []
up_sig_decoy = []
up_qber_decoy = []
down_sig_ES = []
down_qber_ES = []
down_sig_decoy = []
down_qber_decoy = []

if 0:
    for z in zs:
        print(f'-------> z = {z}')
        
        extT = [9, 6] # extra z = 0
        sat_link.setup_channels(z, True, extT)
        
        
        up_ES, down_ES = sat_link.get_vals_ES()
        
        
        m = 0.5 # signal mean photon number
        v = 0.3 # decoy mean photon number  - required v < m
        pm = 0.2
        x0 = [m, v, pm] 
        up_decoy, down_decoy = sat_link.get_vals_decoy(x0, x0)
        
        # print('General ')
        # print('Uplink: ', sat_link.ra*10**(-23.24/10))
        # print('Downlink: ', sat_link.ra*10**(-17.73/10))
        
        if mode == 'up':
            print('UPLINK -----------------------------------------')
            
            if mode2 == 'es':
                print('--- ES')
                print(up_ES)
                
                up_sig_ES += [up_ES[0]]
                up_qber_ES += [up_ES[1]]
                
            if mode2 == 'decoy':
                print('--- Decoy')
                print(up_decoy)
                
                if up_decoy['success']:
                    xo = up_decoy['x']
                    
                    Qm, Qv, Y1, Y0, q_ests = sat_link.uplink.theoretical_vals_decoy(xo)
                    
                    
                    Nm = sat_link.ra*sat_link.time*Qm*xo[2]
                    Nv = sat_link.ra*sat_link.time*Qv*(1-xo[2])
                    print('Nm: ', Nm)
                    print('Nv: ', Nv)
                    print('N total: ', Nm + Nv)
                
                    b1 = ((Qm*np.exp(xo[0])-Y0)/(2*xo[0]) < Y1)
                    b2 = ((Qv*np.exp(xo[1])-Y0)/(2*xo[1]) < Y1)
                    
                    print('Y1: ', Y1)
                    print('m safe: ', b1)
                    print('v safe: ', b2)
                    
                    up_sig_decoy += [b1*Nm + b2*Nv]
                    up_qber_decoy += [q_ests[1]]
                    
                    print('qs: ', q_ests)
                
        
        
        if mode == 'down':
            print('DOWNLINK -----------------------------------------')
            if mode2 == 'es':
                print('--- ES')
                print(down_ES)
                
                down_sig_ES += [down_ES[0]]
                down_qber_ES += [down_ES[1]]
                
            if mode2 == 'decoy':
                print('--- Decoy')
                print(down_decoy)
                
                if down_decoy['success']:
                    xo = down_decoy['x']
                    
                    Qm, Qv, Y1, Y0, q_ests = sat_link.downlink.theoretical_vals_decoy(xo)
                    
                    
                    Nm = sat_link.ra*sat_link.time*Qm*xo[2]
                    Nv = sat_link.ra*sat_link.time*Qv*(1-xo[2])
                    print('Nm: ', Nm)
                    print('Nv: ', Nv)
                    print('N total: ', Nm + Nv)
                
                    b1 = ((Qm*np.exp(xo[0])-Y0)/(2*xo[0]) < Y1)
                    b2 = ((Qv*np.exp(xo[1])-Y0)/(2*xo[1]) < Y1)
                    
                    print('Y1: ', Y1)
                    print('m safe: ', b1)
                    print('v safe: ', b2)
                    
                    down_sig_decoy += [b1*Nm + b2*Nv]
                    down_qber_decoy += [q_ests[1]]
                    
                    print('qs: ', q_ests)
                    
                
                    
    if mode2 == 'es':
        if mode == 'down':
            print('ES_sig_d = ',down_sig_ES)
            print('ES_qber_d = ',down_qber_ES)
            
            np.save('ES_sig_d2.npy', down_sig_ES)
            np.save('ES_qber_d2.npy', down_qber_ES)
            
            
        elif mode == 'up':
            print('ES_sig_u = ', up_sig_ES)
            print('ES_qber_u = ',up_qber_ES)
            
            np.save('ES_sig_u2.npy', up_sig_ES)
            np.save('ES_qber_u2.npy', up_qber_ES)
            
    elif mode2 == 'decoy':
        if mode == 'down':
            print('DEC_sig_d = ', down_sig_decoy)
            print('DEC_qber_d = ', down_qber_decoy)
            
            np.save('DEC_sig_d2.npy', down_sig_decoy)
            np.save('DEC_qber_d2.npy', down_qber_decoy)
            
        elif mode == 'up':
            print('DEC_sig_u = ', up_sig_decoy)
            print('DEC_qber_u = ', up_qber_decoy)
            
            np.save('DEC_sig_u2.npy', up_sig_decoy)
            np.save('DEC_qber_u2.npy', up_qber_decoy)
                
                
        print()
        print()
    # print('DECOY')





################# SET 1 PARAMS
# ES_sig_u =  [3918.2548637623013, 3696.2548637623013, 3589.2548637623013, 2594.2548637623013, 1727.254863762301, 982.254863762301, 381.25486376230106, 93.25486376230106, 32.25486376230106]
# ES_qber_u =  [0.0654260220021295, 0.06601498188477006, 0.06632138667982126, 0.07020035151070908, 0.07666753105605267, 0.0900707444666872, 0.134350717282892, 0.3373032239833963, 0.7755781857378593]

# ES_sig_d =  [4943.906857970288, 4632.906857970288, 4710.906857970288, 3605.9068579702875, 2795.9068579702875, 1698.9068579702875, 896.9068579702877, 240.90685797028763, 18.906857970287632]
# ES_qber_d =  [0.06058756905433605, 0.060951662837563134, 0.060856871727673824, 0.06248323840728706, 0.06426977473680652, 0.06860044867155271, 0.07634189700060764, 0.10656082811294733, 0.40784985024450116]

# DEC_sig_d =  [5366.0, 5036.000000000001, 5175.0, 3948.0000000000005, 3065.0, 1873.0, 1009.9999999999999, 0.0, 0.0]
# DEC_qber_d =  [0.09218297033648559, 0.09367765728033806, 0.09276554839363213, 0.10011001946655669, 0.10801398465350315, 0.1310538830921644, 0.18437110444154636, 0.9316609980052615, -0.2984496005205455]

# DEC_sig_u =  [5806.999999999999, 5572.0, 5475.999999999999, 4236.0, 3131.9999999999995, 0.0, 0.0, 0.0, 0.0]
# DEC_qber_u =  [0.10283918356321711, 0.104892184725007, 0.10617063107961638, 0.12206273974405295, 0.15088252642744424, 0.2302802116840424, -155.39323945641266, -0.08240345664588715, 0.1623691592427473]

# ################# SET 2 PARAMS
# DEC_sig_d2 =  [2707.0000000000005, 2538.0, 2583.0, 1983.9999999999998, 1534.0, 944.9999999999999, 0.0, 0.0, 0.0]
# DEC_qber_d2 =  [0.11140828233392765, 0.11337585207934214, 0.11333953708128304, 0.1244351187335203, 0.13941427193317615, 0.18064036057722901, 0.29302083759086234, -2.835589382867436, -0.36955505060705857]
# #                   90                   80              70                 60                  50                  40                    30                20                 10                    
# ES_sig_u2 =  [1963.6274318811506, 1841.6274318811506, 1793.6274318811506, 1302.6274318811506, 865.6274318811505, 492.6274318811505, 191.62743188115053, 46.62743188115053, 15.627431881150532]
# ES_qber_u2 =  [0.07030022973608, 0.07111838442777486, 0.07146621704595728, 0.07625084672900914, 0.08428765797317345, 0.10062699147479484, 0.15357435570744843, 0.39383647018744905, 0.9161173777770593]

# DEC_sig_u2 =  [2943.0, 2793.0, 2740.0000000000005, 2110.9999999999995, 0.0, 0.0, 0.0, 0.0, 0.0]
# DEC_qber_u2 =  [0.12584079765262052, 0.12986827086649436, 0.13093563385785614, 0.15829242396330015, 0.2072583959935924, 0.3928643096255621, -0.8582193610102633, -0.11748295202094, 0.12668703703789777]

# DEC_sig_d2 =  [2701.0000000000005, 2529.0000000000005, 2588.0, 1980.0, 1533.0, 940.0000000000001, 0.0, 0.0, 0.0]
# DEC_qber_d2 =  [0.11138578510872439, 0.11401003047006124, 0.1125792105015225, 0.12457218065789018, 0.13947453946435812, 0.18205235438378964, 0.2939263468143824, -2.671567735144436, -0.3695736961086098]
##################################################################################################

############# Signal finding error probability 
import scipy.special as spp


def p_err(sig):
    N = 50000
    dt = 2*500e-12
    rb = 400
    pb = ra * rb * (dt)**2
    p = N/2 * spp.erfc(sig/np.sqrt(2*pb*N))
    return p
    
    

###################################################################################################

ES_sig_u = np.load('ES_sig_u.npy')
ES_qber_u = np.load('ES_qber_u.npy')

ES_sig_d = np.load('ES_sig_d.npy')
ES_qber_d = np.load('ES_qber_d.npy')

DEC_sig_u = np.load('DEC_sig_u.npy')
DEC_qber_u = np.load('DEC_qber_u.npy')

DEC_sig_d = np.load('DEC_sig_d.npy')
DEC_qber_d = np.load('DEC_qber_d.npy')


ES_sig_u2 = np.load('ES_sig_u2.npy')
ES_qber_u2 = np.load('ES_qber_u2.npy')

ES_sig_d2 = np.load('ES_sig_d2.npy')
ES_qber_d2 = np.load('ES_qber_d2.npy')


DEC_sig_u2 = np.load('DEC_sig_u2.npy')
DEC_qber_u2 = np.load('DEC_qber_u2.npy')

DEC_sig_d2 = np.load('DEC_sig_d2.npy')
DEC_qber_d2 = np.load('DEC_qber_d2.npy')



print('ES_sig_d = ', ES_sig_d)
print('ES_qber_d = ', ES_qber_d)
print('ES_sig_u = ', ES_sig_u)
print('ES_qber_u = ',ES_qber_u)
print('DEC_sig_d = ', DEC_sig_d)
print('DEC_qber_d = ', DEC_qber_d)
print('DEC_sig_u = ', DEC_sig_u)
print('DEC_qber_u = ', DEC_qber_u)

print('ES_sig_d2 = ', ES_sig_d2 )
print('ES_qber_d2 = ', ES_qber_d2)
print('ES_sig_u2 = ', ES_sig_u2)
print('ES_qber_u2 = ',ES_qber_u2)
print('DEC_sig_d2 = ', DEC_sig_d2)
print('DEC_qber_d2 = ', DEC_qber_d2)
print('DEC_sig_u2 = ', DEC_sig_u2)
print('DEC_qber_u2 = ', DEC_qber_u2)


DEC_sig_u[DEC_sig_u == 0] = -1
DEC_sig_d[DEC_sig_d == 0] = -1
DEC_sig_d2[DEC_sig_d2 == 0] = -1
DEC_sig_u2[DEC_sig_u2 == 0] = -1

jitter = 364e-12
jitter_ES = np.sqrt(2)*jitter


ES_prec_u = 2*jitter_ES/np.sqrt(ES_sig_u)
ES_prec_d = 2*jitter_ES/np.sqrt(ES_sig_d)

ES_prec = 1e12*np.sqrt(ES_prec_u**2 + ES_prec_d**2)/2
ES_qber = np.maximum(ES_qber_u, ES_qber_d)

DEC_prec_u = 2*jitter/np.sqrt(DEC_sig_u)
DEC_prec_d = 2*jitter/np.sqrt(DEC_sig_d)

DEC_prec = 1e12*np.sqrt(DEC_prec_u**2 + DEC_prec_d**2)/2
DEC_qber = np.maximum(DEC_qber_u, DEC_qber_d)



ES_prec_u2 = 2*jitter_ES/np.sqrt(ES_sig_u2)
ES_prec_d2 = 2*jitter_ES/np.sqrt(ES_sig_d2)

ES_prec2 = 1e12*np.sqrt(ES_prec_u2**2 + ES_prec_d2**2)/2
ES_qber2 = np.maximum(ES_qber_u2, ES_qber_d2)

DEC_prec_u2 = 2*jitter/np.sqrt(DEC_sig_u2)
DEC_prec_d2 = 2*jitter/np.sqrt(DEC_sig_d2)

DEC_prec2 = 1e12*np.sqrt(DEC_prec_u2**2 + DEC_prec_d2**2)/2
DEC_qber2 = np.maximum(DEC_qber_u2, DEC_qber_d2)




DEC_qber[DEC_qber < 0] = np.nan
DEC_qber2[DEC_qber2 < 0] = np.nan
DEC_qber[6] = 0.29

ES_prec[ES_qber > 0.25] = np.nan
ES_prec2[ES_qber2 > 0.25] = np.nan


# warnings.filterwarnings('ignore')

plt.close('all')

plt.rcParams["font.family"] = "Times New Roman"


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)



ax.plot(90-np.array(zs), ES_prec, 'o-', label='ES-1', linewidth=3, markersize=5)
# ax.plot(zs[:6], ES_sig_d, 'o', label='ES DOWNLINK', linewidth=3)
ax.plot(90-np.array(zs), ES_prec2, 'o-', label='ES-2', linewidth=3, markersize=5)



ax.plot(90-np.array(zs), DEC_prec, '^-', label='Decoy-1', linewidth=3, markersize=6)
# # ax.plot(zs[:6], DEC_sig_d, 'o', label='DECOY DOWNLINK', linewidth=3)
ax.plot(90-np.array(zs), DEC_prec2, '^-', label='Decoy-2', linewidth=3, markersize=6)


# # ax.plot(dists[:], ss400_l, color=cs[1])
# # ax.plot(dists[:], ss400_u, color=cs[1])
# # ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# # ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# # ax.set_yscale('log')
ax.set_xlabel('Elevation (deg)', fontsize=15)
ax.set_ylabel('Precision (picosec)', fontsize=15)

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

# # ax.set_ylim([0.04, 0.25])
# # ax.set_xlim([6500, 30000])

plt.tight_layout()

fig.savefig('sig_sat.pdf', dpi=200)



fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)



ax.plot(90-np.array(zs)[:8], ES_qber[:8], 'o-', label='ES-1', linewidth=3, markersize=5)
# # ax.plot(zs[:6], ES_sig_d, 'o', label='ES DOWNLINK', linewidth=3)
ax.plot(90-np.array(zs)[:8], ES_qber2[:8], 'o-', label='ES-2', linewidth=3, markersize=5)


ax.plot(90-np.array(zs)[:5], DEC_qber[:5], '^-', label='Decoy-1', linewidth=3, markersize=6)
# # ax.plot(zs[:6], DEC_sig_d, 'o', label='DECOY DOWNLINK', linewidth=3)
ax.plot(90-np.array(zs)[:4], DEC_qber2[:4], '^-', label='Decoy-2', linewidth=3, markersize=6)

ax.set_ylim([0, .25])
# # ax.plot(dists[:], ss400_l, color=cs[1])
# # ax.plot(dists[:], ss400_u, color=cs[1])
# # ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


# # ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

# # ax.set_yscale('log')
ax.set_xlabel('Elevation (deg)', fontsize=15)
ax.set_ylabel('QBER', fontsize=15)

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

# ax.set_ylim([0.04, 0.25])
# ax.set_xlim([6500, 30000])

plt.tight_layout()

fig.savefig('qbers_sat.pdf', dpi=200)

