# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import satellite_link_es as linkes
import satellite_link_decoy as linkdec
import matplotlib.pyplot as plt

plt.close('all')


ra = 20e6
# mode = 'down'
mode = 'up'
mode2 = 'es'
# mode2 = 'decoy'


if mode2 == 'es':
    link = linkes.SatelliteLinkES(ra)
else:
    link = linkdec.SatelliteLinkDecoy(ra)
    
if mode == 'down':
    loss_avg = np.load('pass_loss_down.npy')
    loss_std = np.load('pass_loss_down_std.npy')
if mode == 'up':
    loss_avg = np.load('pass_loss_up.npy')
    loss_std = np.load('pass_loss_up_std.npy')


plt.figure()
plt.title('loss ' + mode)
plt.plot(loss_avg)


l = len(loss_avg)

sig = []
qber = []

qm = []
qv = []
y1 = []
y0 = []
opt_params = []

if 1: 
    for i in range(l):
        lavg = loss_avg[i]
        lstd = loss_std[i]
        
        print('loss avg = ', lavg + 10)
        
        res = link.run_pnt(lavg, lstd)
        
        
        if mode2 == 'es':
            sig += [res[0]]
            qber += [res[1]]
    
        if mode2 == 'decoy':
            opt_p = link.get_opt_params()
            qm += [res[0]]
            qv += [res[1]]
            y1 += [res[2]]
            y0 += [res[3]]
            qber += [res[4][2]]
            opt_params += [list(opt_p)]
       
        # print(signal, qber)
        # print(i)
    
    
# sig = np.array(sig)
# qber = np.array(qber)    

    opt_params = list(opt_params)    
    
    if mode == 'down':
        if mode2 == 'es':
            # print('pass_sig_down = ', sig)
            # print('pass_qber_down = ', qber)
            np.save('pass_sig_down2.npy', sig)
            np.save('pass_qber_down2.npy', qber)
        
        if mode2 == 'decoy':
            # print('pass_qber_decoy_down = ', qber)
            # print('pass_qm_decoy_down = ', qm)
            # print('pass_qv_decoy_down = ', qv)
            # print('pass_y1_decoy_down = ', y1)
            # print('pass_y0_decoy_down = ', y0)
            # print('pass_opt_decoy_down = ', opt_params)
            
            opt_params = np.array(opt_params)
            qm = np.array(qm)
            qv = np.array(qv)
            y1 = np.array(y1)
            y0 = np.array(y0)
            
            m = opt_params[:, 0]
            v = opt_params[:, 1]
            pm = opt_params[:, 2]
            bm = (qm * np.exp(m) - y0)/(2*m) < y1
            bv = (qv * np.exp(v) - y0)/(2*v) < y1
            nm = ra*qm*pm
            nv = ra*qv*(1-pm)
            sig = bm*nm + bv*nv
            
            np.save('pass_sig_down_decoy2.npy', sig)
            np.save('pass_qber_down_decoy2.npy', qber)
            print('pass_sig_down_decoy2.npy ------- pass_qber_down_decoy2.npy')
        
    
    if mode == 'up':
        if mode2 == 'es':
            # print('pass_sig_up = ', sig)
            # print('pass_qber_up = ', qber)
            np.save('pass_sig_up2.npy', sig)
            np.save('pass_qber_up2.npy', qber)
    
        if mode2 == 'decoy':
            # print('pass_qber_decoy_up = ', qber)
            # print('pass_qm_decoy_up = ', qm)
            # print('pass_qv_decoy_up = ', qv)
            # print('pass_y1_decoy_up = ', y1)
            # print('pass_y0_decoy_up = ', y0)
            # print('pass_opt_decoy_up = ', opt_params)
            
            opt_params = np.array(opt_params)
            qm = np.array(qm)
            qv = np.array(qv)
            y1 = np.array(y1)
            y0 = np.array(y0)
            
            m = opt_params[:, 0]
            v = opt_params[:, 1]
            pm = opt_params[:, 2]
            bm = (qm * np.exp(m) - y0)/(2*m) < y1
            bv = (qv * np.exp(v) - y0)/(2*v) < y1
            nm = ra*qm*pm
            nv = ra*qv*(1-pm)
            sig = bm*nm + bv*nv
            
            np.save('pass_sig_up_decoy2.npy', sig)
            np.save('pass_qber_up_decoy2.npy', qber)



sig_down = np.load('pass_sig_down.npy')
qber_down = np.load('pass_qber_down.npy')

sig_up = np.load('pass_sig_up.npy')
qber_up = np.load('pass_qber_up.npy')

sig_decoy_down = np.load('pass_sig_down_decoy.npy')
sig_decoy_up = np.load('pass_sig_up_decoy.npy')

qber_decoy_down = np.load('pass_qber_down_decoy.npy')
qber_decoy_up = np.load('pass_qber_up_decoy.npy')

sig_down2 = np.load('pass_sig_down2.npy')
qber_down2 = np.load('pass_qber_down2.npy')

sig_up2 = np.load('pass_sig_up2.npy')
qber_up2 = np.load('pass_qber_up2.npy')

sig_decoy_down2 = np.load('pass_sig_down_decoy2.npy')
sig_decoy_up2 = np.load('pass_sig_up_decoy2.npy')

qber_decoy_down2 = np.load('pass_qber_down_decoy2.npy')
qber_decoy_up2 = np.load('pass_qber_up_decoy2.npy')


jitter = 400e-12
jitter_es = np.sqrt(2)*jitter


# ES_prec_u = 2*jitter_ES/np.sqrt(ES_sig_u)
# ES_prec_d = 2*jitter_ES/np.sqrt(ES_sig_d)

t = np.arange(len(sig_down))

sig_es = (sig_down + sig_up)/2
sig_es2 = (sig_down2 + sig_decoy_up2)/2

sig_decoy = (sig_decoy_down + sig_decoy_up)/2
sig_decoy2 = (sig_decoy_down2 + sig_decoy_up2)/2

prec_es_d = jitter_es/np.sqrt(sig_down)
prec_es_u = jitter_es/np.sqrt(sig_up)

prec_es_d2 = jitter_es/np.sqrt(sig_down2)
prec_es_u2 = jitter_es/np.sqrt(sig_up2)

prec_decoy_d = jitter/np.sqrt(sig_decoy_down)
prec_decoy_u = jitter/np.sqrt(sig_decoy_up)

prec_decoy_d2 = jitter/np.sqrt(sig_decoy_down2)
prec_decoy_u2 = jitter/np.sqrt(sig_decoy_up2)


prec_es = np.sqrt(prec_es_d**2 + prec_es_u**2)/2
prec_es2 = np.sqrt(prec_es_d2**2 + prec_es_u2**2)/2

prec_decoy = np.sqrt(prec_decoy_d**2 + prec_decoy_u**2)/2
prec_decoy2 = np.sqrt(prec_decoy_d2**2 + prec_decoy_u2**2)/2



qber_es = np.maximum(qber_down, qber_up)
qber_es2 = np.maximum(qber_down2, qber_up2)

qber_decoy = np.maximum(qber_decoy_down, qber_decoy_up)
qber_decoy2 = np.maximum(qber_decoy_down2, qber_decoy_up2)




prec_decoy[qber_decoy >= 0.25] = np.inf
prec_decoy2[qber_decoy2 >= 0.25] = np.inf

prec_es[qber_es >= 0.25] = np.inf
prec_es2[qber_es2 >= 0.25] = np.inf


qber_decoy[prec_decoy == np.inf] = np.nan
qber_decoy2[prec_decoy2 == np.inf] = np.nan


print('Protocol fails = ', sum(np.isnan(qber_decoy))/len(qber_decoy))
print('Protocol fails = ', sum(np.isnan(qber_decoy2))/len(qber_decoy))


# plt.figure()
# plt.title('Signal')
# plt.plot(t, sig_es, label='Signal ES')
# plt.plot(t, sig_es2, label='Signal ES 2')

# plt.plot(t, sig_decoy, label='Signal decoy')
# plt.plot(t, sig_decoy2, label='Signal decoy 2')

# plt.yscale('log')


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)
# ax.set_title('Precision')
ax.plot(t+293, prec_es*1e12, '.', label='ES-1')
ax.plot(t+293, prec_es2*1e12, '.', label='ES-2')


# plt.figure()
ax.plot(t+293, prec_decoy*1e12, '^', label='Decoy-1', markersize=4)
ax.plot(t+293, prec_decoy2*1e12, '^', label='Decoy-2', markersize=4)

ax.legend()
ax.set_xlabel('Orbit time (s)', fontsize=15)
ax.set_ylabel('Precision (ps)', fontsize=15)
plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=15)

plt.tight_layout()

fig.savefig('pass_precision.pdf', dpi=200)



fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)
# ax.set_title('QBER')
ax.plot(t+293, qber_es, '.', label='ES-1')
ax.plot(t+293, qber_es2, '.', label='ES-2')


ax.plot(t+293, qber_decoy, '^', label='Decoy-1', markersize=4)
ax.plot(t+293, qber_decoy2, '^', label='Decoy-2', markersize=4)


ax.set_ylim([0, .25])

# ax.legend()
ax.set_xlabel('Orbit time (s)', fontsize=15)
ax.set_ylabel('QBER', fontsize=15)
# plt.yscale('log')

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
# plt.legend(fontsize=15)

plt.tight_layout()

fig.savefig('pass_qber.pdf', dpi=200)


# Nmd = 20e6*pass_qm_decoy_down *
# pass_ signal_decoy_down = 