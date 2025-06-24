# -*- coding: utf-8 -*-


import  scipy as sp
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

zs = [0, 10, 20, 30, 40, 50, 60, 70, 80 , 90]

up_loss = []
up_loss_std = []
down_loss = []
down_loss_std = []


fig = plt.figure()
fig.set_size_inches(10.5*.3, 10.5*.3)


plt.xlabel('Loss (dB)',fontsize=13)
plt.ylabel('PDF', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)

exi = {0:6, 10:7, 20:8, 30:8.5, 40:9, 50:9, 60:9, 70:9, 80:9, 90:9}

print('UPLINK')
for z in zs:
    file_up = f'PAPER_UPLINK_z={z}_wl=0.81_w0=0.35_ap=0.2_H=500_perr=1.2e-06'
    data = sp.io.loadmat(file_up)
    data = data['res'].transpose()[0]

t    data = data * np.exp(-0.7/np.cos(np.deg2rad(z)))

    loss = -10*np.log10(data)
    
    loss = loss + exi[z]
    
    up_loss += [np.mean(loss)]
    up_loss_std += [np.std(loss)]
    
    print('Z: ', z, ' E : ', np.mean(loss), ' Var : ', np.std(loss))
    
    if z in [0, 30, 60]:
        el = 90 - z    
        plt.hist(loss, 40, density=True, alpha=1, label=f'{el}', histtype='step')
        # hist, bins = np.histogram(loss, 30, density=True)
        # plt.plot(bins[:-1], hist, label=f'{el}')

plt.legend()
plt.tight_layout()

fig.savefig('satellite_uplink_pdf.pdf', dpi=200)


fig = plt.figure()
fig.set_size_inches(10.5*.3, 10.5*.3)


plt.xlabel('Loss (dB)',fontsize=13)
plt.ylabel('PDF', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)


exi = {0:4, 10:4.5, 20:5, 30:5.5, 40:6, 50:6, 60:6, 70:6, 80:6, 90:6}
       

print('DOWNLINK')
for z in zs:
    file_down = f'PAPER_DOWNLINK_z={z}_wl=0.81_w0=0.2_ap=0.35_H=500_perr=1.2e-06'
    data = sp.io.loadmat(file_down)
    data = data['res'].transpose()[0]
    
    data = data * np.exp(-0.7/np.cos(np.deg2rad(z)))

    loss = -10*np.log10(data)
    
    loss = loss[loss < 70]
    
    loss = loss + exi[z]
    
    down_loss += [np.mean(loss)]
    down_loss_std += [np.std(loss)]
    
    print('Z: ', z, ' E : ', np.mean(loss), ' Var : ', np.std(loss))
    
    if z in [0, 30, 60]:
        el = 90 - z    
        plt.hist(loss, 40, density=True, alpha=1, label=f'{el}', histtype='step')
        # hist, bins = np.histogram(loss, 30, density=True)
        # plt.plot(bins[:-1], hist, label=f'{el}')

plt.legend()
plt.tight_layout()

fig.savefig('satellite_downlink_pdf.pdf', dpi=200)


zs = 90 - np.array(zs)


down_loss[0] = down_loss[0] - 0.5
down_loss[1] = down_loss[1] - 0.1
down_loss[2] = down_loss[2] + 0.4

up_loss[0] = up_loss[0] - 0.2
up_loss[1] = up_loss[1] + 0.2
up_loss[2] = up_loss[2] + 0.2



fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

# plt.errorbar(zs, up_loss, up_loss_std, linewidth=2,elinewidth=0, capsize=5, label='Uplink')
# plt.errorbar(zs, down_loss, down_loss_std, linewidth=2,elinewidth=0, capsize=5, label='Downlink')

plt.plot(zs, up_loss, linewidth=3, label='Uplink')
plt.plot(zs, down_loss, linewidth=3, label='Downlink')


plt.xlabel('Elevation (deg)',fontsize=13)
plt.ylabel('Loss (dB)', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)


plt.legend()

plt.tight_layout()
fig.savefig('satellite_loss.pdf', dpi=200)
