# -*- coding: utf-8 -*-


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

els = sp.io.loadmat('data/els.mat')
times = sp.io.loadmat('data/times.mat')


els = els['el_arr'][0]
times = times['secs'][0]

plt.plot(times, els, 'bo')


times = times[330:-1]
times = times - times[0]
els = els[330:-1]








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
    file_up = f'data/PAPER_UPLINK_z={z}_wl=0.81_w0=0.35_ap=0.2_H=500_perr=1.2e-06'
    data = sp.io.loadmat(file_up)
    data = data['res'].transpose()[0]

    data = data * np.exp(-0.7/np.cos(np.deg2rad(z)))

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
    file_down = f'data/PAPER_DOWNLINK_z={z}_wl=0.81_w0=0.2_ap=0.35_H=500_perr=1.2e-06'
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


up_loss_std[-4] += 0.6 


fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.errorbar(zs, up_loss, up_loss_std, linewidth=2,elinewidth=0, capsize=5, label='Uplink')
plt.errorbar(zs, down_loss, down_loss_std, linewidth=2,elinewidth=0, capsize=5, label='Downlink')

# plt.plot(zs, up_loss, linewidth=3, label='Uplink')
# plt.plot(zs, down_loss, linewidth=3, label='Downlink')


plt.xlabel('Elevation (deg)',fontsize=13)
plt.ylabel('Loss (dB)', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)


plt.legend()

plt.tight_layout()
fig.savefig('satellite_loss.pdf', dpi=200)



from scipy.interpolate import interp1d

x = zs
y = down_loss

# Linear interpolation
linear_interp_down = interp1d(x, y, kind='linear')

# Cubic interpolation
cubic_interp_down = interp1d(x, y, kind='cubic')


x_new = np.linspace(x.min(), x.max(), 100)  # 100 evenly spaced points between min and max of x

y_linear = linear_interp_down(x_new)
y_cubic = cubic_interp_down(x_new)


plt.figure(figsize=(8, 5))
plt.plot(x, y, 'o', label='Original data points', markersize=8)
plt.plot(x_new, y_linear, '-', label='Linear interpolation')
plt.plot(x_new, y_cubic, '--', label='Cubic interpolation ')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.title('DOWNLINK')
plt.show()

x = zs
y = up_loss

# Linear interpolation
linear_interp_up = interp1d(x, y, kind='linear')

# Cubic interpolation
cubic_interp_up = interp1d(x, y, kind='cubic')


x_new = np.linspace(x.min(), x.max(), 100)  # 100 evenly spaced points between min and max of x

y_linear = linear_interp_up(x_new)
y_cubic = cubic_interp_up(x_new)


plt.figure(figsize=(8, 5))
plt.plot(x, y, 'o', label='Original data points', markersize=8)
plt.plot(x_new, y_linear, '-', label='Linear interpolation')
plt.plot(x_new, y_cubic, '--', label='Cubic interpolation')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.title('UPLINK')
plt.show()




pass_loss_up = linear_interp_up(els)
pass_loss_down = linear_interp_down(els)

plt.figure()
plt.plot(times, pass_loss_up)
plt.plot(times, pass_loss_down)


np.save('pass_loss_down', pass_loss_down)
np.save('pass_loss_up', pass_loss_up)



x = zs
y = down_loss_std

# Linear interpolation
linear_interp_down = interp1d(x, y, kind='linear')

# Cubic interpolation
cubic_interp_down = interp1d(x, y, kind='cubic')


x_new = np.linspace(x.min(), x.max(), 100)  # 100 evenly spaced points between min and max of x

y_linear = linear_interp_down(x_new)
y_cubic = cubic_interp_down(x_new)


plt.figure(figsize=(8, 5))
plt.plot(x, y, 'o', label='Original data points', markersize=8)
plt.plot(x_new, y_linear, '-', label='Linear interpolation')
plt.plot(x_new, y_cubic, '--', label='Cubic interpolation ')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.title('DOWNLINK - STD')
plt.show()


x = zs
y = up_loss_std

# Linear interpolation
linear_interp_up = interp1d(x, y, kind='linear')

# Cubic interpolation
cubic_interp_up = interp1d(x, y, kind='cubic')


x_new = np.linspace(x.min(), x.max(), 100)  # 100 evenly spaced points between min and max of x

y_linear = linear_interp_up(x_new)
y_cubic = cubic_interp_up(x_new)


plt.figure(figsize=(8, 5))
plt.plot(x, y, 'o', label='Original data points', markersize=8)
plt.plot(x_new, y_linear, '-', label='Linear interpolation')
plt.plot(x_new, y_cubic, '--', label='Cubic interpolation ')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.title('UPLINK - STD')
plt.show()


pass_loss_up_std = linear_interp_up(els)
pass_loss_down_std = linear_interp_down(els)

plt.figure()
plt.plot(times, pass_loss_up_std)
plt.plot(times, pass_loss_down_std)


np.save('pass_loss_down_std', pass_loss_down_std)
np.save('pass_loss_up_std', pass_loss_up_std)
