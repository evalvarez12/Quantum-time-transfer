# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt

#### functions

def h(x):
    return -x*np.log2(x) - (1-x)*np.log2(1-x)


def Delta(n, e, e_pa):
    return 7 * np.sqrt(np.log2(2/e)/n) + 2*np.log2(1/e_pa)/n
    


#### Parameters
ra = 100e6 # source rate
loss = 50 # total loss dB
t = 10**(-loss/10)
time =  7*60 # protocol time 

pz = .5 # probability of choosing  X basis for Alice and Bob
px = 1 - pz

e = 1e-5
e_ec = 1e-10
Q = 0.05
eta = 0.1
p_d = 1e-5

e_bar = (e - e_ec)/3 
e_pe = (e - e_ec)/3
e_pa = (e - e_ec)/3





def key_rate(ra, loss) : 
    
    t = 10**(-loss/10)
    
    N = time*t*ra # number of signals received by Bob


    n = N*pz**2
    m = N*px**2

    e_xu = Q + np.sqrt((2*np.log(1/e_pe) + 2*np.log(m+1))/(2*m))

    Dn = Delta(n, e, e_pa)

    leak = 1.05*h(e)

    K = ra*t*pz**2*(1 - h(e_xu) - Dn - leak)
    return K


def key_rate_INGAS(ra, loss) : 
    
    # 20 % efficiency of the detectors
    t = 10**(-loss/10)*.2
    
    N = time*t*ra # number of signals received by Bob


    n = N*pz**2
    m = N*px**2

    e_xu = Q + np.sqrt((2*np.log(1/e_pe) + 2*np.log(m+1))/(2*m))

    Dn = Delta(n, e, e_pa)

    leak = 1.05*h(e)

    K = ra*t*pz**2*(1 - h(e_xu) - Dn - leak)
    return K


K = key_rate(ra, loss)
print('Key rate (bits/sec): ', K)



plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

# loss = np.array([35, 40, 45, 50, 55, 60, 65])
loss = np.linspace(35, 67, 200)

ras = [40e6, 70e6, 100e6, 130e6]

fig = plt.figure()



K1 = key_rate(ras[0], loss)
K2 = key_rate(ras[1], loss)
K3 = key_rate(ras[2], loss)
K4 = key_rate(ras[3], loss)

K1i = key_rate_INGAS(ras[0], loss)
K2i = key_rate_INGAS(ras[1], loss)
K3i = key_rate_INGAS(ras[2], loss)
K4i = key_rate_INGAS(ras[3], loss)


plt.plot(loss, K1, label='40 MHz - SNSPDS')
plt.plot(loss, K2, label='70 MHz - SNSPDS')
plt.plot(loss, K3, label='100 MHz - SNSPDS')
plt.plot(loss, K4, label='130 MHz - SNSPDS')

plt.gca().set_prop_cycle(None)

plt.plot(loss, K1i,'--', label='40 MHz - InGaAs')
plt.plot(loss, K2i,'--', label='70 MHz - InGaAs')
plt.plot(loss, K3i,'--', label='100 MHz - InGaAs')
plt.plot(loss, K4i,'--', label='130 MHz - InGaAs')


plt.ylim([.5, 10000])
plt.yscale('log')

plt.xlabel('Loss (dB)', fontsize=13)
plt.ylabel('Key rate (bits/sec)', fontsize=13)

plt.xticks(fontsize=13, rotation=0)
plt.yticks(fontsize=13, rotation=0)
plt.legend(fontsize=13)
fig.set_size_inches(18.5*.42, 10.5*.42)
fig.savefig('key_rate.pdf', dpi=200)

