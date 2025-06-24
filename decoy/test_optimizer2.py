# -*- coding: utf-8 -*-


# -*- coding: utf-8 -*-

import optimizer
import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

ra = 20e6 # pulse rate
bg_rate = 500e3 # day background rate
# bg_rate = 10e3 # night background rate

jitter = 364e-12 # total time jitter
e_d = 0.05 # detector QBER
time_filter = 2 # 2 sigma time filter on the coincidence window

channel_loss = 11.5 # in dB
detector_loss = 5 
eta_db = channel_loss + detector_loss # total loss
eta = 10**(-eta_db/10) # loss
time = 1 # total time of prototocol

eta = eta*0.6826
eta_es = 10**(- (channel_loss + 2*detector_loss)/10)*0.6826



channel = {
    'ra': ra,
    'channel_loss': channel_loss,
    'detector_loss': detector_loss,
    'bg_rate': bg_rate,
    'jitter': jitter,
    'time_filter': time_filter,
    'time': time, 
    'e_d': e_d}
opt = optimizer.Optimizer(channel)

m = .6 # signal mean photon number
v = 0.3 # decoy mean photon number  - required v < m


pm = 0.4 # probability of sigal produced

x = [m, v, pm]

opt.set_initial_conditions(x)
opt.set_confidence(0.99)
# opt.set_confidence(0.5)


opt.optimize2()

print('Num bg photons: ', opt.bg_prob*ra*time)

print("x: ", opt.x)

print(opt.theoretical_vals(opt.x, opt.conf))
# print(opt.theoretical_vals(opt.x, 0.5))


m, v, pm = x

# print('Signal: ', ra*time*eta*(pm*(1-np.exp(-m)) + (1-pm)*(1-np.exp(-v))))
print('Signal: ', ra*time*(pm*(opt.bg_prob + 1 - np.exp(-eta*m)) + (1-pm)*(opt.bg_prob + 1 - np.exp(-eta*v))))
print('Signal/2: ', .5*ra*time*(pm*(opt.bg_prob + 1 - np.exp(-eta*m)) + (1-pm)*(opt.bg_prob + 1 - np.exp(-eta*v))))



print('------------- REF values')
mr = 0.55
vr = 0.3
pmr = 0.6

xr = [mr, vr, pmr]
print(xr)
print(opt.theoretical_vals(xr, opt.conf))
# print(opt.theoretical_vals(xr, .5))
# print(opt.theoretical_vals2(confidence))

print('Entangled source coincidences: ', ra*time*(eta_es + opt.bg_prob))
# print('Signal ref: ', ra*time*eta*(pmr*(1-np.exp(-mr)) + (1-pmr)*(1-np.exp(-vr))))
print('Signal ref: ', ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))
print('Signal ref/2: ', .5*ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))


print('------------- REF values')
mr = 0.55
vr = 0.3
pmr = 0.7

xr = [mr, vr, pmr]
print(xr)
print(opt.theoretical_vals(xr, opt.conf))
# print(opt.theoretical_vals(xr, .5))
# print(opt.theoretical_vals2(confidence))

print('Entangled source coincidences: ', ra*time*(eta_es + opt.bg_prob))
# print('Signal ref: ', ra*time*eta*(pmr*(1-np.exp(-mr)) + (1-pmr)*(1-np.exp(-vr))))
print('Signal ref: ', ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))
print('Signal ref/2: ', .5*ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))

print('------------- REF values')
mr = 0.55
vr = 0.3
pmr = 0.8

xr = [mr, vr, pmr]
print(xr)
print(opt.theoretical_vals(xr, opt.conf))
# print(opt.theoretical_vals(xr, .5))
# print(opt.theoretical_vals2(confidence))

print('Entangled source coincidences: ', ra*time*(eta_es + opt.bg_prob))
# print('Signal ref: ', ra*time*eta*(pmr*(1-np.exp(-mr)) + (1-pmr)*(1-np.exp(-vr))))
print('Signal ref: ', ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))
print('Signal ref/2: ', .5*ra*time*(pmr*(opt.bg_prob + 1 - np.exp(-eta*mr)) + (1-pmr)*(opt.bg_prob + 1 - np.exp(-eta*vr))))





########################################################################################
########################################################################################


conf = 0.99
vs = [0.15, 0.3, 0.47, 0.5]
pm = 0.6

ms1 = np.linspace(vs[0]+0.01, 0.99,200)
ms2 = np.linspace(vs[1]+0.01, 0.99,200)
ms3 = np.linspace(vs[2]+0.01, 0.99,200)
ms4 = np.linspace(vs[3]+0.01, 0.99,200)

ns1 = []
ns2 = []
ns3 = []
ns4 = []

qs1 = []
qs2 = []
qs3 = []
qs4 = []


for i in range(len(ms1)):
    x1 = [ms1[i], vs[0], pm] 
    n1, q1 = opt.theoretical_vals(x1, opt.conf)
    ns1 += [n1]
    qs1 += [q1]

    x2 = [ms2[i], vs[1], pm] 
    n2, q2 = opt.theoretical_vals(x2, opt.conf)
    ns2 += [n2]
    qs2 += [q2]

    x3 = [ms3[i], vs[2], pm] 
    n3, q3 = opt.theoretical_vals(x3, opt.conf)
    ns3 += [n3]
    qs3 += [q3]
    
    x4 = [ms4[i], vs[3], pm] 
    n4, q4 = opt.theoretical_vals(x4, opt.conf)
    ns4 += [n4]
    qs4 += [q4]
    

plt.close('all')




ms0 = np.linspace(0, 1, 200)
# Security threshold
prob_ph = 1- np.exp(-ms0)
# sec_lim1 = np.ones(200)*(ra*time*eta*pm*prob_ph)/2

plt.figure()
plt.plot(ms1, ns1, label=vs[0])
plt.plot(ms2, ns2, label=vs[1])
plt.plot(ms3, ns3, label=vs[2])
plt.plot(ms4, ns4, label=vs[3])

# plt.plot(ms0, sec_lim, 'r--')


sig1 = np.ones(200)*ra*time*(opt.bg_prob + pm*(1 -  np.exp(-eta*ms1)) + (1-pm)*(1 - np.exp(-eta*vs[0])))
sig2 = np.ones(200)*ra*time*(opt.bg_prob + pm*(1 -  np.exp(-eta*ms2)) + (1-pm)*(1 - np.exp(-eta*vs[1])))
sig3 = np.ones(200)*ra*time*(opt.bg_prob + pm*(1 -  np.exp(-eta*ms3)) + (1-pm)*(1 - np.exp(-eta*vs[2])))
sig4 = np.ones(200)*ra*time*(opt.bg_prob + pm*(1 -  np.exp(-eta*ms4)) + (1-pm)*(1 - np.exp(-eta*vs[3])))



plt.gca().set_prop_cycle(None)

plt.plot(ms1, sig1/2, '--')
plt.plot(ms2, sig2/2, '--')
plt.plot(ms3, sig3/2, '--')
plt.plot(ms4, sig4/2, '--')

plt.legend()
# # plt.ylim([0, 5e6])
# plt.xlim([0, 1])

# plt.figure()
# plt.plot(ms1, ns1-sig1/2, label=vs[0])
# plt.plot(ms2, ns2-sig2/2, label=vs[1])
# plt.plot(ms3, ns3-sig3/2, label=vs[2])
# plt.plot(ms4, ns4-sig4/2, label=vs[3])

# # plt.plot(ms0, sec_lim, 'r--')
# plt.legend()
# plt.xlim([0, 1])

plt.figure()
plt.plot(ms1, qs1, label=vs[0])
plt.plot(ms2, qs2, label=vs[1])
plt.plot(ms3, qs3, label=vs[2])
plt.plot(ms4, qs4, label=vs[3])
plt.legend()
plt.ylim([-0.001, 0.25])
plt.xlim([0, 1])


###############################################################################

conf = 0.99
vs = 0.3
pm = [0.3, 0.5, 0.6, 0.7]

ms1 = np.linspace(vs + 0.01, 0.99,200)
ms2 = np.linspace(vs + 0.01, 0.99,200)
ms3 = np.linspace(vs + 0.01, 0.99,200)
ms4 = np.linspace(vs + 0.01, 0.99,200)

ns1 = []
ns2 = []
ns3 = []
ns4 = []

qs1 = []
qs2 = []
qs3 = []
qs4 = []


for i in range(len(ms1)):
    x1 = [ms1[i], vs, pm[0]] 
    n1, q1 = opt.theoretical_vals(x1, opt.conf)
    ns1 += [n1]
    qs1 += [q1]

    x2 = [ms2[i], vs, pm[1]] 
    n2, q2 = opt.theoretical_vals(x2, opt.conf)
    ns2 += [n2]
    qs2 += [q2]

    x3 = [ms3[i], vs, pm[2]] 
    n3, q3 = opt.theoretical_vals(x3, opt.conf)
    ns3 += [n3]
    qs3 += [q3]
    
    x4 = [ms4[i], vs, pm[3]] 
    n4, q4 = opt.theoretical_vals(x4, opt.conf)
    ns4 += [n4]
    qs4 += [q4]
    

# plt.close('all')




ms0 = np.linspace(0, 1, 200)
# Security threshold
prob_ph = 1- np.exp(-ms0)
# sec_lim1 = np.ones(200)*(ra*time*eta*pm*prob_ph)/2

plt.figure()
plt.plot(ms1, ns1, label=pm[0])
plt.plot(ms2, ns2, label=pm[1])
plt.plot(ms3, ns3, label=pm[2])
plt.plot(ms4, ns4, label=pm[3])

# plt.plot(ms0, sec_lim, 'r--')


sig1 = np.ones(200)*ra*time*(opt.bg_prob + pm[0]*(1 -  np.exp(-eta*ms1)) + (1-pm[0])*(1 - np.exp(-eta*vs)))
sig2 = np.ones(200)*ra*time*(opt.bg_prob + pm[1]*(1 -  np.exp(-eta*ms2)) + (1-pm[1])*(1 - np.exp(-eta*vs)))
sig3 = np.ones(200)*ra*time*(opt.bg_prob + pm[2]*(1 -  np.exp(-eta*ms3)) + (1-pm[2])*(1 - np.exp(-eta*vs)))
sig4 = np.ones(200)*ra*time*(opt.bg_prob + pm[3]*(1 -  np.exp(-eta*ms4)) + (1-pm[3])*(1 - np.exp(-eta*vs)))



plt.gca().set_prop_cycle(None)

plt.plot(ms1, sig1/2, '--')
plt.plot(ms2, sig2/2, '--')
plt.plot(ms3, sig3/2, '--')
plt.plot(ms4, sig4/2, '--')




# plt.legend()
# # plt.ylim([0, 5e6])
# plt.xlim([0, 1])

# plt.figure()
# plt.plot(ms1, ns1-sig1/2, label=vs[0])
# plt.plot(ms2, ns2-sig2/2, label=vs[1])
# plt.plot(ms3, ns3-sig3/2, label=vs[2])
# plt.plot(ms4, ns4-sig4/2, label=vs[3])

# # plt.plot(ms0, sec_lim, 'r--')
# plt.legend()
# plt.xlim([0, 1])

plt.figure()
plt.plot(ms1, qs1, label=pm[0])
plt.plot(ms2, qs2, label=pm[1])
plt.plot(ms3, qs3, label=pm[2])
plt.plot(ms4, qs4, label=pm[3])
plt.legend()
plt.ylim([-0.001, 0.25])
plt.xlim([0, 1])





