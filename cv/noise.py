# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 10:21:13 2023

@author: vil034
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from os.path import dirname, join

plt.close('all')

distances = [2000, 4000, 6000, 8000, 10000]
cns = ['4.67e-15', '5e-14', '2.47e-13']

data_dir = '../../Laser_propagation/Horizontal/'


Tmeans_sqrt = np.zeros((len(distances), len(cns)))
Tmeans = np.zeros((len(distances), len(cns)))

           
               
for i in range(len(distances)):
    for j in range(len(cns)):
        mat_name = join(data_dir, 'HORIZONTAL_CV_d=%d_cn=%s' % (distances[i],cns[j]))
        print(mat_name)
        
        mat = sio.loadmat(mat_name)
        data = mat['res']
        
        
        
        Tmeans_sqrt[i, j] = np.average(np.sqrt(data))
        Tmeans[i, j] = np.average(data)
        
        
Tf = Tmeans_sqrt**2
VarT = Tmeans - Tf

Vmod = 20

e_line = 0.1
e_det = 0.
e0 = e_line + e_det/Tf

e = VarT/Tf * Vmod + Tmeans*e0

# e = Tf*e

plt.figure()
plt.plot(distances, e[:,0], label=cns[0])
plt.plot(distances, e[:,1], label=cns[1])
plt.plot(distances, e[:,2], label=cns[2])
plt.legend()
plt.xlabel('Distance')
plt.ylabel(r'$\epsilon$')


plt.figure()
plt.plot(distances, e0[:,0], label=cns[0])
plt.plot(distances, e0[:,1], label=cns[1])
plt.plot(distances, e0[:,2], label=cns[2])
plt.legend()
plt.xlabel('Distance')
plt.ylabel(r'$\epsilon$')

############################## INDIVIDUAL DATA
d = distances[-1]
c = cns[0]

mat_name = join(data_dir, 'HORIZONTAL_CV_d=%s_cn=%s' % (d,c))
# print(mat_name)

mat = sio.loadmat(mat_name)
# print(mat)

data = mat['res']

plt.figure()
plt.hist(data, bins=200)