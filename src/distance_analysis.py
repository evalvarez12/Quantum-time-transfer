 # -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 12:18:46 2023

@author: vil034
"""

import numpy as np


def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def get_distances_A_to_B(ta, tb, qa, qb, n, offset_rough, limit = 100):
    # dists = np.zeros(len(ta)*2*n)
    dists = []
    
    # index of shorlisted times in A
    idx_b = []
    
    for i in range(len(ta)):
        t = ta[i]
        idx = find_nearest_idx(tb, t - offset_rough)
        
        a = max(0,idx - n) 
        b = min(idx+n+1,len(tb))
        
        tt = t - tb[a:b]
        
        idx_b += [*range(a,b)]

        dists += list(tt[np.abs(tt)<limit])
        # dists[i*2*n: (i+1)*2*n - 1] = tb[idx-n:idx+n] - t
        
    return dists, idx_b


def get_distances_A_to_B_noq(ta, tb, n, offset_rough, limit):
    # dists = np.zeros(len(ta)*2*n)
    dists = []
    
    
    for i in range(len(ta)):
        t = ta[i]
        idx = find_nearest_idx(tb, t + offset_rough)
        
        a = max(0,idx - n) 
        b = min(idx+n+1,len(tb))
        
        tt = tb[a:b] - t

        dists += list(tt[np.abs(tt)<limit])
        # dists[i*2*n: (i+1)*2*n - 1] = tb[idx-n:idx+n] - t
        
    return dists



def get_distances_B_to_A(ta, tb, qa, qb, n, offset_rough, limit = 100):
    # dists = np.zeros(len(ta)*2*n)
    dists = []
    
    # index of shorlisted times in A
    idx_a = []
    
    for i in range(len(tb)):
        t = tb[i]
        idx = find_nearest_idx(ta, t - offset_rough)
        
        a = max(0,idx - n) 
        b = min(idx+n+1,len(ta))
        
        tt = t - ta[a:b]
        
        idx_a += [*range(a,b)]

        dists += list(tt[np.abs(tt)<limit])
        # dists[i*2*n: (i+1)*2*n - 1] = tb[idx-n:idx+n] - t
        
    return dists, idx_a


# def get_distances_moving(ta, tb, qa, qb, n, offset_rough, velocity_rough):
#     # velocity_rough must be in units of m/s
    
#     # dists = np.zeros(len(ta)*2*n)
#     dists = np.zeros((len(tb), 2*n))
    
#     drift = (tb*velocity_rough/299792458).astype(int)
    
#     for i in range(len(tb)):
#         t = tb[i]
#         idx = find_nearest_idx(ta, t - offset_rough - drift[i])
        
#         a = max(0,idx - n) 
#         b = min(idx+n,len(ta))
        
#         tt = t - ta[a:b]
#         # print('len tt', len(tt))
#         # idx_a += [*range(a,b)]

#         dists[i] = np.concatenate((tt, [0]*(8 - len(tt))))
#         # dists[i*2*n: (i+1)*2*n - 1] = tb[idx-n:idx+n] - t
        
#     return dists



def get_distances_moving2(ta, tb, qa, qb, n, offset_rough, velocity_rough, limit):
    # velocity_rough must be in units of m/s
    
    # dists = np.zeros(len(ta)*2*n)
    dists = []
    tb_idx = []
    qas = []
    
    drift = (tb*velocity_rough/299792458).astype(int)
    
    for i in range(len(tb)):
        t = tb[i]
        idx = find_nearest_idx(ta, t - offset_rough - drift[i])
        
        a = max(0,idx - n) 
        b = min(idx+n,len(ta))
        
        tt = t - ta[a:b]
        qt = qa[a:b]
        # print('len tt', len(tt))
        # idx_a += [*range(a,b)]
        
        mask = np.abs(tt - offset_rough)<limit
        a = list(tt[mask])
        if a:
            # dists += [list(tt[np.abs(tt-offset_rough)<limit])]
            dists += [a]
            tb_idx += [i]
            qas += [qt[mask]]
        
    return dists, tb_idx, qas



def get_distances_B_to_A_noq(ta, tb, n, offset_rough, limit):
    # dists = np.zeros(len(ta)*2*n)
    dists = []
    
    # index of shortlisted times in A
    idx_a = []
    
    for i in range(len(tb)):
        t = tb[i]
        idx = find_nearest_idx(ta, t - offset_rough)
        a = max(0,idx - n) 
        b = min(idx+n+1,len(ta))
        
        tt = t - ta[a:b]
        
        idx_a += [*range(a,b)]

        dists += list(tt[np.abs(tt)<limit])
        # dists[i*2*n: (i+1)*2*n - 1] = tb[idx-n:idx+n] - t
        
    return dists, idx_a


def get_distances_B_to_A_simple(ta, tb, window):
    # dists = np.zeros(len(ta)*2*n)
    dists = np.array([])
    
    for i in range(len(tb)):
        t = tb[i]
        
        taa = ta[np.abs(t - ta)<window]
        
        dists = np.append(dists, t - taa)
        
        
    return dists
        


def get_coincidences_A_to_B_noq(ta, tb, offset_rough, window):
    # dists = np.zeros(len(ta)*2*n)
    coinc = 0
    
    
    for i in range(len(ta)):
        t = ta[i]
        
        if (np.abs(t - tb) < window).any():
            coinc += 1
        
    return coinc       



def get_coincidences_A_to_B_wtags(ta, tb, tags, offset_rough, window):
    # dists = np.zeros(len(ta)*2*n)
    coinc = 0
    
    tagged = 0
    for i in range(len(ta)):
        t = ta[i]
        
        in_window = np.abs(t - tb) < window
        if in_window.any():
            coinc += 1
            
            inds = np.where(in_window)[0]
            
            tagged += tags[np.random.choice(inds)]
        
    return coinc, tagged