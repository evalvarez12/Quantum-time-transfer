# -*- coding: utf-8 -*-


import numpy as np

def beam_waist(w0, lamb, L):
    zr = np.pi*w0**2/lamb
    return w0 * np.sqrt(1 + (L/zr))


def link_length(h0, H, z):
    R = 6371e3;
    l = np.sqrt(((R+h0)*np.cos(z))**2 + (H+R)**2 - (R+h0)**2) - (R+h0)*np.cos(z);
    return l

# Downlink
w0 = 0.3
det = 1.2
H = 500e3
lamb = 810e-9
z = np.deg2rad(60)

L = link_length(0, H, z)
w = beam_waist(w0, lamb, L)

w2 = w0 + 2*L*np.tan(15e-6/2)

print(w,w2)
print((det/w)**2)
print((det/w2)**2)
