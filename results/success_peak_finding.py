# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 09:57:00 2023

@author: vil034
"""

import numpy as np
import math
from scipy import special


def p_peak(val, p, N, dt):
    c = math.comb(N, val)
    
    return c*(p**val)*((1-p)**(N-val))


rate_a = 250e6
rate_b = 5e5

N = 2**25
total_time = .05
dt = total_time/N

loss = 30
t = 10**(-loss/10)

signal = rate_a*total_time*t

p = rate_a*rate_b*(dt**2)

noise_avg = p*N
noise_var = N*p*(1-p)


prob = .5*special.erfc((signal - noise_avg)/np.sqrt(2*noise_var))
print(prob*N)



# print(p_peak(signal, p, N, dt))



