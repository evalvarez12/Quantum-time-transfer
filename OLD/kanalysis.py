#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:23:09 2023

@author: eduardo
"""

import numpy as np
from bisect import bisect

def find_peak(alice_times, bob_times, n):
    taus = []
    for ai in alice_times:
        index = bisect(bob_times, ai)
        bi = bob_times[index-n:index+n]
        taus = taus + list(ai-bi)
        
        
    return np.array(taus)

def find_peak_reverse(alice_times, bob_times, n):
    taus = []
    for bi in bob_times:
        index = bisect(alice_times, bi)
        ai = alice_times[index-n:index+n]
        taus = taus + list(bi-ai)
        
        
    return np.array(taus)

