#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:20:51 2023

@author: eduardo
"""

import numpy as np
import matplotlib.pyplot as plt
import gauss_fit as gf

mean = 0.5
std = 0.03
sample = 5000

bins = 1000

integration_time = 3.5

background_counts = 500
time_data = np.random.normal(mean, std, sample)

background_data = np.random.rand(background_counts) * integration_time

time_data = np.concatenate((time_data, background_data))

hist_data=np.histogram(time_data, bins)


x_data = hist_data[1][1:]
y_data = hist_data[0]


params = [0.5, .1]
popt = gf.fit(x_data, y_data, params, plot=False)

print("Before tampering")
# print('fitted_amplitude= {}counts'.format(popt[0]))
print('time difference = {}s'.format(popt[1]))
print('time resolution = {}s'.format(popt[2]))

og_offset = popt[1]

# Tapering percentage
p = 0.05

sample_eve = int(sample*p)

time_data_eve = np.random.normal(mean + 3*std, std, sample_eve)

time_data = np.concatenate((time_data, time_data_eve))

hist_data=np.histogram(time_data, bins)


x_data = hist_data[1][:-1]
y_data = hist_data[0]


params = [0.5, .1]
popt = gf.fit(x_data, y_data, params, plot=True)

print("Tapered")
# print('fitted_amplitude= {}counts'.format(popt[0]))
print('time difference = {}s'.format(popt[1]))
print('time resolution = {}s'.format(popt[2]))

final_offset = popt[1]

print('Tamper effect = {}s'.format(np.abs(og_offset - final_offset)))
