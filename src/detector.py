#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 10:27:14 2023

@author: eduardo
"""

import numpy as np
import detector

def remove_dead_times_w_quantum(times, quantum, dead_time):
    if dead_time == 0:
        return times, quantum
    
    cleaned_times, cleaned_quantum = dead_time_erase_w_quantum(times, quantum, dead_time)
    if len(times) != len(cleaned_times) :
        cleaned_times, cleaned_quantum = remove_dead_times_w_quantum(cleaned_times, cleaned_quantum, dead_time)
    
    return cleaned_times, cleaned_quantum

def remove_dead_times_w_quantum_w_tags(times, quantum, tags, dead_time):
    if dead_time == 0:
        return times, quantum, tags
    
    cleaned_times, cleaned_quantum, cleaned_tags = dead_time_erase_w_quantum_w_tags(times, quantum, tags, dead_time)
    if len(times) != len(cleaned_times) :
        cleaned_times, cleaned_quantum, cleaned_tags = remove_dead_times_w_quantum_w_tags(cleaned_times, cleaned_quantum, cleaned_tags, dead_time)
    
    return cleaned_times, cleaned_quantum, cleaned_tags


def dead_time_erase_w_quantum(times, quantum, dead_time):
    mask = np.zeros(len(times)-1)
    
    # Mark the times that are within dead time of the count before
    mask[ (times[1:] - times[:-1]) <= dead_time] = 1

    # The first count is always valid
    mask = np.insert(mask, 0, 0)    
    # Mark only the  marked times that precede a valid count
    mask2 = np.zeros(len(times)-1)
    mask2[ (mask[:-1] - mask[1:]) == -1 ] = 1
    
    # The first count is always valid
    mask2 = np.insert(mask2, 0, 0)    
    return times[mask2 == 0], quantum[mask2 == 0]

def dead_time_erase_w_quantum_w_tags(times, quantum, tags, dead_time):
    mask = np.zeros(len(times)-1)
    
    # Mark the times that are within dead time of the count before
    mask[ (times[1:] - times[:-1]) <= dead_time] = 1

    # The first count is always valid
    mask = np.insert(mask, 0, 0)    
    # Mark only the  marked times that precede a valid count
    mask2 = np.zeros(len(times)-1)
    mask2[ (mask[:-1] - mask[1:]) == -1 ] = 1
    
    # The first count is always valid
    mask2 = np.insert(mask2, 0, 0)    
    return times[mask2 == 0], quantum[mask2 == 0], tags[mask2 == 0]


def remove_dead_times(times, dead_time):
    cleaned_times = dead_time_erase(times, dead_time)
    if len(times) != len(cleaned_times) :
        cleaned_times = remove_dead_times(cleaned_times, dead_time)
    
    return cleaned_times


def dead_time_erase(times, dead_time):
    mask = np.zeros(len(times)-1)
    
    # Mark the times that are within dead time of the count before
    mask[ (times[1:] - times[:-1]) <= dead_time] = 1

    # The first count is always valid
    mask = np.insert(mask, 0, 0)    
    # Mark only the  marked times that precede a valid count
    mask2 = np.zeros(len(times)-1)
    mask2[ (mask[:-1] - mask[1:]) == -1 ] = 1
    
    # The first count is always valid
    mask2 = np.insert(mask2, 0, 0)    
    return times[mask2 == 0]



def time_diff_reverse(a_times, b_times, n, n0):
    # Calculate the time differences based on alice times
    # for alice time finde the nearest bob time and go from there
    time_diffs = np.array([])
    
    for i in a_times:
        idx = np.searchsorted(b_times, i, 'left')
        try:
            ts = b_times[idx + n0 - (n-1) : idx + n0 + n]
        except IndexError:
            print('BREAKING', i)
            break
        time_diffs = np.concatenate((time_diffs, ts-i))
    return time_diffs
        
        
def time_diff_direct(alice_times, bob_times, n, n0):
    return
    # Calculate the time differences based on bob times
    # concatenate all the times together and identify the bob times
    # for each bob time finde the nearest bob time and go from there



# Debbuging code 

# t = np.array([1, 2, 3, 4, 5])

# print(dead_time_erase(t, 2))

# print(remove_dead_times(t, 2))

# print('-------------------------------')

# t = np.array([0, 10, 20, 35, 40, 50])

# print(dead_time_erase(t, 6))

# print(remove_dead_times(t, 6))

# print('-------------------------------')

# t = np.array([0, 10, 15, 21, 45, 50])

# print(dead_time_erase(t, 10))

# print(remove_dead_times(t, 10))
