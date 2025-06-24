
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 12:05:42 2023

@author: vil034
"""

import numpy as np
import detector




def generate_time_data(time_shift, total_protocol_time_sec, loss=True, background_arg=0, jitter=True, verbose=False):
    ############ System parameters
    # All time is set in unit of microsecond 1e-6
    
    # detector_time_res = 50e-6  # time resolution of the detector - 50ps
    
    start_time = 0 
    # total_protocol_time_sec = .001
    total_protocol_time = total_protocol_time_sec * 1e6 # protocol duration in microsecond
    entangled_pair_rate = 3 # 3 MHz source rate
    N_created_pairs = int(entangled_pair_rate * total_protocol_time)
    
    # pair_generation_time = 1/entangled_pair_rate
    
    c = 299792458 # speed light m/s
    distance_alice_bob = 100000 # 10 km separation Alice and Bob
    light_travel_time = distance_alice_bob/c # light travel time in seconds
    light_travel_time = light_travel_time*1e6 # light travel time in microsecond
    
    light_travel_time = time_shift
    
    
    # Loss
    detector_loss = 4
    channel_loss = 25
    
    # Jitter
    detector_jitter = 350e-6
    atmosphere_turbulence_jitter = 20e-6
    time_tag_jitter = 42e-6 # Could be 8e-6
    
    
    detector_dead_time = 30e-3
    
    
    background_counts_radiation_rate_sec =  500e3 # In seconds
    dark_count_rate_sec = 500 # In seconds
     
    N_background_counts = (background_counts_radiation_rate_sec + dark_count_rate_sec)*total_protocol_time_sec
    N_background_counts = int(N_background_counts)
    
    
    if background_arg:
        N_background_counts = background_arg

    
    background_count_times = np.random.rand(N_background_counts)*total_protocol_time + start_time
    
    
    ######################## Creating the time samples
    
    
    # Alice and Bob parameters
    alice_loss_dB = detector_loss # Loss on Alice detection devices - in dB
    alice_transmissivity = 10**(-alice_loss_dB/10)
    
    bob_loss_dB = detector_loss + channel_loss
    bob_transmissivity = 10**(-bob_loss_dB/10)
    
    alice_jitter = detector_jitter + time_tag_jitter
    bob_jitter = detector_jitter + atmosphere_turbulence_jitter + time_tag_jitter
    
    if not jitter:    
        alice_jitter = 0
        bob_jitter = 0
    
    if not loss:
        alice_transmissivity = 1
        bob_transmissivity = 1
    
    entangled_pairs_time = np.random.rand(N_created_pairs)*total_protocol_time + start_time    
    entangled_pairs_time = np.sort(entangled_pairs_time)
    
    # Alice times 
    photons_loss_roll_alice = np.random.rand(N_created_pairs)
    alice_times = entangled_pairs_time[photons_loss_roll_alice < alice_transmissivity] 
    

    # Bob entangled pair times
    photons_loss_roll_bob = np.random.rand(N_created_pairs)
    bob_pair_times = entangled_pairs_time[photons_loss_roll_bob < bob_transmissivity]

    # Count the coincidences
    # coincidences = 0
    # for ti in bob_pair_times:
    #     if ti in alice_times:
    #         coincidences += 1
    
    # Add time jitter to Bob source times
    jitter_roll_alice = np.random.normal(0, alice_jitter, len(alice_times))
    alice_times = alice_times + jitter_roll_alice
    
    jitter_roll_bob = np.random.normal(light_travel_time, bob_jitter, len(bob_pair_times))
    bob_pair_times = bob_pair_times + jitter_roll_bob
    
    
    bob_times = np.concatenate((bob_pair_times, background_count_times))
    bob_times = np.sort(bob_times)
    
    # Remove counts that fall within the detector dead time - see detector.py
    bob_times = detector.remove_dead_times(bob_times, detector_dead_time)
    alice_times = detector.remove_dead_times(alice_times, detector_dead_time)
    
    # Approximated - real number is less because of detector dead time
    N_coincidences = int(np.floor(N_created_pairs*alice_transmissivity*bob_transmissivity)) 
    
    if verbose:
        print('Protocol time: ', total_protocol_time)
        print('Entangled pairs generated: ', N_created_pairs)
        
        print('Alice jitter: ', alice_jitter)
        print('Bob jitter: ', bob_jitter)
        print('Total jitter: ', alice_jitter + bob_jitter)
        
        print('Alice pairs detected: ', len(alice_times))
        print('Bob pairs detected: ', len(bob_pair_times))
        print('Bob background counts: ', N_background_counts)
        print('Bob total counts: ', len(bob_times))
        # print('Total coincidences: ', coincidences) 
        # print('Total coincidences predicted: ', N_coincidences) 

        
        
    return alice_times, bob_times, N_coincidences, N_background_counts

