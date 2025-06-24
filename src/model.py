# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:15:38 2023

@author: vil034
"""

import numpy as np
import detector

class Model:
    def __init__(self, source='entangled', tags=False):
        # Defaul model parameters corresponsing to a realistic setup
        self.distance = 100000 # 10 km separation Alice and Bob
        
        # Source rate
        self.entangled_pair_rate = 3
        
        # Loss
        self.detector_loss = 4
        self.channel_loss = 25
        
        # Jitter
        self.detector_jitter = 350e-6
        self.atmosphere_turbulence_jitter = 20e-6
        self.time_tag_jitter = 42e-6
        
        self.detector_dead_time = 30e-3
        
        self.background_counts_radiation_rate_sec =  5e6 # In seconds
        self.dark_count_rate_sec = 500 # In seconds
        
        self.qber = .01 
        
        self.source = source
        
        self.tags = tags
        self.tags_data = []
        
        
    def set_loss(self, loss, detector_loss = 0):
        self.detector_loss = detector_loss
        self.channel_loss = loss
        
        
    def set_background_rate(self, background_counts_rate, dark_count_rate_sec=0):
        self.background_counts_radiation_rate_sec =  background_counts_rate * 1e6 # In seconds
        self.dark_count_rate_sec = 0
        
        
    def set_qber(self, qber):
        self.qber = qber
        
    def set_detector_dead_time(self, dead_time):
        self.detector_dead_time = dead_time
        
    def set_jitter(self, jitter):
        # Jitter
        self.detector_jitter = 0
        self.atmosphere_turbulence_jitter = 0
        self.time_tag_jitter = jitter
        

    def set_signal_rate(self, rate):
        self.entangled_pair_rate = rate

    def generate_time_data(self, time_shift, start_time, total_protocol_time_sec, verbose=False):
        ############ System parameters
        # All time is set in unit of microsecond 1e-6
        
        # detector_time_res = 50e-6  # time resolution of the detector - 50ps
        
        # total_protocol_time_sec = .001
        total_protocol_time = total_protocol_time_sec * 1e6 # protocol duration in microsecond
     
        N_created_pairs = int(self.entangled_pair_rate * total_protocol_time)
        
        # pair_generation_time = 1/entangled_pair_rate
        
        # c = 299792458 # speed light m/s

        # light_travel_time = self.distance/c # light travel time in seconds
        # light_travel_time = light_travel_time*1e6 # light travel time in microsecond
        
        light_travel_time = time_shift
        
         
        N_background_counts = (self.background_counts_radiation_rate_sec + self.dark_count_rate_sec)*total_protocol_time_sec
        N_background_counts = int(N_background_counts)
        
        
  
        background_count_times = np.random.rand(N_background_counts)*total_protocol_time + start_time
        
        
        ######################## Creating the time samples
        
        
        # Alice and Bob parameters
        alice_loss_dB = self.detector_loss # Loss on Alice detection devices - in dB
        alice_transmissivity = 10**(-alice_loss_dB/10)
        
        bob_loss_dB = self.detector_loss + self.channel_loss
        bob_transmissivity = 10**(-bob_loss_dB/10)
        
        alice_jitter = np.sqrt(self.detector_jitter**2 + self.time_tag_jitter**2)
        bob_jitter = np.sqrt(self.detector_jitter**2 + self.atmosphere_turbulence_jitter**2 + self.time_tag_jitter**2)
        
        # Entangled pairs are generated and on photon is measured by Alices keeping a timestamp
        N_created_pairs = int(self.entangled_pair_rate * total_protocol_time)
            
        entangled_pairs_time = np.random.rand(N_created_pairs)*total_protocol_time + start_time    
        entangled_pairs_time = np.sort(entangled_pairs_time)
        
        # Create quantum information of the entangled pairs
        entangled_pairs_quantum = np.random.choice([1, -1], len(entangled_pairs_time))
        
        
        # Alice times 
        photons_loss_roll_alice = np.random.rand(N_created_pairs)
        mask = photons_loss_roll_alice < alice_transmissivity
        alice_times = entangled_pairs_time[mask]
        
        # Alice quantum information
        alice_quantum = entangled_pairs_quantum[mask]
        
       
                                         
        # Bob entangled pair times & quantum info
        photons_loss_roll_bob = np.random.rand(N_created_pairs)
        mask = photons_loss_roll_bob < bob_transmissivity
        bob_pair_times = entangled_pairs_time[mask]

        
        bob_quantum = entangled_pairs_quantum[mask]
        
        # Add qber error to Bobs quantum info
        errors = np.ones(len(bob_quantum))
        error_roll = np.random.rand(len(bob_quantum))
        errors[error_roll < self.qber] *= -1
        bob_quantum *= errors.astype(int)
    
        # Count the coincidences
        # coincidences = 0
        # for ti in bob_pair_times:
        #     if ti in alice_times:
        #         coincidences += 1
        
        # Add time jitter to Bob source times
        jitter_roll_alice = np.random.normal(0, alice_jitter, len(alice_times))
        alice_times = alice_times + jitter_roll_alice
        
        jitter_roll_bob = np.random.normal(0, bob_jitter, len(bob_pair_times))
        bob_pair_times = bob_pair_times + jitter_roll_bob + light_travel_time
        
    
        
        bob_times = np.concatenate((bob_pair_times, background_count_times))
        sort_ind = np.argsort(bob_times)
        bob_times = bob_times[sort_ind]
        
        # bob_quantum_thermal = np.zeros(N_background_counts)
        bob_quantum_thermal = np.random.choice([1, -1], N_background_counts) 
        bob_quantum = np.concatenate((bob_quantum, bob_quantum_thermal))
        bob_quantum = bob_quantum[sort_ind]
        
        if self.tags:
            bob_pair_tags = np.zeros_like(bob_pair_times)
            bob_background_tags = np.ones_like(background_count_times)
            bob_tags = np.concatenate((bob_pair_tags, bob_background_tags))
            bob_tags = bob_tags[sort_ind]
            
            bob_times, bob_quantum, bob_tags = detector.remove_dead_times_w_quantum_w_tags(bob_times, bob_quantum, bob_tags, self.detector_dead_time)
            
            self.tags_data = bob_tags
        else:
            bob_times, bob_quantum = detector.remove_dead_times_w_quantum(bob_times, bob_quantum, self.detector_dead_time)
            
            
       
        # Remove counts that fall within the detector dead time - see detector.py
        alice_times, alice_quantum = detector.remove_dead_times_w_quantum(alice_times, alice_quantum, self.detector_dead_time)
        
        
        
        if verbose:
            print('Protocol time: ', total_protocol_time)
            print('Entangled pairs generated: ', N_created_pairs)
            
            print('Alice jitter: ', alice_jitter)
            print('Bob jitter: ', bob_jitter)
            
            total_jitter = np.sqrt(alice_jitter**2 +bob_jitter**2)
            print('Total jitter: ', total_jitter)
            
            print('Alice pairs detected: ', len(alice_times))
            print('Bob pairs detected: ', len(bob_pair_times))
            print('Bob background counts: ', N_background_counts)
            print('Bob total counts: ', len(bob_times))
            
            # Approximated - real number is less because of detector dead time
            # print('Num signal photons: ', N_created_pairs*bob_transmissivity* alice_transmissivity)
            # total_background_rate = (self.background_counts_radiation_rate_sec + self.dark_count_rate_sec)*1e-6
            # print('Num background noise: ', len(alice_times)*total_background_rate*2*total_jitter)
            # print('S: ', N_created_pairs*bob_transmissivity* alice_transmissivity/(total_background_rate*4*total_jitter))             
            # print('Total coincidences: ', coincidences) 
             
            
            
            
            print('--------------------------------')
    
            
            
        return alice_times, bob_times, alice_quantum, bob_quantum
    
    
    