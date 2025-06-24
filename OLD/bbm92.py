# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 16:01:22 2023

@author: vil034
"""

import numpy as np

def chooseBasis(numBits):
    """Return Alice and Bob's randomly chosen mstment axes for the specified
       number of qubits in the BBM92 protocol:
           A chooses from (0, pi/4, pi/2) with equal probability,
           B chooses from (pi/4, pi/2, 3pi/4) with equal probability.
    """
    choicesA = [0, 1]
    choicesB = [0, 1]

    basesA = np.random.choice(choicesA, numBits)
    basesB = np.random.choice(choicesB, numBits)
    
    return basesA, basesB

def measureEntangledState(basisA, basisB, qber=0.0):
    """Return Alice and Bob's measurement results on a pair of maximally
    entangled qubits. basis[A,B] contain Alice and Bob's axes of mstment.
    """
    
    l = len(basisA)
    # Alice measures either basis state with equal probability
    # -1 or 1
    resultsA = np.random.choice([-1, 1], l)

    # If Alice and Bob chose the same basis, Bob's result is
    # perfectly anti-correlated with Alice's. Otherwise Bob's
    # results are random.
    
    r = np.random.choice([0, 1], l)
    r = (-1) ** (np.abs(basisA - basisB)*r)

    resultsB = -resultsA * r

    if qber:
        sample_err = np.random.rand(2, l)
        resultsA[sample_err[0] < qber ] *= -1 
        resultsB[sample_err[1] < qber ] *= -1 


    return resultsA, resultsB



def matchData(dataA, dataB, basesA, basesB):
    keyA = dataA[basesA == basesB]
    keyB = dataB[basesA == basesB]
    
    return keyA, keyB


def runProtocol(numSignals, errorRate=0.0, verbose=True):
    """Simulation of BBM92 entanglement-based protocol for quantum key distribution."""

    if verbose:
        print("\n=====BBM92 protocol=====\n%d initial signals" % (numSignals))
        if errorRate: print("channel noise: %d" % (errorRate) )

    # An entangled state is generated in the singlet state
    #     +0.7071 |0> -0.7071 |1>,
    #  one particleis sent to Alice and the other to Bob.
    # Alice and Bob each choose randomly their basis for each measurement from [Z, X]

    basesA, basesB = chooseBasis(numSignals)
    
    if verbose:
        print('Alice basis: ', str(basesA))
        print('Bob basis:   ', str(basesB))
    

    dataA, dataB = measureEntangledState(basesA, basesB, errorRate)

    if verbose:
        print('Measurements are made')
        print('Alice data: ', str(dataA))
        print('Bob data:   ', str(dataB))
        
    
    keyA, keyB = matchData(dataA, dataB, basesA, basesB)
    
    if verbose:
        print('Alice and Bob anounce their bases')
        print('Alice raw key: ', str(keyA))
        print('Bob raw key:   ', str(keyB))
        
    print('Estimation of the QBER still needs to be done')
    
    return keyA