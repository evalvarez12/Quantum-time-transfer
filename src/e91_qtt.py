# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:10:00 2023

@author: vil034
"""

import qutip as qt
import numpy as np

def angle_to_qubit(theta):
    return np.cos(theta)*qt.basis(2, 0) + np.sin(theta)*qt.basis(2,1)

# In phyisical space - angles of qubit space / 2
anglesA = [0, np.pi/4, 3*np.pi/8]
anglesB = [np.pi/4, 3*np.pi/8, 5*np.pi/8]

statesA = [[angle_to_qubit(i), angle_to_qubit(i+np.pi/2)] for i in anglesA]
statesB = [[angle_to_qubit(i), angle_to_qubit(i+np.pi/2)] for i in anglesB]