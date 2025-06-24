# -*- coding: utf-8 -*-


import numpy as np


x = 0.006
n = 1e6

x = int(x*n)

e = 0.003
e2 = 0.003


ml = x - np.sqrt(n*np.log(1/e)/2)

if (1/e2)**(1/ml) > np.exp(1/3):
    print('Condition not satisfied upper bound : ', x)
    # return 0

D = np.sqrt(2*x*np.log(1/(e2**(3/2))))

xu = x + D
print(xu)