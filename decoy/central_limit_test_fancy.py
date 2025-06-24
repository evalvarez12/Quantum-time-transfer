# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


plt.close('all')

Nl = int(5e6)

##### Lognormal dist
# m = -3
# s = .7

# mean = np.exp(m + s**2/2)
# var = (np.exp(s**2) - 1) * np.exp(2*m + s**2)
# print('mean: ', mean)
# print('var: ', var)
# print('-----------------------')

# loss = np.random.lognormal(m, s, Nl)
# loss = loss[loss < 1]

# lossdb = -10*np.log10(loss)
# plt.hist(lossdb, 100)
# plt.show()


# f = lambda n : np.random.lognormal(m, s, n)

##### Normal dist
mean = 0.5
var = 0.3
print('mean: ', mean)
print('var: ', var)
print('-----------------------')

loss = np.random.normal(mean, var, Nl)
loss = loss[loss < 1]
loss = loss[loss > 0]
lossdb = -10*np.log10(loss)
plt.hist(lossdb, 100)
plt.show()

f = lambda n : np.random.normal(mean, var, n)


###############
n_it = 1
n = 1000000
ps = []
for i in range(n_it):

    # loss_sample = np.random.choice(loss, n)
    
    # loss = np.random.lognormal(m, s, n)
    loss = f(n)
    loss = loss[loss < 1]
    loss = loss[loss > 0]
    
    n = len(loss)
    
    print('Mean loss: ', np.mean(loss))
    print('Var loss: ', np.std(loss)**2)
    
    clicks = 0
    
    roll = np.random.rand(n)
    clicks = roll < loss
    
    p = sum(clicks)/n
    # print('Total clicks: ', sum(clicks))
    # print('Probability: ', p)
    
    # print('Distance from pop mean: ', p - mean)
    # print('Distance variance expected: ', np.sqrt(p*(1-p)/N))
    
    
    ps += [p]

      
ps = np.array(ps)  

plt.figure()
plt.hist(ps, 100)

print('Average: ', np.mean(ps))
print('Var: ', np.std(clicks)**2)
print(np.mean(ps)*(1-np.mean(ps)))
print(np.mean(ps)*(1-np.mean(ps)) + var)

# print(np.sqrt(np.mean(ps)*(1-np.mean(ps))/(n_it)))

######### Fancy analysis M. Curty et al.
def g(x,y):
    return np.sqrt(2*x*np.log(1/y))

e = 10e-4/3
e1 = 10e-4/3
e2 = 10e-4/3

x = sum(clicks)

ml = x - np.log(n*np.log(1/e)/2)

print('Condition 1: ', (2/e1)**(1/ml), np.exp(3/(4*np.sqrt(2))**2), (2/e1)**(1/ml) < np.exp(3/(4*np.sqrt(2))**2))
print('Condition 2: ', (1/e2)**(1/ml), np.exp(1/3), (1/e2)**(1/ml) < np.exp(1/3))

D1 = g(x, e1**4/16)
D2 = g(x, e2**(3/2))

print('x:', x)

print('D1: ', D1)
print('D2: ', D2)
print('Total err prob: ', e + e1 + e2)


print(sum(loss))

print(x/n, D1/n, D2/n)