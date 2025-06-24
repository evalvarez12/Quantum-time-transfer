# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


nm10 = [78, 95, 99, 100, 100]
nm5 = [66, 93, 98, 99, 100]

rate = [100, 300, 500, 700, 900]

plt.figure()

plt.plot(rate, nm5, 'o', linewidth=3, label='5 nm filter')
plt.plot(rate, nm10, 'o', linewidth=3, label='10 nm filter')
plt.xlabel('Source rate (kHz)',fontsize=13)
plt.ylabel('Probability of success', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.legend()
plt.show()