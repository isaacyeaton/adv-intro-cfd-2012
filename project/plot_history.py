"""Tool to plot the history of a run using the history.dat file.
"""

import numpy as np
import matplotlib.pyplot as plt

f = 'MMS_PJ_129x129.dat'
d = np.genfromtxt(f, skip_header=2)
resp, resu, resv = d[:, 2], d[:, 3], d[:, 4]
n, t = d[:, 0], d[:, 1]

fig, ax = plt.subplots()
ax.semilogy(n, resp, label='pressure')
ax.semilogy(n, resu, label='u')
ax.semilogy(n, resv, label='v')

ax.legend(loc='upper right', frameon=False)
ax.set_xlabel('iteration')
ax.set_ylabel('relative iterative residual')


plt.show()
