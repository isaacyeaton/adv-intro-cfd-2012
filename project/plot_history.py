"""Tool to plot the history of a run using the history.dat file.
"""

import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt('history.dat', skip_header=2)
resp, resu, resv = d[:, 2], d[:, 3], d[:, 4]
n, t = d[:, 0], d[:, 1]

fig, ax = plt.subplots()
ax.semilogy(n, resp, label='pressure')
ax.semilogy(n, resu, label='u')
ax.semilogy(n, resv, label='v')

ax.set_xlabel('iteration')
ax.set_ylabel('relative iterative residual')


plt.show()
