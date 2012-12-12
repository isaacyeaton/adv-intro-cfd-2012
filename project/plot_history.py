"""Tool to plot the history of a run using history.tec files.
   This file is used for make some plots in the report of iterative
   residuals.
"""

import numpy as np
import matplotlib.pyplot as plt


## Where to save things
save_dir = '/home/isaac/Desktop/adv-intro-cfd-2012/project/report/figs/'
save_figures = True
gr = 1.61803398875  # golden ratio


def save_fig(fig, name):
    """Save a figure."""
    
    if save_figures:
        fig.savefig(save_dir + name + '.eps', bbox_inches='tight', pad_inches=0)
    else:
        print('No saving figure {0}'.format(name))
        
    
def load_data(f):
    """Load in the data file.
    """
    
    d = np.genfromtxt(f, skip_header=2)
    resp, resu, resv = d[:, 2], d[:, 3], d[:, 4]
    n, t = d[:, 0], d[:, 1]
    
    return np.column_stack([n, resp, resu, resv])


##### Iterative residual for MMS

base = './out_MMS/'
fbase_sg = base + 'history_SGS_{0}x{0}_Re=10_cfl=0.9.tec'
fbase_pj = base + 'history_PJ_{0}x{0}_Re=10_cfl=0.9.tec'

fig, ax = plt.subplots()
for node, color in zip([17, 33, 65, 129, 257], ['b', 'g', 'r', 'k', 'm']):
    
    d_sg = load_data(fbase_sg.format(node))
    d_pj = load_data(fbase_pj.format(node))
    
    ax.semilogy(d_sg[:, 0], d_sg[:, 2], color, lw=2, label='{0}x{0}: SGS'.format(node))
    ax.semilogy(d_pj[:, 0], d_pj[:, 2], color+'--', lw=2, label='{0}x{0}: PJ'.format(node))  

ax.set_xlim(0, 1.5e5)
ax.set_ylim(1e-12, 1)
ax.legend(loc='upper right', frameon=False, fancybox=True, ncol=2,
          title='Driven cavity MMS: Re=10, CFL=0.9', prop={'size':10})
ax.set_xlabel('iteration')
ax.set_ylabel('relative iterative residual')
#fig.set_size_inches(5*gr, 5)
save_fig(fig, 'MMS_history')


###### Time deriv pre-conditioning kappa

kappa_base = './out_MMS/history_rkappa={0}_SGS_65x65_Re=10_cfl=0.5.tec'
kappa_vals = ['.01', '.1', '.5', '.9']
fig, ax = plt.subplots()
for kappa, sty in zip(kappa_vals, ['b', 'g', 'r--', 'k']):
    d = load_data(kappa_base.format(kappa))
    ax.semilogy(d[:, 0], d[:, 1], sty, lw=2, label=r'$\kappa$ = {0}'.format(kappa))
    #ax.semilogy(d[-1, 0], d[-1:, 1], sty, ms=8, lw=2, label=r'$\kappa$ = {0}'.format(kappa))
ax.legend(loc='upper right', frameon=False, fancybox=True, ncol=2,
          title='Driven cavity MMS: Re=10, \nCFL=0.5, 65x65 nodes', prop={'size':12})
ax.set_xlabel('iteration')
ax.set_ylabel('relative iterative residual')
ax.set_ylim(1e-11, 1e1)
#save_fig(fig, 'MMS_kappa')

#kappa_base = './out_MMS/history_rkappa={0}_SGS_65x65_Re=10_cfl=0.5.tec'
#kappa_vals = ['.01', '.1', '.5', '.9']
#fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=True)
#for kappa, sty in zip(kappa_vals, ['bo', 'gs', 'r*', 'm^']):
#    d = load_data(kappa_base.format(kappa))
#    
#    # p
#    ax1.semilogy(d[:, 0], d[:, 1], sty[0], lw=2, label=r'$\kappa$ = {0}'.format(kappa))
#    ax1.semilogy(d[-1, 0], d[-1:, 1], sty, ms=8, lw=2, label=r'$\kappa$ = {0}'.format(kappa))
#    
#    # u
#    ax2.semilogy(d[:, 0], d[:, 2], sty[0], lw=2, label=r'$\kappa$ = {0}'.format(kappa))
#    ax2.semilogy(d[-1, 0], d[-1:, 2], sty, ms=8, lw=2, label=r'$\kappa$ = {0}'.format(kappa))
#    
#    # v
#    ax3.semilogy(d[:, 0], d[:, 3], sty[0], lw=2, label=r'$\kappa$ = {0}'.format(kappa))
#    ax3.semilogy(d[-1, 0], d[-1:, 3], sty, ms=8, lw=2, label=r'$\kappa$ = {0}'.format(kappa))
#
#ax1.legend(loc='upper right', frameon=False, fancybox=True, ncol=2,
#          title='Driven cavity MMS: Re=10, pressure,\nCFL=0.5, 65x65 nodes', prop={'size':11})
#ax2.legend(loc='upper right', frameon=False, fancybox=True, ncol=2,
#          title='Driven cavity MMS: Re=10, u-velocity,\nCFL=0.5, 65x65 nodes', prop={'size':11})
#ax3.legend(loc='upper right', frameon=False, fancybox=True, ncol=2,
#          title='Driven cavity MMS: Re=10, v-velocity,\nCFL=0.5, 65x65 nodes', prop={'size':11})
#ax1.set_xlabel('iteration')
#ax1.set_ylabel('relative iterative residual')
#ax1.set_ylim(1e-12, 1e-1)

plt.show()







