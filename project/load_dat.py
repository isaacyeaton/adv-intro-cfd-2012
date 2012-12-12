
import numpy as np
import matplotlib.pyplot as plt
from streamplot import streamplot

def run_num(string):
    """Run number from index."""

    n = np.int(string.split('=')[-1].split('"')[0])
    return n


def breakup_file(s):
    """Break-up the file for where different iterations occur.
    
    Note: needs a little more work to be finished up.'
    """
    
    dstart, iterations = [], []
    for cnt, line in enumerate(s):
        
        if line.startswith('DATAPACKING'):
            dstart.append(cnt+1)
    
    return dstart
    

def load_tec(fname):
    """Load in a tecplot file.
    """
    
    def _reshape(dset):
        """Reshape the data set."""
        return dset.reshape(node, node)
    
    if not fname.endswith('.tec'):
        fname += '.tec'
    
    s = open(fname, 'r').read().split('\n')
    nlines = len(s)
    node = np.int32(s[3].split(' ')[1])  # 'I= 33 J= 33'
    dstart = nlines - node**2 - 1
    ss = s[dstart:-1]  # last line might have some junk/CR
    
    nrows, ncols = len(ss), len(ss[0].split(' '))
    d = np.zeros((nrows, ncols))
    
    for i, line in enumerate(ss):
        d[i] = np.array(line.split(' ')).astype(np.float64)
    
    out = np.zeros((node, node, ncols))

    for i in range(ncols):
        out[:, :, i] = _reshape(d[:, i])
    
    return node, out


def save_fig(fig, name):
    """Save a figure."""
    
    if save_figures:
        fig.savefig(save_dir + name + '.eps', bbox_inches='tight', pad_inches=0)
    else:
        print('No saving figure {0}'.format(name))


def plot_mms(dd):
    """Make contour plots of the manufactured solution."""
    
    cmap = plt.cm.RdBu_r

    # U-velocity
    fig = plt.figure()
    cax = plt.contourf(dd[:, :, 0], dd[:, :, 1], dd[:, :, 3], np.r_[.4:1.05:.025], cmap=cmap)
    cbar = fig.colorbar(cax, orientation='vertical')
    cbar.set_label('u, m/s')
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    #if save_figures:
    #    plt.savefig(save_dir + 'MMS_u.eps', bbox_inches='tight', pad_inches=0)
    
    # P
    plt.figure()
    cax = plt.contourf(dd[:, :, 0], dd[:, :, 1], dd[:, :, 2], np.r_[.2:1.25:.025], cmap=cmap)
    cbar = plt.colorbar()
    cbar.set_label(r'p, N/m$^2$')
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    #if save_figures:
    #    plt.savefig(save_dir + 'MMS_p.eps', bbox_inches='tight', pad_inches=0)
    
    # V-velocity
    plt.figure()
    cax = plt.contourf(dd[:, :, 0], dd[:, :, 1], dd[:, :, 4], np.r_[.5:.73:.005], cmap=cmap)
    cbar = plt.colorbar()
    cbar.set_label('v, m/s')
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    #if save_figures:
    #    plt.savefig(save_dir + 'MMS_v.eps', bbox_inches='tight', pad_inches=0)


def ooa(de2, de1, r=2.):
    """Observed order of accuracy.

    Inputs:
        de2 (float) - DE error of coarse mesh
        de1 (float) - DE error of fine mesh
        r (float, default=2.) - grid refinement factor

    Outputs:
        phat (float) - observed order of accuracy
    """
    
    return np.log(de2 / de1) / np.log(r)


def find_phat(DE, r=2.):
    """Find the observed order of accuracy.
    
    Note: requires coarse to fine mesh
    
    Inputs:
        DE (array) - discretization error using L2 and Linf norms
        r (float, default=2.) - refinement factor

    Outputs:
        phat (array) - observed order of accuracy
    """

    nruns = DE.shape[0] - 1
    phat = np.zeros((nruns, DE.shape[1]))
    for i in range(nruns):
        # rows = L1, L2, and Linf norms
        phat[i] = ooa(DE[i, 0], DE[i+1, 0]), \
                  ooa(DE[i, 1], DE[i+1, 1]), \
                  ooa(DE[i, 2], DE[i+1, 2])

    return phat



save_dir = '/home/isaac/Desktop/adv-intro-cfd-2012/project/report/figs/'
save_figures = True
gr = 1.61803398875  # golden ratio

### Make contour plots of manufactured solution
n, dd = load_tec('./out_MMS/cavity_SGS_257x257_Re=10_cfl=0.9.tec')
plot_mms(dd)
###

files = '''
cavity_SGS_129x129_Re=10_cfl=0.9.tec
cavity_SGS_17x17_Re=10_cfl=0.9.tec
cavity_SGS_257x257_Re=10_cfl=0.9.tec
cavity_SGS_33x33_Re=10_cfl=0.9.tec
cavity_SGS_65x65_Re=10_cfl=0.9.tec'''
files = files.split('\n')

f = './out_cavity/cavity_SGS_33x33_Re=100_cfl=.05'
f = './out_cavity/cavity_SGS_129x129_Re=500_cfl=0.4'
f = './out_MMS/' + files[2]
n, dd = load_tec(f)


##### DE norms for MMS ##################################################

mms_temp = './out_MMS/cavity_SGS_{0}x{0}_Re=10_cfl=0.9.tec'
nodes = [5, 9, 17, 33, 65, 129, 257]
DEp = np.zeros((len(nodes), 3))  # L1, L2, Linf
DEu = np.zeros((len(nodes), 3))
DEv = np.zeros((len(nodes), 3))

L1 = lambda xx: np.sum(np.abs(xx.flatten()))/ xx.size
L2 = lambda xx: np.sqrt(np.sum(np.abs(xx.flatten()**2))/ xx.size)
Linf = lambda xx: np.linalg.norm(xx.flatten(), np.inf)

#L1 = lambda xx: np.linalg.norm(xx.flatten(), 1) / xx.size
#L2 = lambda xx: np.linalg.norm(xx.flatten(), 2) / xx.size
#Linf = lambda xx: np.linalg.norm(xx.flatten(), np.inf)

for i, node in enumerate(nodes):
    n, dd = load_tec(mms_temp.format(node))
    
    # normalized DE errors (dividing by exact value)
    dp = dd[:, :, 8] / dd[:, :, 5]
    du = dd[:, :, 9] / dd[:, :, 6]
    dv = dd[:, :, 10] / dd[:, :, 7]
    
    DEp[i] = L1(dp), L2(dp), Linf(dp)
    DEu[i] = L1(du), L2(du), Linf(du)
    DEv[i] = L1(dv), L2(dv), Linf(dv)

# grid refinement factor
hh = np.asarray(nodes, dtype=np.float64) - 1
h = hh[-1] / hh

fig, ax = plt.subplots()
ax.loglog(h, DEp[:, 0], 'bo-', lw=2, label='L1 p')
ax.loglog(h, DEp[:, 1], 'bs--', lw=2, label='L2 p')
ax.loglog(h, DEp[:, 2], 'b*:', lw=2, label='Linf p')
ax.loglog(h, DEu[:, 0], 'go-', lw=2, label='L1 u')
ax.loglog(h, DEu[:, 1], 'gs--', lw=2, label='L2 u')
ax.loglog(h, DEu[:, 2], 'g*:', lw=2, label='Linf u')
ax.loglog(h, DEv[:, 0], 'ro-', lw=2, label='L1 v')
ax.loglog(h, DEv[:, 1], 'rs--', lw=2, label='L2 v')
ax.loglog(h, DEv[:, 2], 'r*:', lw=2, label='Linf v')

ax.legend(loc='lower right', frameon=False, fancybox=True, ncol=3,
          title='Driven cavity MMS: Re=10, \nCFL=0.9, finest mesh=257x257',
          prop={'size':10})
ax.set_xlim(1, h.max())
ax.set_xlabel('grid refinement factor, h')
ax.set_ylabel('DE norms')
#save_fig(fig, 'MMS_DE_norms')


##### Observed order of accuracy phat

phat_p = find_phat(DEp)
phat_u = find_phat(DEu)
phat_v = find_phat(DEv)

fig, ax = plt.subplots()
ax.semilogx(h[1:], phat_p[:, 0], 'bo-', lw=2, label='L1 p')
ax.semilogx(h[1:], phat_p[:, 1], 'bs--', lw=2, label='L2 p')
ax.semilogx(h[1:], phat_p[:, 2], 'b*:', lw=2, label='Linf p')
ax.semilogx(h[1:], phat_u[:, 0], 'go-', lw=2, label='L1 u')
ax.semilogx(h[1:], phat_u[:, 1], 'gs--', lw=2, label='L2 u')
ax.semilogx(h[1:], phat_u[:, 2], 'g*:', lw=2, label='Linf u')
ax.semilogx(h[1:], phat_v[:, 0], 'ro-', lw=2, label='L1 v')
ax.semilogx(h[1:], phat_v[:, 1], 'rs--', lw=2, label='L2 v')
ax.semilogx(h[1:], phat_v[:, 2], 'r*:', lw=2, label='Linf v')
ax.legend(loc='lower right', frameon=False, fancybox=True, ncol=3,
          title='Driven cavity MMS: Re=10, \nCFL=0.9, finest mesh=257x257',
          prop={'size':10})
ax.set_xlim(1, h[1:].max())
ax.set_ylim(0, 4)
ax.set_xlabel('grid refinement factor, h')
ax.set_ylabel('observed order of accuracy, p')
#save_fig(fig, 'MMS_OOA')


def find_gci(d2, d1, p=2, Fs=3, r=2):
    """Grid Convergence Index using Roache (1994) method.

    Inputs:
        x2 (array) - grid of the coarse mesh
        x1 (array) - grid of the fine mesh
        f2 (array) - solution on coarse mesh
        f1 (array) - solution on fine meash
        p (float) - order of accuracy
        FS (default=3) - factor of safety
        r (default=2) - grid refinement factor
    
    Outputs:
        gci (float) - grid convergence index of fine grid
    """
    
    x2, f2 = d2[:, 0], d2[:, 1]
    x1, f1 = d1[:, 0], d1[:, 1]
    
    assert len(f1) > len(f2)
    
    # common grid points for both meshes
    f1_mask = np.in1d(x1, x2)
    
    c = np.float(Fs) / (r**p - 1)
    gci = c * np.abs((f2 - f1[f1_mask]) / f1[f1_mask])
    
    return x2, gci



##### Effect of C4 constant #################################################

c4_base = './out_MMS/cavity_C4={0}_SGS_129x129_Re=10_cfl=0.5.tec'
c4_vals = ['.008', '.01', '.04', '.0625']
vlevels = 1e3*np.linspace(-0.00025, 2.56e-05, 20)
vlevels = 1e3*np.linspace(-0.00031, 3.045e-05, 20)
vlevels = np.linspace(-2*.125, .15, 41)
#max, min = [], []  # use to get contour levels
cmap = plt.cm.PuOr_r
for val in c4_vals:
    n, d = load_tec(c4_base.format(val))
    
    plt.figure()
    plt.contourf(d[:, :, 0], d[:, :, 1], 1e3*(d[:, :, 8]/d[:, :, 5]), vlevels,
                 cmap=cmap)
    cbar = plt.colorbar()
    cbar.set_label(r'1000 $\times$ DE-p, N/m$^2$')
    plt.xlim(0, .006)
    plt.ylim(.050-.006, .050)
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    CS = plt.contour(d[:, :, 0], d[:, :, 1], 1e3*(d[:, :, 8]/d[:, :, 5]), levels=[0],
                color='k', lw=4)
    plt.clabel(CS, inline=1, fontsize=10)
    if False:
        plt.savefig(save_dir + 'MMS_C4_{0}.eps'.format(val.strip('.')), 
                    bbox_inches='tight', pad_inches=0)


##### Contour plots for the cavity at different Re ###########################

def plot_cavity(dd, Re):
    """Make contour plots of the cavity."""
    
    cmap = plt.cm.RdBu_r

    # U-velocity
    fig = plt.figure()
    cax = plt.contourf(dd[:, :, 0], dd[:, :, 1], dd[:, :, 3], 
                       np.r_[-.4:1.15:.05], cmap=cmap)
    cbar = fig.colorbar(cax, orientation='vertical')
    cbar.set_label('u, m/s')
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    if save_figures:
        plt.savefig(save_dir + 'cavityU_Re={0}.eps'.format(Re), 
                    bbox_inches='tight', pad_inches=0)
    
    # P
    plt.figure()
    cax = plt.contourf(dd[:, :, 0], dd[:, :, 1], dd[:, :, 2],
                       np.linspace(-3, 3, 41), cmap=cmap)
    cbar = plt.colorbar()
    cbar.set_label(r'p, N/m$^2$')
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    if Re == 100 and save_figures:
        plt.savefig(save_dir + 'cavityP_Re={0}.eps'.format(Re),
                    bbox_inches='tight', pad_inches=0)
    
    # V-velocity
    plt.figure()
    cax = plt.contourf(dd[:, :, 0], dd[:, :, 1], dd[:, :, 4],
                       np.r_[-.65:.525:.05], cmap=cmap)
    cbar = plt.colorbar()
    cbar.set_label('v, m/s')
    plt.xlabel('x, m')
    plt.ylabel('y, m')
    if save_figures:
        plt.savefig(save_dir + 'cavityV_Re={0}.eps'.format(Re),
        bbox_inches='tight', pad_inches=0)



re_files = ['cavity_SGS_257x257_Re=100_cfl=1.2',
            'cavity_SGS_257x257_Re=500_cfl=0.4',
            'cavity_SGS_257x257_Re=1000_cfl=0.4']
res = [100, 500, 1000]
for idx, f in enumerate(re_files):
    n, dd = load_tec('./out_cavity/' + f)
    plot_cavity(dd, res[idx])


## try to make a streamplot
#n, dd = load_tec('./out_cavity/' + re_files[0])
#x, y = dd[0, :, 0], dd[:, 0, 1]
#u, v = dd[:, :, 3], dd[:, :, 4]
#x = np.around(x, decimals=4)
#y = np.around(x, decimals=4)
#plt.figure()
#streamplot(x, y, u, v, density=5)
    

## Show the figures
plt.show()