import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from load_hyades import hya_load
from teff_bv import teff2bv_orig

ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#99CC99', '#9933FF', '#CC0000']
plotpar = {'axes.labelsize': 15,
           'text.fontsize': 15,
           'legend.fontsize': 15,
           'xtick.labelsize': 16,
           'ytick.labelsize': 16,
           'text.usetex': True}
pl.rcParams.update(plotpar)

pars = [0.5189, 0.7725, 0.601, 0.4]
n, a, b, c = pars
a_surf = np.linspace(.65, 15, 100)*1000
c_surf = np.linspace(.4, 1.4, 100)
a_surf, c_surf = np.meshgrid(a_surf, c_surf)
p_surf = a_surf**n * a * (c_surf - c)**b
a_surf /= 1000.

# load hyades data
p, c, a, p_err, c_err, a_err = hya_load()
a /= 1000.

# load praesepe data
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/praesepe.txt').T
c = np.concatenate((c, (data[5]-data[6])))
c_err = np.concatenate((c_err, np.ones_like(data[5])*.1))
p = np.concatenate((p, 1./data[3]))
p_err = np.concatenate((p_err, (1./data[3])*(data[4]/data[3])))
a = np.concatenate((a, np.ones_like(data[5])*.588))
a_err = np.concatenate((a_err, np.ones_like(data[5])*.137))

# load my data
# data = np.genfromtxt("/Users/angusr/Python/Gyro/data/data.txt").T
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/new_matched.txt").T
x = data[10] > 1.
period = data[1][x]
teff = data[3][x]
feh = data[5][x]
logg = data[10][x]
age = data[13][x]
bv = teff2bv_orig(teff, logg, feh)

pl.clf()
fig = pl.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(c_surf, a_surf, p_surf, alpha = 0.0)
ax.scatter(c, a, p, color='k')
ax.scatter(bv, age, period, color=ocols[0])
ax.set_zlabel('$P_{rot}~\mathrm{(days)}$')
ax.set_xlabel('$\mathrm{B-V}$')
ax.set_ylabel('$\mathrm{Age~(Gyr)}$')
pl.show()
pl.savefig('3d_meibom')
