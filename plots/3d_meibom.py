import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from load_hyades import hya_load
from teff_bv import teff2bv_orig

ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#99CC99', '#9933FF', '#CC0000']
plotpar = {'axes.labelsize': 30,
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

# load my data
# data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_all_astero.txt")
# data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_all_astero_no_precise.txt")
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/all_astero_plusgarcia.txt")
t = data[1]
a = data[3]
p = data[6]
logg = data[8]
feh = data[11]
bv = teff2bv_orig(t, logg, feh)

# load precise data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/vandt.txt", skip_header=1).T
pt = data[1]
pa = data[3]
pp = data[6]
plogg = data[8]
pfeh = data[11]
pbv = teff2bv_orig(pt, plogg, pfeh)

# load cluster and field stars data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt").T
cbv = data[0]
cp = data[2]
ca = data[4]

pl.clf()
fig = pl.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(c_surf, a_surf, p_surf, alpha = 0.0)
ax.scatter(bv, a, p, color='k', zorder=1)
# ax.scatter(pbv, pa, pp, color=ocols[2], zorder=3)
s1 = ax.scatter(cbv[:-5], ca[:-5], cp[:-5], color='r', zorder=2)
s1.set_edgecolors = s1.set_facecolors = lambda *args:None
s2 = ax.scatter(cbv[-5:-1], ca[-5:-1], cp[-5:-1], color='k', zorder=2)
s2.set_edgecolors = s2.set_facecolors = lambda *args:None
s3 = ax.scatter(cbv[-1], ca[-1], cp[-1], color='k', zorder=2)
s3.set_edgecolors = s3.set_facecolors = lambda *args:None

ax.set_zlabel('$P_{rot}~\mathrm{(days)}$')
ax.set_xlabel('$\mathrm{B-V}$')
ax.set_ylabel('$\mathrm{Age~(Gyr)}$')
pl.show()
pl.savefig('3d_meibom')
