import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv
import pretty5

plotpar = {'axes.labelsize': 20, 'text.fontsize': 20,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#99CC99', '#9933FF', '#CC0000']

def log_period_model(par, log_a, bv):
    return np.log10(par[0]) + par[1] * log_a + par[2] * np.log10(bv - par[3]) # colour

def age_model(par, p, bv):
    return (1./par[0]) * p**(1./par[1]) * (bv-par[3])**(-par[2]/par[0])

def distance(pars, a_obs, bv, a_err):
    model = age_model(pars, a_obs, bv)
    return ((a_obs-model))**2/a_err

def iso_calc(pars, age):
    x = np.linspace(1.3, .4, 1000)
    y = 10**log_period_model(pars, np.log10(age*1000), x)
    return x, y

data = np.genfromtxt("/Users/angusr/Python/Gyro/data/matched_data.txt").T

# make up bv_errs
c_err = .01

# remove subs
g = data[9] > 4.

p1 = data[1][g]
p_err1 = data[2][g]
a1 = data[3][g]
a_errp1 = data[4][g]
a_errm1 = data[5][g]
a_err1 = .5 * (a_errp1 + a_errm1)
t1 = data[12][g]
t_err1 = data[13][g]
logg1 = data[9][g]
logg_err1 = data[10][g]
bv1 = teff2bv(t1, logg1, np.ones_like(t1)*-.02, t_err1, logg_err1, \
        np.ones_like(t1)*.001, error=False)
bv_err1 = np.ones_like(bv1)*c_err

# add clusters
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt").T
bv2 = data[0]
bv_err2 = data[1]
p2 = data[2]
p_err2 = data[3]
a2 = data[4]
a_err2 = data[5]
a_errp2 = data[5]
a_errm2 = data[5]
logg2 = data[6]

logg_err2 = data[7]

bv = np.concatenate((bv1, data[0]))
bv_err = np.concatenate((bv_err1, data[1]))
p = np.concatenate((p1, data[2]))
p_err = np.concatenate((p_err1, data[3]))
a = np.concatenate((a1, data[4]))
a_err = np.concatenate((a_err1, data[5]))
a_errp = np.concatenate((a_errp1, data[5]))
a_errm = np.concatenate((a_errm1, data[5]))
logg = np.concatenate((logg1, data[6]))
logg_err = np.concatenate((logg_err1, data[7]))

# ages = [.65, 1., 2., 3., 4., 5., 8., max(a)]
ages = [.65, 1., 3., 5., 7., 9., 11.,max(a)]
pars = [.7725, .5189, .601, .4] # Barnes
pars_err = [.0070, .011, .024, 0.]
pars2 = [.566, .407, .325, .495] # MH
pars2_err = [0.008, 0.021, 0.024, 0.010]
# pars3 = np.array([.6, .5, .601, .4])
# pars3 = [0.6459, 0.52115295, 0.52449312, .45]
pars3 = [0.63951838, 0.4811916, 0.48499184, .4]
pars3_err = np.array([.03, .03, .03, .00])

for i, age in enumerate(ages):
#     sig = .5
    sig = 1.
    l11 = (a1-(a_errm1*sig) < age) * (age < a1+(a_errp1*sig)) # cool astero
    l12 = (a2-(a_errm2*sig) < age) * (age < a2+(a_errp2*sig)) # cool clusters
    l21 = (a1-(a_errm1*sig) < age) * (age < a1+(a_errp1*sig)) * (bv1 < .4) # hot astero
    l22 = (a2-(a_errm2*sig) < age) * (age < a2+(a_errp2*sig)) * (bv2 < .4) # hot clusters

    dist = distance(pars3, a1, bv1, a_errp1)
    dist[np.isnan(dist)] = 0

    # Plot data
    pl.clf()
    cm = pl.cm.get_cmap('Greys')
    pl.errorbar(bv1[l11], p1[l11], xerr=bv_err1[l11], yerr=p_err1[l11], color='k', \
            fmt='o', mec='k', capsize=0, markersize=5, ecolor='0.8', zorder=4)
    pl.scatter(bv1[l11], p1[l11], c=distance(pars3, a1[l11], bv1[l11], a_errp1[l11]), \
            cmap=cm, marker='o', s=40, zorder=5, edgecolors='None')
    pl.errorbar(bv1[l21], p1[l21], xerr=bv_err1[l21], yerr=p_err1[l21], color='.75', \
            fmt='o', mec='.75', capsize=0, markersize=5, ecolor='0.8', zorder=4)
    pl.errorbar(bv2[l12], p2[l12], xerr=bv_err2[l12], yerr=p_err2[l12], color='k', \
            fmt='+', mec='k', capsize=0, markersize=5, ecolor='0.8', zorder=4)
    pl.errorbar(bv2[l22], p2[l22], xerr=bv_err2[l22], yerr=p_err2[l22], color='.75', \
            fmt='+', mec='.75', capsize=0, markersize=5, ecolor='0.8', zorder=4)
    pl.colorbar()

    # Add Isochrones
    xs, ys = iso_calc(pars, ages[i])
    pl.plot(xs, ys, color = ocols[i], linestyle='--', linewidth = 2, \
            label = '$%s~\mathrm{Gyr}$~ \
            $\mathrm{(Barnes~2007)}$' %ages[i])
    xs, ys = iso_calc(pars2, ages[i])
    pl.plot(xs, ys, color = ocols[i], linestyle='-.', linewidth = 2, \
            label = '$%s~\mathrm{Gyr}$~ \
            $\mathrm{(M\&H~2008)}$' %ages[i])
    xs, ys = iso_calc(pars3, ages[i])
    pl.plot(xs, ys, color = ocols[i], linestyle='-', linewidth = 2, \
            label = '$%s~\mathrm{Gyr}$~ \
            $\mathrm{Angus~\emph{et~al.}~(in~prep)}$' %ages[i])
    xs, ys1 = iso_calc(pars3-pars3_err, ages[i])
    xs, ys2 = iso_calc(pars3+pars3_err, ages[i])
    pl.fill_between(xs, ys1, ys2, facecolor=ocols[i], alpha=0.3, edgecolor='None')

    pl.xlabel("$\mathrm{B-V}$")
    pl.ylabel("$\mathrm{P_{rot} (days)}$")
    pl.xlim(.2, .8)
    pl.ylim(0, 70)
    pl.legend(loc='upper left')
    pl.savefig("p_vs_bv%s"%i)
