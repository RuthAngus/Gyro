import numpy as np

import matplotlib.pyplot as pl
from teff_bv import teff2bv_orig, teff2bv_err
import pretty5

plotpar = {'axes.labelsize': 20, 'text.fontsize': 20,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#9933FF', '#CC0000', '#9933FF', '#99cc99', '#CC0000']

lw = 1

def log_period_model(par, log_a, bv):
    return np.log10(par[0]) + par[1] * log_a + par[2] * np.log10(bv - par[3]) # colour

def age_model(par, p, bv):
    return (1./par[0]) * p**(1./par[1]) * (bv-par[3])**(-par[2]/par[0])

def iso_calc(pars, age):
    x = np.linspace(1.6, .4, 10000)
    y = 10**log_period_model(pars, np.log10(age*1000), x)
    return x, y

# data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_all_astero.txt")
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_irfm.txt")

# make up bv_errs
c_err = .04

# remove subs
subgiant = 4.2
g = data[8] > subgiant

p1 = data[6][g]
p_err1 = data[7][g]
t1 = data[1][g]
t_err1 = data[2][g]
a1 = data[3][g]
a_errp1 = data[4][g]
a_errm1 = data[5][g]
a_err1 = .5 * (a_errp1 + a_errm1)
logg1 = data[8][g]
logg_errp1 = data[9][g]
logg_errm1 = data[10][g]
logg_err1 = .5*(logg_errp1+logg_errm1)
feh1 = data[11][g]
feh_err1 = data[12][g]
bv1, bv_err1 = teff2bv_err(t1, logg1, feh1, t_err1, logg_err1, feh_err1)

# add clusters
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt").T
aaa = data[4]
l = (aaa!=1.1) * (aaa!=.588)
bv2 = data[0][l]
bv_err2 = data[1][l]
p2 = data[2][l]
p_err2 = data[3][l]
a2 = data[4][l]
a_err2 = data[5][l]
a_errp2 = data[5][l]
a_errm2 = data[5][l]
logg2 = data[6][l]
logg_err2 = data[7][l]

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

ages = [.625, 2, 4.568, 8, 10]
pars = [.7725, .5189, .601, .4] # Barnes
pars_err = [.0070, .011, .024, 0.]
pars2 = [.407, .566, .325, .495] # MH
pars2_err = [0.008, 0.021, 0.024, 0.010]
pars3_err = np.array([.03, .03, .03, .00])
params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parametersACHF45irfm2.txt').T
pars3 = np.zeros(4)
err = np.zeros((2, 4))
pars3[:3] = params[0][:3]
pars3[-1] = .45; err[:,3] = .0
err[0,:3] = params[1][:3]
err[1,:3] = params[2][:3]
pars3_err = np.zeros(4)
for i in range(4):
    pars3_err[i] = .5*sum(err[:,i])

for i, age in enumerate(ages):
    sig = 1.
    l11 = (a1-(a_errm1*sig) < age) * (age < a1+(a_errp1*sig)) # cool astero
    l12 = (a2-(a_errm2*sig) < age) * (age < a2+(a_errp2*sig)) # cool clusters
    l21 = (a1-(a_errm1*sig) < age) * (age < a1+(a_errp1*sig)) * (bv1 < pars3[-1]) # hot astero
    l22 = (a2-(a_errm2*sig) < age) * (age < a2+(a_errp2*sig)) * (bv2 < pars3[-1]) # hot clusters
    sun = a2==4.568

    # Plot data
    pl.clf()
    pl.errorbar(bv1[l11], p1[l11], xerr=bv_err1[l11], yerr=p_err1[l11], color='k', \
            fmt='o', mec='k', capsize=0, markersize=3, ecolor='.7', zorder=3)
    if age == 4.568:
        pl.errorbar(bv2[sun], p2[sun], xerr=bv_err2[sun], yerr=p_err2[sun], color='r', \
                fmt='o', mec='r', capsize=0, markersize=4, ecolor='.7', zorder=3)
    pl.errorbar(bv2[l12], p2[l12], xerr=bv_err2[l12], yerr=p_err2[l12], color='r', \
            fmt='.', mec='r', capsize=0, markersize=8, ecolor='0.7', zorder=0)

    # Add Isochrones
    xs, ys = iso_calc(pars, ages[i])
    pl.plot(xs, ys, color='k', linestyle='-.', linewidth=lw, \
            label='$%s~\mathrm{Gyr}$~$\mathrm{Barnes~(2007)}$' %ages[i], zorder=0)
    xs, ys = iso_calc(pars2, ages[i])
    pl.plot(xs, ys, color='k', linestyle='--', linewidth=lw, \
            label='$%s~\mathrm{Gyr}$~$\mathrm{M\&H~(2008)}$' %ages[i], zorder=0)
    xs, ys = iso_calc(pars3, ages[i])
    pl.plot(xs, ys, color='k', linestyle='-', linewidth=lw, \
            label = '$%s~\mathrm{Gyr}$~ \
            $\mathrm{Angus~\emph{et~al.}~(2014)}$' %ages[i], zorder=0)
    xs, ys1 = iso_calc(pars3-pars3_err, ages[i])
    xs, ys2 = iso_calc(pars3+pars3_err, ages[i])
    pl.fill_between(xs, ys1, ys2, facecolor='0.5', alpha=0.3, edgecolor='None', \
            zorder=0)

    pl.xlabel("$\mathrm{B-V}$")
    pl.ylabel("$\mathrm{P_{rot} (days)}$")
    pl.xlim(.2, 1.)
    pl.ylim(0, 50)
    pl.legend(loc='upper left')
    pl.savefig("p_vs_bv%s"%i)
