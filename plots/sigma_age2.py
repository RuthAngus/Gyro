import numpy as np

import matplotlib.pyplot as pl
from teff_bv import teff2bv_orig
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

def distance(pars, a_obs, p_obs, bv_obs, a_err):
    model = age_model(pars, p_obs, bv_obs)
    return ((a_obs-model))**2/a_err

def distance2(age_obs, age_model, a_err):
    x = (age_obs-age_model)/a_err
    y = 1./(1+np.exp((x**2-.05)/1.))
    pl.clf()
    pl.plot(x, y, 'k.')
    pl.savefig('function')
    return 1./(1+np.exp((x**2-.05)/.5))

def distance3(age_obs, age_model, a_err):
    x = (age_obs-age_model)/a_err
    return 1./(1+np.exp((x**2-1)/10.))

def iso_calc(pars, age):
    x = np.linspace(1.6, .4, 10000)
    y = 10**log_period_model(pars, np.log10(age*1000), x)
    return x, y

data = np.genfromtxt("/Users/angusr/Python/Gyro/data/no_dup.txt").T

# make up bv_errs
c_err = .01

# remove subs
subgiant = 3.9
g = data[10] > subgiant

# p1 = data[1][g]
# p_err1 = data[2][g]
# a1 = data[3][g]
# a_errp1 = data[4][g]
# a_errm1 = data[5][g]
# a_err1 = .5 * (a_errp1 + a_errm1)
# t1 = data[12][g]
# t_err1 = data[13][g]
# logg1 = data[9][g]
# logg_err1 = data[10][g]
# bv1 = teff2bv(t1, logg1, np.ones_like(t1)*-.02, t_err1, logg_err1, \
#         np.ones_like(t1)*.001, error=False)

p1 = data[1][g]
p_err1 = data[2][g]
t1 = data[3][g]
t_err1 = data[4][g]
a1 = data[13][g]
a_errp1 = data[14][g]
a_errm1 = data[15][g]
a_err1 = .5 * (a_errp1 + a_errm1)
logg1 = data[10][g]
logg_err1 = data[11][g]
feh1 = data[5][g]
feh_err1 = data[6][g]
bv1 = teff2bv_orig(t1, logg1, feh1)
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

ages = [.65, 1, 2, 3, 4.57, 7, 10, 13]

pars = [.7725, .5189, .601, .4] # Barnes
pars_err = [.0070, .011, .024, 0.]

pars2 = [.407, .566, .325, .495] # MH
pars2_err = [0.008, 0.021, 0.024, 0.010]

pars4 = [0.63951838, 0.4811916, 0.48499184, .4] # .4
pars45 = [0.65713966, 0.47360685, 0.46631881, .45] # 45
pars5 = [0.41229508, 0.54640667, 0.46986845, .5]
# sun = [0.40117544, 0.55100903, 0.40028554, .45] # sun
pars3_err = np.array([.03, .03, .03, .00])

# params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameterssun.txt').T # sun fixed
# params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parametersgyro.txt').T
# params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters45_2.txt').T
params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters45.txt').T # sun free
# params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parametersgyrosun.txt').T
pars3 = np.zeros(4)
err = np.zeros((2, 4))
pars3[:3] = params[0][:3]
pars3[-1] = .45; err[:,3] = .0
err[0,:3] = params[1][:3]
err[1,:3] = params[2][:3]
pars3_err = np.zeros(4)
for i in range(4):
#     pars3_err[i] = .5*sum(err[:,i])
    pars3_err[i] = min(err[:,i])

for i, age in enumerate(ages):
#     sig = .5
    sig = 1.
    l11 = (a1-(a_errm1*sig) < age) * (age < a1+(a_errp1*sig)) # cool astero
    l12 = (a2-(a_errm2*sig) < age) * (age < a2+(a_errp2*sig)) # cool clusters
    l21 = (a1-(a_errm1*sig) < age) * (age < a1+(a_errp1*sig)) * (bv1 < .4) # hot astero
    l22 = (a2-(a_errm2*sig) < age) * (age < a2+(a_errp2*sig)) * (bv2 < .4) # hot clusters

    dist = distance2(a1, age, a_errp1)
    dist2 = distance2(a2, age, a_errp2)

    # Plot data
    pl.clf()
    cm = pl.cm.get_cmap('Greys')
    pl.errorbar(bv1[l11], p1[l11], xerr=bv_err1[l11], yerr=p_err1[l11], color='k', \
            fmt='o', mec='k', capsize=0, markersize=2, ecolor='.7', zorder=3)#, alpha=.05)
#             fmt='o', mec='k', capsize=0, markersize=6, ecolor='k', zorder=3)#, alpha=.05)
#     pl.scatter(bv1, p1, c=dist, cmap=cm, marker='o', s=30, zorder=2, \
#             edgecolors='None', alpha=.6)
#     pl.scatter(bv2, p2, c=dist2, cmap=cm, marker='^', s=40, zorder=1, \
#             edgecolors='None', alpha=.6)
    pl.errorbar(bv2[l12], p2[l12], xerr=bv_err2[l12], yerr=p_err2[l12], color='k', \
            fmt='^', mec='k', capsize=0, markersize=3, ecolor='0.7', zorder=0)#, alpha=.1)
#             fmt='^', mec='k', capsize=0, markersize=6, ecolor='0.4', zorder=0)#, alpha=.1)

    # Add Isochrones
    xs, ys = iso_calc(pars, ages[i])
#     pl.plot(xs, ys, color=ocols[i], linestyle='-.', linewidth=lw, \
    pl.plot(xs, ys, color='k', linestyle='-.', linewidth=lw, \
            label='$%s~\mathrm{Gyr}$~$\mathrm{(Barnes~2007)}$' %ages[i], zorder=0)
    xs, ys = iso_calc(pars2, ages[i])
#     pl.plot(xs, ys, color=ocols[i], linestyle='--', linewidth=lw, \
    pl.plot(xs, ys, color='k', linestyle='--', linewidth=lw, \
            label='$%s~\mathrm{Gyr}$~$\mathrm{(M\&H~2008)}$' %ages[i], zorder=0)
    xs, ys = iso_calc(pars3, ages[i])
#     pl.plot(xs, ys, color = ocols[i], linestyle='-', linewidth=2, \
    pl.plot(xs, ys, color='k', linestyle='-', linewidth=lw, \
            label = '$%s~\mathrm{Gyr}$~ \
            $\mathrm{Angus~\emph{et~al.}~(in~prep)}$' %ages[i], zorder=0)
    xs, ys1 = iso_calc(pars3-pars3_err, ages[i])
    xs, ys2 = iso_calc(pars3+pars3_err, ages[i])
#     pl.fill_between(xs, ys1, ys2, facecolor=ocols[i], alpha=0.3, edgecolor='None', \
    pl.fill_between(xs, ys1, ys2, facecolor='0.5', alpha=0.3, edgecolor='None', \
            zorder=0)

    pl.xlabel("$\mathrm{B-V}$")
    pl.ylabel("$\mathrm{P_{rot} (days)}$")
    pl.xlim(.2, 1.)
#     pl.xlim(.2, 1.52)
    pl.ylim(0, 60)
    pl.legend(loc='upper left')
    pl.savefig("p_vs_bv%s"%i)
