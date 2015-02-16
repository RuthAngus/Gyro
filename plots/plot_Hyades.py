import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv_orig, teff2bv_err
import pretty5
import sys

plotpar = {'axes.labelsize': 18, 'text.fontsize': 10,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#9933FF', '#CC0000', '#9933FF', '#99cc99', '#CC0000']

ms = 8
lw = 1
c = ocols[0]
c = 'k'

def log_period_model(par, log_a, bv):
    return np.log10(par[0]) + par[1] * log_a + par[2] * np.log10(bv - par[3]) # colour

def age_model(par, p, bv):
    return ((p/(par[0]*(bv-par[3])**(par[2])))**(1./par[1]))/1000.

def log_age_model(par, log_p, bv):
    return (log_p - np.log10(par[0]) - par[2]*np.log10(bv - par[3])) / par[1]

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
    x = np.linspace(1.8, .1, 10000)
    y = 10**log_period_model(pars, np.log10(age*1000), x)
    return x, y

data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_irfm.txt")

# remove subs
subgiant = 4.1
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
logg_err1 = .5*(data[9][g]+data[10][g])
feh1 = data[11][g]
feh_err1 = data[12][g]
flag1 = data[13][g]
bv1, bv_err1 = teff2bv_err(t1, logg1, feh1, t_err1, logg_err1, feh_err1)

# add clusters
data = np.genfromtxt("/Users/angusr/Python/Gyro/code/clusters.txt").T
l = (data[4] != 1.1) * (data[4] != .588)
bv2 = data[0][l]; bv_err2 = data[1][l]
p2 = data[2][l]; p_err2 = data[3][l]
a2 = data[4][l]; a_err2 = data[5][l]; a_errp2 = data[5][l]; a_errm2 = data[5][l]
logg2 = data[6][l]; logg_err2 = data[7][l]
flag2 = data[8][l]

# combine astero and cluster
bv = np.concatenate((bv1, bv2))
bv_err = np.concatenate((bv_err1, bv_err2))
p = np.concatenate((p1, p2))
p_err = np.concatenate((p_err1, p_err2))
a = np.concatenate((a1, a2))
a_err = np.concatenate((a_err1, a_err2))
a_errp = np.concatenate((a_errp1, a_errp2))
a_errm = np.concatenate((a_errm1, a_errm2))
logg = np.concatenate((logg1, logg2))
logg_err = np.concatenate((logg_err1, logg_err2))
flag = np.concatenate((flag1, flag2))

age = 4.568
pars = [.7725, .5189, .601, .4] # Barnes
pars_err = [.0070, .011, .024, 0.]

pars2 = [.407, .566, .325, .495] # MH
pars2_err = [0.008, 0.021, 0.024, 0.010]

# fnames
fnames = ['A', 'H', 'P', 'N', 'C', 'F', 'V']

fname = "HF45"

ck = 0

npars = 3
if ck >= 0:
    npars = 4

params = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters%s.txt'%fname).T
pars3 = np.zeros(npars)
err = np.zeros((2, npars))

pars3[:npars-1] = params[0][:npars-1]

# set position of colour singularity
print fname
if ck <0 :
    if fname.find('4') >= 0:
        pars3[-1] = .4; err[:,3] = .0
    elif fname.find('5') >= 0:
        pars3[-1] = .5; err[:,3] = .0
    if fname.find('5') >= 0 and fname.find('4') >= 0:
        pars3[-1] = .45; err[:,3] = .0
    if fname.rfind('5') != fname.find('5'):
        pars3[-1] = .55; err[:,3] = .0
pars3[-1] = .45

err[0,:npars] = params[1][:npars]
err[1,:npars] = params[2][:npars]
pars3_err = np.zeros(npars)
for i in range(npars):
    pars3_err[i] = .5*sum(err[:,i])
pars3_err[-1] = 0.

# turn astero target flags to '3' and take 3 away
flag[flag==2] = 9
flag -= 3
flag[flag<0] = 0

# select star group
flist = []
for i in range(len(fnames)):
    if fname.find(fnames[i]) >= 0:
        print fnames[i], 'yes'
        flist.append(i)

l = (np.sum([flag == i for i in flist], axis=0)) == 1

# ages
ages = [age, .625, 1.1, .588, .5, age, age]

# field star data
DIR = '/Users/angusr/Python/Gyro/data'
scobv, ebv, scop, ep, scoa, ea = np.genfromtxt('%s/18sco.txt'%DIR, skip_header=2).T
cygt, et, cygp, ep, cyga, ea, cygfeh, efeh, cygbv, ebv = \
        np.genfromtxt('%s/16Cygb.txt'%DIR, skip_header=2).T
acbv, ebv, acp, ep, aca, ea = np.genfromtxt('%s/alphacen.txt'%DIR, skip_header=2).T

pl.clf()
pl.subplot(2, 1, 1)
bv -= 0.45
pl.errorbar(bv[l], p[l], xerr=bv_err[l], yerr=p_err[l], fmt='k.', capsize=0,
           ecolor='.7', markersize=6)

lw = 1
styles = ['-', '-', '-', '-', '-', '--', '-', '-', '-']
# plot isochrones
for i in range(len(fnames)):
    if fname.find(fnames[i]) >= 0:
        print fnames[i], styles[i]
        xs, ys = iso_calc(pars3, ages[i])
        xs -= .45
        pl.plot(xs, ys, color = c, linestyle=styles[i], linewidth=lw, \
                label = '$%s~\mathrm{Gyr}$' %ages[i], zorder=0)
pl.legend(loc='upper left')
# pl.text(1.45, 2, '$%s~\mathrm{Gyr}$'%ages[1])

pl.xlabel('$\mathrm{B-V-}~c$')
# pl.ylim(0,70)
pl.ylabel('$\mathrm{Period~(days)}$')
# pl.xlim(.2,1.8)
pl.subplots_adjust(hspace=0.001)
pl.plot(.65-.45, 26.09, 'ro', markersize=ms, mec='r')
# pl.xlim(.2, 1.8)
pl.xlim(10**-1.6, 10**0.1)
pl.ylim(10**.6, 10**1.8)
pl.loglog()

pl.subplot(2, 1, 2)
pl.errorbar(a[l], p[l], xerr=(a_errp[l], a_errm[l]), yerr=p_err[l], fmt='k.',
           capsize=0, ecolor='.7', markersize=ms)
xs = np.linspace(0, 20, 1000)
pl.plot(xs, 10**log_period_model(pars3, np.log10(xs*1000), .65), 'k',
        label = '$\mathrm{B-V}$=0.65')
pl.plot(4.568, 26.09, 'ro', markersize=ms, mec='r')
# pl.text(1.2, 50, '$4.568~\mathrm{Gyr}$')
pl.xlabel('$\mathrm{Age~(Gyr)}$')
pl.ylabel('$\mathrm{Period~(days)}$')
# pl.ylim(0,70)
# pl.xlim(-5, 20)
pl.legend(loc='upper left')
pl.xlim(3e-1,20)
pl.ylim(10**.5, 10**1.9)
pl.loglog()
pl.savefig("/Users/angusr/Python/Gyro/gyro_paper/show%s.pdf"%fname)
