import numpy as np

import matplotlib.pyplot as pl
from teff_bv import teff2bv_orig, teff2bv_err
import pretty5
import sys
from colours import plot_colours
cols = plot_colours()

plotpar = {'axes.labelsize': 18, 'text.fontsize': 10,
           'legend.fontsize': 12,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#9933FF', \
        '#CC0000', '#9933FF', '#99cc99', '#CC0000']

ms = 8
lw = 1
c = ocols[0]
c = 'k'

def log_period_model(par, log_a, bv):
    return np.log10(par[0]) + par[1] * log_a + par[2] * np.log10(bv - par[3]) # colour

def iso_calc(pars, age):
    x = np.linspace(1.8, .1, 10000)
    y = 10**log_period_model(pars, np.log10(age*1000), x)
    return x, y

data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_irfm.txt")

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
logg_err1 = .5*(data[9][g]+data[10][g])
feh1 = data[11][g]
feh_err1 = data[12][g]
flag1 = data[13][g]
bv1, bv_err1 = teff2bv_err(t1, logg1, feh1, t_err1, logg_err1, feh_err1)

# add clusters
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt").T
l = (data[4] != 1.1) * (data[4] != .588)
bv2 = data[0][l]; bv_err2 = data[1][l]
p2 = data[2][l]; p_err2 = data[3][l]
a2 = data[4][l]; a_err2 = data[5][l]; a_errp2 = data[5][l]
a_errm2 = data[5][l]
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

fname = sys.argv[1]
# fname = "ACHF45_2"

ck = 0

npars = 3
if ck >= 0:
    npars = 4

params = np.genfromtxt('/Users/angusr/Python/Gyro/code/parameters%s.txt'
                       % fname).T
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

sig = 1
per = .2
b = age*per
# cool stars
# l = (a-(a_errm*sig) < age) * (age < a+(a_errp*sig))
l = (age-b < a) * (a < age+b)
# hot stars
# l2 = (a-(a_errm*sig) < age) * (age < a+(a_errp*sig)) * (bv<.45)
l2 = (age-b < a) * (a < age+b) * (bv<.45)

# Plot data
pl.clf()
sun = a==4.568
bv -= 0.45
pl.errorbar(bv[l], p[l], xerr=bv_err[l], yerr=p_err[l], color='k', \
        fmt='o', mec='k', capsize=0, markersize=5, ecolor='.6')
pl.errorbar(bv[l2], p[l2], xerr=bv_err[l2], yerr=p_err[l2], color='k', \
        fmt='o', mec='.7', capsize=0, markersize=5, ecolor='.6')
pl.errorbar(bv[sun], p[sun], xerr=bv_err[sun], yerr=p_err[sun], color='r', \
        fmt='o', mec='r', capsize=0, markersize=10, ecolor='.6')

from prob_uncerts import sigmas
samples = sigmas(fname)
for i in range(np.shape(samples)[1]):
    par_samp = np.zeros(4)
    par_samp[:3] = samples[:, i]
    par_samp[3] = 0.45
    xs, ys = iso_calc(par_samp, age)
    xs -= 0.45
    pl.plot(xs, ys, color=c, zorder=0, alpha=.1)

# Add Isochrones
xs, ys = iso_calc(pars, age)
xs -= 0.4
c = "k"
pl.plot(xs, ys, color=c, linestyle='-.', linewidth=lw, \
        label='$%s~\mathrm{Gyr}$~$\mathrm{(Barnes~2007)}$' %age, zorder=0)
xs, ys = iso_calc(pars2, age)
xs -= pars2[-1]
pl.plot(xs, ys, color=c, linestyle='--', linewidth=lw, \
        label='$%s~\mathrm{Gyr}$~$\mathrm{(M\&H~2008)}$' %age, zorder=0)
xs, ys = iso_calc(pars3, age)
xs -= 0.45
pl.plot(xs, ys, color = c, linestyle='-', linewidth=lw, \
        label = '$%s~\mathrm{Gyr}$~ \
        $\mathrm{Angus~\emph{et~al.}~(2014)}$' %age, zorder=0)
xs, ys1 = iso_calc(pars3-pars3_err, age)
xs -= 0.45
xs, ys2 = iso_calc(pars3+pars3_err, age)
xs -= 0.45
# pl.fill_between(xs, ys1, ys2, facecolor=c, alpha=0.1, edgecolor='None', \
#         zorder=0)

pl.xlabel("$\mathrm{B-V-}~c$")
# pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$\mathrm{P_{rot} (days)}$")
# pl.xlim(.2, 1.)
pl.xlim(10**-3.5, 1.)
pl.ylim(10**.6, 10**1.8)
# pl.ylim(0, 60)
pl.loglog()
pl.legend(loc='upper left')
pl.savefig("/Users/angusr/Python/Gyro/gyro_paper/p_vs_bv_solar.pdf")
bv += 0.45
colour = .65
b = .1
sig = 2
l = (bv-(bv_err*sig) < colour) * (colour < bv+(bv_err*sig))
l = (colour-b < bv) * (bv < colour+b)
# l = (5778-300 < t1) * (t1 < 5778+300) # test colour conversion
sun = a==4.568
# print l[-5:]

lw = .5
# Plot data
pl.clf()
pl.errorbar(a[l], p[l], xerr=a_err[l], yerr=p_err[l], color='k', \
        fmt='o', mec='k', capsize=0, markersize=5, ecolor='.6', zorder=4)
pl.errorbar(a[-5:][l[-5:]], p[-5:][l[-5:]], xerr=a_err[-5:][l[-5:]],
            yerr=p_err[-5:][l[-5:]], color="b", fmt='o', mec="b",
            capsize=0, markersize=7, ecolor='.6', zorder=4)
pl.errorbar(a[sun], p[sun], xerr=a_err[sun], yerr=p_err[sun], color="r", \
        fmt='o', mec=cols.pink, capsize=0, markersize=7, ecolor='.6', zorder=4)
xs = np.linspace(0, 20, 100)
pl.plot(xs, 10**log_period_model(pars, np.log10(xs*1000), .65), 'k-.',
        label='$\mathrm{B-V}=0.65~\mathrm{Barnes~(2007)}$', zorder=3,
        linewidth=2, color=cols.lightblue)
pl.plot(xs, 10**log_period_model(pars2, np.log10(xs*1000), .65), 'k--',
        label = '$\mathrm{B-V}=0.65~\mathrm{M\&H~(2008)}$', zorder=3,
        linewidth=2, color=cols.lightblue)

pl.plot(xs, 10**log_period_model(pars3, np.log10(xs*1000), .65), 'k',\
        label = '$\mathrm{B-V}=0.65~\mathrm{Angus~\emph{et~al.}~(2014)}$',
        linewidth=1, zorder=1)
ys1 = 10**log_period_model(pars3+pars3_err, np.log10(xs*1000), .65)
ys2 = 10**log_period_model(pars3-pars3_err, np.log10(xs*1000), .65)

from prob_uncerts import sigmas
samples = sigmas(fname)
for i in range(np.shape(samples)[1]):
    par_samp = np.zeros(4)
    par_samp[:3] = samples[:, i]
    par_samp[3] = 0.45
    pl.plot(xs, 10**log_period_model(par_samp, np.log10(xs*1000), .65), 'k',
            alpha=.05, zorder=0)

pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.ylabel("$\mathrm{P_{rot}~(days)}$")
pl.xlim(3e-1,20)
pl.ylim(10**.6, 10**1.8)
pl.legend(loc='upper left')
pl.loglog()
pl.savefig("/Users/angusr/Python/Gyro/gyro_paper/p_vs_a_solar.pdf")

# ages
ages = [age, .625, 1.1, .588, .5, age, age]
pl.clf()
pl.subplot(2, 1, 1)
# bv -= 0.45
pl.errorbar(bv[l], p[l], xerr=bv_err[l], yerr=p_err[l], fmt='k.', capsize=0,
           ecolor='.7', markersize=6)

lw = 1
styles = ['-', '-', '-', '-', '-', '--', '-', '-', '-']
# plot isochrones
for i in range(len(fnames)):
    if fname.find(fnames[i]) >= 0:
        print fnames[i], styles[i]
        xs, ys = iso_calc(pars3, ages[i])
#         xs -= .45
        pl.plot(xs, ys, color = c, linestyle=styles[i], linewidth=lw, \
                label = '$%s~\mathrm{Gyr}$' %ages[i], zorder=0)
# pl.legend(loc='upper left')

pl.xlabel('$\mathrm{B-V-}~c$')
pl.ylim(0,70)
pl.ylabel('$\mathrm{Period~(days)}$')
pl.subplots_adjust(hspace=.3)
# pl.plot(.65-.45, 26.09, 'ro', markersize=6, mec='r')
pl.plot(.65, 26.09, 'ro', markersize=6, mec='r')
pl.xlim(.2, 1.8)
# pl.xlim(10**-1.6, 10**0.1)
# pl.ylim(10**.6, 10**1.8)
# pl.loglog()

pl.subplot(2, 1, 2)
pl.errorbar(a[l], p[l], xerr=(a_errp[l], a_errm[l]), yerr=p_err[l], fmt='k.',
           capsize=0, ecolor='.7', markersize=ms)
xs = np.linspace(0, 20, 1000)
pl.plot(xs, 10**log_period_model(pars3, np.log10(xs*1000), .65), 'k',
        label = '$\mathrm{B-V}$=0.65')
pl.plot(4.568, 26.09, 'ro', markersize=ms, mec='r')
pl.xlabel('$\mathrm{Age~(Gyr)}$')
pl.ylabel('$\mathrm{Period~(days)}$')
pl.ylim(0,70)
pl.xlim(-5, 20)
# pl.legend(loc='upper left')
# pl.xlim(3e-1,20)
# pl.ylim(10**.5, 10**1.9)
# pl.loglog()
pl.savefig("/Users/angusr/Python/Gyro/gyro_paper/show%s.pdf" % fname)
