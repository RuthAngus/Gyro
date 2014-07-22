import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv_orig, bv2teff, teff2bv_err
from pretty5 import plotting

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)

# ocols = ['#FF9933','#66CCCC','#FF33CC','#3399FF','#CC0066','#99CC99','#9933FF','#CC0000']
ocols = ['#FF9933','#009999','#FF33CC','#3399FF','#CC0066','#99CC99','#9933FF','#CC0000']

colour=True
kbreak = .45

# Mamajek and Hillenbrand
amh = 0.407 #pm 0.021
bmh = 0.325 #pm 0.024
cmh = 0.495 #pm 0.010
nmh = 0.566 #pm 0.008

# Barnes
amh = 0.7725 #pm 0.021
bmh = 0.601 #pm 0.024
cmh = 0.4 #pm 0.010
nmh = 0.5189 #pm 0.008

lab = "$\mathrm{M\&H~(2008)}$"
lab = "$\mathrm{Barnes~(2007)}$"

def iso(age):
    if colour==True:
        xp = np.linspace(.41, 1.2, 100)
        yp = age**nmh * amh * (xp-cmh)**bmh
    else:
        xp = np.linspace(6249, 3000, 100)
        yp = age**nmh * amh * (6250-xp)**bmh
    return xp, yp

data = np.genfromtxt('/Users/angusr/Python/Gyro/data/garcia_all_astero.txt')
KID = data[0]
teff = data[1]
age = data[3]
age_errp = data[4]
age_errm = data[5]
teff_err = data[2]
period = data[6]
period_err = data[7]
age_err = .5*(age_errp+age_errm)
logg = data[8]
logg_errp = data[9]
logg_errm = data[10]
logg_err = .5*(logg_errp+logg_errm)
feh = data[11]
feh_err = data[12]
bv, bv_err = teff2bv_err(teff, logg, feh, teff_err, logg_err, feh_err)

# cool dwarfs
pl.clf()
pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], mec=ocols[1], \
    ecolor='.7', capsize=0, fmt='.', markersize=8)
xp, yp = iso(1000.)
pl.plot(xp, yp, linestyle='--', color='.75', label=lab)
xp, yp = iso(2000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(5000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(10000.)
pl.plot(xp, yp, linestyle='--', color='.75')
pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(.2,1.)
pl.legend(loc="upper left")
pl.ylim(0,70)
pl.savefig("cool_stars1")

# load travis and victor data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/vandt.txt", skip_header=1).T
vtKID = data[0]
vtt = data[1]
vtt_err = data[2]
vta = data[3]
vta_errp = data[4]
vta_errm = data[5]
vta_err = .5*(vta_errp+vta_errm)
vtperiod = data[6]
vtperiod_err = data[7]
vtlogg = data[8]
vtlogg_errp = data[9]
vtlogg_errm = data[10]
vtlogg_err = .5*(vtlogg_errp+vtlogg_errm)
vtfeh = data[11]
vtfeh_err = data[12]
vtbv, vtbv_err = teff2bv_err(vtt, vtlogg, vtfeh, vtt_err, vtlogg_err, vtfeh_err)

# cool dwarfs with precise
pl.clf()
pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], mec=ocols[1], \
    ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(vtbv, vtperiod, yerr=vtperiod_err, xerr=vtbv_err, color=ocols[2],\
        mec=ocols[2], ecolor='.7', capsize=0, fmt='.', markersize=8)
xp, yp = iso(1000.)
pl.plot(xp, yp, linestyle='--', color='.75', label=lab)
xp, yp = iso(2000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(5000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(10000.)
pl.plot(xp, yp, linestyle='--', color='.75')
pl.legend(loc="upper left")
pl.ylim(0,70)
pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(.2, 1.)
pl.savefig("cool_stars11")

# p vs bv paper
# load all astero data
data = np.genfromtxt('../data/garcia_all_astero.txt')
ta = data[1]
ta_err = data[2]
pa = data[6]
pa_err = data[7]
logga = data[8]
logga_err = .5*(data[9]+data[10])
feha = data[11]
feha_err = data[12]
bva, bva_err = teff2bv_err(ta, logga, feha, ta_err, logga_err, feha_err)
pl.clf()
pl.errorbar(bva, pa, yerr=pa_err, xerr=bva_err, color='k', mec='k', \
    ecolor='.7', capsize=0, fmt='.', markersize=4)
pl.ylim(0,70)
pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(.2, 1.)
pl.savefig("p_vs_bv_paper")
# p vs t paper
pl.clf()
pl.errorbar(ta, pa, yerr=pa_err, xerr=ta_err, color='k', mec='k', \
    ecolor='.7', capsize=0, fmt='.', markersize=4)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
pl.ylim(0,70)
pl.savefig("p_vs_t_paper")

# load clusters
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/clusters.txt').T
cbv = data[0]
cbv_err = data[1]
cp = data[2]
cp_err = data[3]
ca = data[4]
ca_err = data[5]
print cbv

# cool stars with precise and cluster
pl.clf()
pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], mec=ocols[1], \
        ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(vtbv, vtperiod, yerr=vtperiod_err, xerr=vtbv_err, color=ocols[2],\
        mec=ocols[2], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[:-5], cp[:-5], yerr=cp_err[:-5], xerr=cbv_err[:-5], color=ocols[3], \
        mec=ocols[3], ecolor='.7', capsize=0, fmt='.', markersize=8)
xp, yp = iso(1000.)
pl.plot(xp, yp, linestyle='--', color='.75', label=lab)
xp, yp = iso(2000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(5000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(10000.)
pl.plot(xp, yp, linestyle='--', color='.75')
pl.legend(loc="upper left")
pl.ylim(0,70)
pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(.2, 1.)
pl.savefig("cool_stars12")

# cool stars with precise, cluster and field
pl.clf()
pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], mec=ocols[1], \
        ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(vtbv, vtperiod, yerr=vtperiod_err, xerr=vtbv_err, color=ocols[2],\
        mec=ocols[2], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[:-5], cp[:-5], yerr=cp_err[:-5], xerr=cbv_err[:-5], color=ocols[3], \
        mec=ocols[3], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[-5:], cp[-5:], yerr=cp_err[-5:], xerr=cbv_err[-5:], color=ocols[4], \
        mec=ocols[4], ecolor='.7', capsize=0, fmt='.', markersize=8)
xp, yp = iso(1000.)
pl.plot(xp, yp, linestyle='--', color='.75', label=lab)
xp, yp = iso(2000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(5000.)
pl.plot(xp, yp, linestyle='--', color='.75')
xp, yp = iso(10000.)
pl.plot(xp, yp, linestyle='--', color='.75')
pl.legend(loc="upper left")
pl.ylim(0,70)
pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(.2, 1.)
pl.savefig("cool_stars13")

# hot dwarfs
l = bv < kbreak
p, t, p_err, t_err, b, b_err = period[l], teff[l], period_err[l], teff_err[l], \
        bv[l], bv_err[l]
pl.clf()
pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], fmt='.', \
        mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~stars}$", \
        markersize=8)
pl.errorbar(t, p, yerr=p_err, xerr=t_err, color=ocols[0], mec=ocols[0], fmt='.', \
        ecolor='.7', capsize=0, label="$\mathrm{Hot~stars}$", markersize=8)
pl.errorbar(vtbv, vtperiod, yerr=vtperiod_err, xerr=vtbv_err, color=ocols[2],\
        mec=ocols[2], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[:-5], cp[:-5], yerr=cp_err[:-5], xerr=cbv_err[:-5], color=ocols[3], \
        mec=ocols[3], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[-5:], cp[-5:], yerr=cp_err[-5:], xerr=cbv_err[-5:], color=ocols[4], \
        mec=ocols[4], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.legend(loc="upper left")
pl.savefig("cool_stars2")

# subs
l2 = logg < 4.
p2, t2, p_err2, t_err2, b2, b_err2 = period[l2], teff[l2], period_err[l2], teff_err[l2], \
        bv[l2], bv_err[l2]
pl.clf()
pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], fmt='.', \
        mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~dwarfs}$", \
        markersize=8)
pl.errorbar(t, p, yerr=p_err, xerr=t_err, color=ocols[0], mec=ocols[0], fmt='.', \
        ecolor='.7', capsize=0, label="$\mathrm{Hot~dwarfs}$", markersize=8)

pl.errorbar(t2, p2, yerr=p_err2, xerr=t_err2, color=ocols[2], mec=ocols[2], fmt='.', \
        ecolor='.7', capsize=0, label="$\mathrm{Subgiants}$", markersize=8)
pl.errorbar(vtbv, vtperiod, yerr=vtperiod_err, xerr=vtbv_err, color=ocols[2],\
        mec=ocols[2], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[:-5], cp[:-5], yerr=cp_err[:-5], xerr=cbv_err[:-5], color=ocols[3], \
        mec=ocols[3], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[-5:], cp[-5:], yerr=cp_err[-5:], xerr=cbv_err[-5:], color=ocols[4], \
        mec=ocols[4], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.legend(loc="upper left")
pl.savefig("cool_stars3")

# with isochrones
p2, t2, p_err2, t_err2, b2, b_err2 = period[l2], teff[l2], period_err[l2], teff_err[l2], \
        bv[l2], bv_err[l2]
pl.clf()
pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], fmt='.', \
        mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~dwarfs}$", \
        markersize=8)
pl.errorbar(t, p, yerr=p_err, xerr=t_err, color=ocols[0], mec=ocols[0], fmt='.', \
        ecolor='.7', capsize=0, label="$\mathrm{Hot~dwarfs}$", markersize=8)
pl.errorbar(t2, p2, yerr=p_err2, xerr=t_err2, color=ocols[2], mec=ocols[2], fmt='.', \
        ecolor='.7', capsize=0, label="$\mathrm{Subgiants}$", markersize=8)
pl.errorbar(vtbv, vtperiod, yerr=vtperiod_err, xerr=vtbv_err, color=ocols[2],\
        mec=ocols[2], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[:-5], cp[:-5], yerr=cp_err[:-5], xerr=cbv_err[:-5], color=ocols[3], \
        mec=ocols[3], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.errorbar(cbv[-5:], cp[-5:], yerr=cp_err[-5:], xerr=cbv_err[-5:], color=ocols[4], \
        mec=ocols[4], ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.legend(loc="upper left")
pl.ylim(0,70)
pl.savefig("cool_stars31")
