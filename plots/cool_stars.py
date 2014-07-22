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
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/no_dup.txt').T

# data = np.genfromtxt('/Users/angusr/Python/Gyro/data/old_data.txt').T
# data = np.genfromtxt('/Users/angusr/Python/Gyro/data/new_matched.txt').T
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/p_errs.txt').T
# data = np.genfromtxt('/Users/angusr/Python/Gyro/data/no_dup.txt').T
# data = np.genfromtxt('/Users/angusr/Python/Gyro/data/all_astero.txt').T
KID = data[0]
period = data[1]
period_err = data[2]
teff = data[3]
teff_err = data[4]
age = data[13]
age_err = data[14]
logg = data[10]
logg_errp = data[11]
logg_errm = data[12]
logg_err = .5*(logg_errp+logg_errm)
# bv = teff2bv_orig(teff, logg, np.ones_like(teff)*-.2)
# bv_err = np.ones_like(bv)*0.01
feh = np.ones_like(teff)*-.2
feh_err = np.ones_like(teff)*.3
bv, bv_err = teff2bv_err(teff, logg, feh, teff_err, logg_err, feh_err)

# # load travis and victor data
# data = np.genfromtxt("/Users/angusr/Python/Gyro/data/vandt.txt", skip_header=1).T
# vtKID = data[0]
# vtt = data[1]
# vtt_err = data[2]
# vta = data[3]
# vta_errp = data[4]
# vta_errm = data[5]
# vta_err = .5*(vta_errp+vta_errm)
# vtp = data[6]
# vtp_err = data[7]
# vtlogg = data[8]
# vtlogg_errp = data[9]
# vtlogg_errm = data[10]
# vtlogg_err = .5*(vtlogg_errp+vtlogg_errm)
# vtfeh = data[11]
# vtfeh_err = data[12]
# vtbv, vtbv_err = teff2bv_err(vtt, vtlogg, vtfeh, vtt_err, vtlogg_err, vtfeh_err)

# # remove travis and victor stars
# for i, kid in enumerate(vtKID):
#     l = KID==kid
#     if KID[l]:
#         age[l] = 0; age_err[l] = 0
#         bv[l] = 0; bv_err[l] = 0
#         period[l] = 0; period_err[l] = 0
#         teff[l] = 0; teff_err[l] = 0
#         logg[l] = 0; logg_err[l] = 0
#         logg_errp[l] = 0; logg_errm[l] = 0
#         feh[l] = 0; feh_err[l] = 0
# l = age>0
# age = age[l]; age_err = age_err[l]
# bv = bv[l]; bv_err = bv_err[l]
# period = period[l]; period_err = period_err[l]
# teff = teff[l]; teff_err = teff_err[l]
# logg = logg[l]; logg_err = logg_err[l]
# logg_errp = logg_errp[l]; logg_errm = logg_errm[l]
# feh = feh[l]; feh_err = feh_err[l]

# # add travis and victor stars
# print len(period), len(teff)
# age = np.concatenate((age, vta))
# age_err = np.concatenate((age_err, vta_err))
# period = np.concatenate((period, vtp))
# period_err = np.concatenate((period_err, vtp_err))
# teff = np.concatenate((teff, vtt))
# teff_err = np.concatenate((teff_err, vtt_err))
# bv = np.concatenate((bv, vtbv))
# bv_err = np.concatenate((bv_err, vtbv_err))
# logg = np.concatenate((logg, vtlogg))
# logg_err = np.concatenate((logg_err, vtlogg_err))
# logg_errp = np.concatenate((logg_errp, vtlogg_errp))
# logg_errm = np.concatenate((logg_errm, vtlogg_errm))
# feh = np.concatenate((feh, vtfeh))
# feh_err = np.concatenate((feh_err, vtfeh_err))

# cool dwarfs
pl.clf()
if colour==True:
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
else: pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], mec=ocols[1], \
        ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True:
    pl.xlim(.2,1.)
    pl.legend(loc="upper left")
    pl.ylim(0,70)
pl.savefig("cool_stars1")

# add amy's data
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/extra_amy_matched.txt').T
print len(period), len(data[1])
period = np.concatenate((period, data[1]))
period_err = np.concatenate((period_err, data[2]))
print len(period), len(teff), len(data[3])
teff = np.concatenate((teff, data[3]))
print len(teff)
teff_err = np.concatenate((teff_err, data[4]))
age = np.concatenate((age, data[13]))
age_err = np.concatenate((age_err, data[14]))
logg = np.concatenate((logg, data[10]))
logg_errp = np.concatenate((logg_errp, data[11]))
logg_errm = np.concatenate((logg_errm, data[12]))
logg_err = np.concatenate((logg_err, .5*(data[11]+data[12])))
feh = np.concatenate((feh, np.ones_like(data[2])*-.2))
feh_err = np.concatenate((feh_err, np.ones_like(data[2])*.3))
bva, bv_erra = teff2bv_err(data[3], data[10], np.ones_like(data[2])*.3, data[4], .5*(data[11]+data[12]), np.ones_like(data[2])*.3)
bv = np.concatenate((bv, bva))
bv_err = np.concatenate((bv_err, bv_erra))

# cool dwarfs with Amy
pl.clf()
if colour==True:
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
    pl.legend(loc="upper left")
    pl.ylim(0,70)
else: pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], mec=ocols[1], \
        ecolor='.7', capsize=0, fmt='.', markersize=8)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.savefig("cool_stars11")

# p vs bv paper
# load all astero data
data = np.genfromtxt('../data/all_astero.txt', skip_header=1).T
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
print len(teff)

# # cool dwarfs with clusters
# data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/hyades.txt',skip_header=2).T
# data3 = np.genfromtxt('/Users/angusr/Python/Gyro/data/praesepe.txt').T
# period = np.concatenate((period, data1[2], 1./data3[3]))
# period_err = np.concatenate((period_err, data1[3], (1./data3[3])*(data3[4]/data3[3])))
# bv2 = data1[0]
# bv3 = data3[5]-data3[6]
# bv = np.concatenate((bv, bv2, bv3))
# bv_err = np.concatenate((bv_err, data1[1], np.ones_like(data3[5])*.01))
# age = np.concatenate((age, data1[4], np.ones_like(data3[3])*.59))
# age_err = np.concatenate((age_err, data1[5], np.ones_like(data3[3])*.05))
# teff = np.concatenate((teff, bv2teff(bv2, np.ones_like(bv2)*4.5, np.ones_like(bv2)*-.2), bv2teff(bv3, np.ones_like(bv3)*4.5, np.ones_like(bv3)*-.2)))
# teff_err = np.concatenate((teff_err, np.ones_like(bv2)*100, np.ones_like(bv3)*100))

# cool dwarfs with clusters
data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/clusters.txt')
data1 = data1[:-1]
data1 = data1.T
period = np.concatenate((period, data1[2]))
period_err = np.concatenate((period_err, data1[3]))
bv = np.concatenate((bv, data1[0]))
bv_err = np.concatenate((bv_err, data1[1]))
age = np.concatenate((age, data1[4]))
age_err = np.concatenate((age_err, data1[5]))
teff = np.concatenate((teff, bv2teff(data1[0], np.ones_like(data1[0])*4.5, \
        np.ones_like(data1[0])*-.2)))
teff_err = np.concatenate((teff_err, np.ones_like(data1[0])*100))

data2 = np.genfromtxt('/Users/angusr/Python/Gyro/data/sun.txt',skip_header=2).T
bv_sun = data2[0]
bv_sun_err = data2[1]
period_sun = data2[2]
period_sun_err = data2[3]
age_sun = data[4]
age_sun_err = data[5]
teff_sun = 5778.
teff_sun_err = 10.

pl.clf()
if colour==True:
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
    pl.legend(loc="upper left")
    pl.ylim(0,70)
else:
    pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], mec=ocols[1], \
        ecolor='.7', capsize=0, fmt='.', markersize=8)
    xp, yp = iso(1000.)
    pl.plot(xp, yp, linestyle='--', color='.75', label=lab)
    xp, yp = iso(2000.)
    pl.plot(xp, yp, linestyle='--', color='.75')
    xp, yp = iso(5000.)
    pl.plot(xp, yp, linestyle='--', color='.75')
    xp, yp = iso(10000.)
    pl.plot(xp, yp, linestyle='--', color='.75')
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.savefig("cool_stars12")

pl.clf()
if colour==True:
    pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], mec=ocols[1], \
            ecolor='.7', capsize=0, fmt='.', markersize=8)
    pl.errorbar(bv_sun, period_sun, yerr=period_sun_err, xerr=bv_sun_err, \
            color='k', mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)
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
else:
    pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], mec=ocols[1], \
        ecolor='.7', capsize=0, fmt='.', markersize=8)
    pl.errorbar(teff_sun, period_sun, yerr=period_sun_err, xerr=teff_sun_err, color='k', \
            mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)


pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.savefig("cool_stars13")

# hot dwarfs
# l = teff > 6250
l = bv < kbreak
p, t, p_err, t_err, b, b_err = period[l], teff[l], period_err[l], teff_err[l], \
        bv[l], bv_err[l]
pl.clf()
if colour==False:
    pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], fmt='.', \
            mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~stars}$", \
            markersize=8)
    pl.errorbar(t, p, yerr=p_err, xerr=t_err, color=ocols[0], mec=ocols[0], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Hot~stars}$", markersize=8)
    pl.errorbar(teff_sun, period_sun, yerr=period_sun_err, xerr=teff_sun_err, color='k', \
            mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)
else:
    pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], fmt='.', \
            mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~stars}$", \
            markersize=8)
    pl.errorbar(b, p, yerr=p_err, xerr=b_err, color=ocols[0], mec=ocols[0], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Hot~stars}$", markersize=8)
    pl.errorbar(bv_sun, period_sun, yerr=period_sun_err, xerr=bv_sun_err, \
            color='k', mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)
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

pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.legend(loc="upper left")
pl.savefig("cool_stars2")

# subs
l2 = logg < 3.9
p2, t2, p_err2, t_err2, b2, b_err2 = period[l2], teff[l2], period_err[l2], teff_err[l2], \
        bv[l2], bv_err[l2]
pl.clf()
if colour==False:
    pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], fmt='.', \
            mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~dwarfs}$", \
            markersize=8)
    pl.errorbar(t, p, yerr=p_err, xerr=t_err, color=ocols[0], mec=ocols[0], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Hot~dwarfs}$", markersize=8)
    pl.errorbar(t2, p2, yerr=p_err2, xerr=t_err2, color=ocols[2], mec=ocols[2], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Subgiants}$", markersize=8)
    pl.errorbar(teff_sun, period_sun, yerr=period_sun_err, xerr=teff_sun_err, color='k', \
            mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)
else:
    pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], fmt='.', \
            mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~dwarfs}$", \
            markersize=8)
    pl.errorbar(b, p, yerr=p_err, xerr=b_err, color=ocols[0], mec=ocols[0], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Hot~dwarfs}$", markersize=8)
    pl.errorbar(b2, p2, yerr=p_err2, xerr=b_err2, color=ocols[2], mec=ocols[2], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Subgiants}$", markersize=8)
    pl.errorbar(bv_sun, period_sun, yerr=period_sun_err, xerr=bv_sun_err, \
            color='k', mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)

pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.legend(loc="upper left")
pl.savefig("cool_stars3")

# with isochrones
l2 = logg < 3.9
p2, t2, p_err2, t_err2, b2, b_err2 = period[l2], teff[l2], period_err[l2], teff_err[l2], \
        bv[l2], bv_err[l2]
pl.clf()
if colour==False:
    pl.errorbar(teff, period, yerr=period_err, xerr=teff_err, color=ocols[1], fmt='.', \
            mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~dwarfs}$", \
            markersize=8)
    pl.errorbar(t, p, yerr=p_err, xerr=t_err, color=ocols[0], mec=ocols[0], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Hot~dwarfs}$", markersize=8)
    pl.errorbar(t2, p2, yerr=p_err2, xerr=t_err2, color=ocols[2], mec=ocols[2], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Subgiants}$", markersize=8)
    pl.errorbar(teff_sun, period_sun, yerr=period_sun_err, xerr=teff_sun_err, color='k', \
            mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)
else:
    pl.errorbar(bv, period, yerr=period_err, xerr=bv_err, color=ocols[1], fmt='.', \
            mec=ocols[1], ecolor='.7', capsize=0, label="$\mathrm{Cool~dwarfs}$", \
            markersize=8)
    pl.errorbar(b, p, yerr=p_err, xerr=b_err, color=ocols[0], mec=ocols[0], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Hot~dwarfs}$", markersize=8)
    pl.errorbar(b2, p2, yerr=p_err2, xerr=b_err2, color=ocols[2], mec=ocols[2], fmt='.', \
            ecolor='.7', capsize=0, label="$\mathrm{Subgiants}$", markersize=8)
    pl.errorbar(bv_sun, period_sun, yerr=period_sun_err, xerr=bv_sun_err, \
            color='k', mec='k', ecolor='.7', capsize=0, fmt='.', markersize=10)
    xp, yp = iso(1000.)
    pl.plot(xp, yp, linestyle='--', color='.75', label=lab)
    xp, yp = iso(2000.)
    pl.plot(xp, yp, linestyle='--', color='.75')
    xp, yp = iso(5000.)
    pl.plot(xp, yp, linestyle='--', color='.75')
    xp, yp = iso(10000.)
    pl.plot(xp, yp, linestyle='--', color='.75')

pl.xlabel("$T_{eff}~\mathrm{(K)}$")
if colour==True: pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlim(7000, 4500)
if colour==True: pl.xlim(.2, 1.)
pl.legend(loc="upper left")
pl.ylim(0,70)
pl.savefig("cool_stars31")
