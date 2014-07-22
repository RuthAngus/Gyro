import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv_err

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)

ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066',
'#99CC99', '#9933FF', '#CC0000', '#99CC00']

# load asteroseismic data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_all_astero.txt")
KID = data[0]
t1 = data[1]
t1_err = data[2]
a1 = data[3]
a1_errp = data[4]
a1_errm = data[5]
a1_err = .5*(a1_errp+a1_errm)
p1 = data[6]
p1_err = data[7]
logg1 = data[8]
logg1_errp = data[9]
logg1_errm = data[10]
logg1_err = .5*(logg1_errp+logg1_errm)
feh1 = data[11]
feh1_err = data[12]
bv1, bv1_err = teff2bv_err(t1, logg1, feh1, t1_err, logg1_err, feh1_err)

# p1 = data[1]
# p1_err = data[2]
# t1 = data[3]
# t1_err = data[4]
# feh1 = data[5]
# feh1_err = data[6]
# logg1 = data[10]
# logg1_errp = data[11]
# logg1_errm = data[12]
# logg1_err = .5*(logg1_errp+logg1_errm)
# a1 = data[13]
# a1_errp = data[14]
# a1_errm = data[15]
# a1_err = .5*(a1_errp+a1_errm)
# bv1, bv1_err = teff2bv_err(t1, logg1, feh1, t1_err, logg1_err, feh1_err)

# load travis and victor data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/vandt.txt", skip_header=1).T
vtKID = data[0]
vtt = data[1]
vtt_err = data[2]
vta = data[3]
vta_errp = data[4]
vta_errm = data[5]
vta_err = .5*(vta_errp+vta_errm)
vtp = data[6]
vtp_err = data[7]
vtlogg = data[8]
vtlogg_errp = data[9]
vtlogg_errm = data[10]
vtlogg_err = .5*(vtlogg_errp+vtlogg_errm)
vtfeh = data[11]
vtfeh_err = data[12]
vtbv, vtbv_err = teff2bv_err(vtt, vtlogg, vtfeh, vtt_err, vtlogg_err, vtfeh_err)

l = (t1>0)*(logg1>0)
hot = t1[l]>6250
sub = logg1[l]<4.
pl.clf()
pl.errorbar(t1[l], logg1[l], xerr=t1_err[l], yerr=(logg1_errp[l], logg1_errm[l]), \
        fmt='k.', capsize=0, ecolor='.7', mec='k')
pl.errorbar(t1[l][hot], logg1[l][hot], xerr=t1_err[l][hot], yerr=(logg1_errp[l][hot], \
        logg1_errm[l][hot]), fmt='r.', capsize=0, ecolor='.7', mec='r')
pl.errorbar(t1[l][sub], logg1[l][sub], xerr=t1_err[l][sub], yerr=(logg1_errp[l][sub], \
        logg1_errm[l][sub]), fmt='b.', capsize=0, ecolor='.7', mec='b')
hot = vtt>6250
sub = vtlogg<4.
pl.errorbar(vtt, vtlogg, xerr=vtt_err, yerr=(vtlogg_errp, vtlogg_errm), \
        capsize=0, ecolor='.7', mec='k', fmt="k^", markersize=4)
pl.errorbar(vtt[hot], vtlogg[hot], xerr=vtt_err[hot], yerr=(vtlogg_errp[hot], \
        vtlogg_errm[hot]), capsize=0, ecolor='.7', mec='r', fmt="r^", markersize=4)
pl.errorbar(vtt[sub], vtlogg[sub], xerr=vtt_err[sub], yerr=(vtlogg_errp[sub], \
        vtlogg_errm[sub]), capsize=0, ecolor='.7', mec='b', fmt="b^", markersize=4)
pl.xlabel("$T_{eff}~\mathrm{(K)}$")
pl.ylabel("$\mathrm{log}~g$")
pl.ylim(pl.gca().get_ylim()[::-1])
pl.xlim(pl.gca().get_xlim()[::-1])
pl.savefig("logg_vs_t_paper")

# load cluster data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt").T
bv2 = data[0]
bv2_err = data[1]
p2 = data[2]
p2_err = data[3]
a2 = data[4]
a2_err = data[5]

# load travis data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/travis_shortlist.txt",skip_header=1).T
tKID = data[0]
ta = data[5]
ta_err = data[6]
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/travis_actual_periods.txt").T
tp = data[1]
tp_err = data[2]
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/travis_teffs.txt",skip_header=1).T
tt = data[1]
tt_err = data[2]
tlogg = data[3]
tlogg_err = data[4]
tfeh = data[5]
tfeh_err = data[6]
tbv, tbv_err = teff2bv_err(tt, tlogg, tfeh, tt_err, tlogg_err, tfeh_err)

# load Victor data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/Victor_params.txt",skip_header=1).T
vKID = data[0]
vfeh = data[1]
vfeh_err = data[2]
vt = data[3]
vt_err = data[4]
vlogg = data[11]
vlogg_errp = data[12]
vlogg_errm = data[13]
vlogg_err = .5*(vlogg_errp+vlogg_errm)
va = data[14]
va_errp = data[15]
va_errm = data[16]
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/Victor_p_errs.txt").T
vp = data[1]
vp_err = data[2]
vbv, vbv_err = teff2bv_err(vt, vlogg, vfeh, vt_err, vlogg_err, vfeh_err)

# remove travis' stars
for i, kid in enumerate(tKID):
    l = KID==kid
    if KID[l]:
        a2[l] = 0; a2_err[l] = 0
        bv2[l] = 0; bv2_err[l] = 0
        p2[l] = 0; p2_err[l] = 0
l = a2>0
a2 = a2[l]; a2_err = a2_err[l]
bv2 = bv2[l]; bv2_err = bv2_err[l]
p2 = p2[l]; p2_err = p2_err[l]

# remove victors' stars
for i, kid in enumerate(vKID):
    l = KID==kid
    if KID[l]:
        a2[l] = 0; a2_err[l] = 0
        bv2[l] = 0; bv2_err[l] = 0
        p2[l] = 0; p2_err[l] = 0
l = a2>0
a2 = a2[l]; a2_err = a2_err[l]
bv2 = bv2[l]; bv2_err = bv2_err[l]
p2 = p2[l]; p2_err = p2_err[l]

# plot p_vs_bv for all stars
sun = a2 == 4.568
pl.clf()
# astero stars
pl.errorbar(bv1, p1, xerr=bv1_err, yerr=p1_err, fmt='k.', capsize=0, ecolor='.7')
# cluster + field
pl.errorbar(bv2, p2, xerr=bv2_err, yerr=p2_err, fmt='.', color='r', capsize=0, ecolor='.7')
# sun
pl.errorbar(bv2[sun], p2[sun], xerr=bv2_err[sun], yerr=p2_err[sun], fmt='*', color='k', capsize=0, ecolor='.7',\
        markersize = 8, mec='k')
# Travis and victor
pl.errorbar(tbv, tp, xerr=tbv_err, yerr=tp_err, fmt='.', color='b', capsize=0, ecolor='.7')
pl.errorbar(vbv, vp, xerr=vbv_err, yerr=vp_err, fmt='.', color='b', capsize=0, ecolor='.7')
pl.xlabel("$\mathrm{B-V}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.ylim(0,70)
pl.savefig("p_vs_bv_paper2")

# plot p_vs_a for all stars
# a1 = np.log10(a1)
# a2 = np.log10(a2)
# p1 = np.log10(p1)
# p2 = np.log10(p2)
pl.clf()
l = a1>0
a1 = a1[l]
p1 = p1[l]
a1_err = a1_err[l]
p1_err = p1_err[l]
pl.errorbar(a1, p1, xerr=a1_err, yerr=p1_err, fmt='k.', capsize=0, ecolor='.7')
pl.errorbar(a2, p2, xerr=a2_err, yerr=p2_err, fmt='.', color='r', capsize=0, ecolor='.7')
pl.errorbar(a2[sun], p2[sun], xerr=a2_err[sun], yerr=p2_err[sun], fmt='*', color='k', capsize=0, ecolor='.7',\
        markersize = 8, mec='k')
pl.errorbar(ta, tp, xerr=ta_err, yerr=tp_err, fmt='.', color='b', capsize=0, ecolor='.7')
# pl.plot(a1, p1, 'k.')
sun = 10**a2==4.568
# pl.plot(a2, p2, '.', color='.5')
# pl.plot(a2[sun], p2[sun], '.', color='r')
# print a2[sun]
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.ylim(0,70)
pl.xlim(0,15)
# pl.loglog()
pl.savefig("p_vs_a_paper2")

# data = np.empty((len(tKID)+len(vKID), 14))
# # ["KID", "t", "t_err", "a", "a_errp", "a_errm", "p", "p_err", "logg", "logg_errp", "logg_errm", "feh", "feh_err"]
# data[0:] = np.concatenate((tKID, vKID))
# data[1:] = np.concatenate((tt, vt))
# data[2:] = np.concatenate((tt_err, vt_err))
# data[3:] = np.concatenate((ta, va))
# data[4:] = np.concatenate((ta_err, va_errp))
# data[5:] = np.concatenate((ta_err, va_errm))
# data[6:] = np.concatenate((tp, vp))
# data[7:] = np.concatenate((tp_err, vp_err))
# data[8:] = np.concatenate((tlogg, vlogg))
# data[9:] = np.concatenate((tlogg_err, vlogg_errp))
# data[10:] = np.concatenate((tlogg_err, vlogg_errm))
# data[11:] = np.concatenate((tfeh, vfeh))
# data[12:] = np.concatenate((tfeh_err, vfeh_err))
#
# np.savetxt("/Users/angusr/Python/Gyro/data/vandt.txt", data.T)
