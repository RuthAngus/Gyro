import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv_err, bv2teff, teff2bv_orig

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
# ocols = ['#FF9933','#66CCCC','#FF33CC','#3399FF','#CC0066','#99CC99','#9933FF','#CC0000']
ocols = ['#FF9933','#009999','#FF33CC','#3399FF','#CC0066','#99CC99','#9933FF','#CC0000']
ms = 4
kbreak = .45
subgiants = 3.9

def model(pars, log_age, bv):
    return np.log10(pars[0]) + pars[1]*log_age + pars[2]*np.log10(bv-pars[3])

data = np.genfromtxt('/Users/angusr/Python/Gyro/data/garcia_all_astero.txt')
teff = data[1]
teff_err = data[2]
age = data[3]
age_errp = data[4]
age_errm = data[5]
age_err = .5*(age_errp+age_errm)
period = data[6]
period_err = data[7]
logg = data[8]
logg_errp = data[9]
logg_errm = data[10]
logg_err = .5*(logg_errp+logg_errm)
feh = data[11]
feh_err = data[11]
bv, bv_err = teff2bv_err(teff, logg, feh, teff_err, logg_err, feh_err)
a_errp = age_errp
a_errm = age_errm
p_err = period_err

pl.clf()
pl.errorbar(age, period, xerr = (a_errp, a_errm), yerr = p_err, \
            color = ocols[1], fmt = 'o', mec = ocols[1], capsize = 0, \
    markersize = ms, ecolor = '0.8')
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 70)
pl.savefig('cool_stars5')

# add Amy
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/extra_amy_matched.txt').T
period = np.concatenate((period, data[1]))
period_err = np.concatenate((period_err, data[2]))
teff = np.concatenate((teff, data[3]))
teff_err = np.concatenate((teff_err, data[4]))
age = np.concatenate((age, data[13]))
age_err = np.concatenate((age_err, data[14]))
logg = np.concatenate((logg, data[10]))
bv = teff2bv_orig(teff, logg, np.ones_like(teff)*-.2)
bv_err = np.ones_like(bv)*.01

a_errp = age_err
a_errm = age_err
p_err = period_err

pl.clf()
pl.errorbar(age, period, xerr = (a_errp, a_errm), yerr = p_err, \
            color = ocols[1], fmt = 'o', mec = ocols[1], capsize = 0, \
    markersize = ms, ecolor = '0.8')
print len(age)
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 70)
pl.savefig('cool_stars51')

pl.clf()
pl.errorbar(age, period, xerr=(a_errp, a_errm), yerr=p_err, \
            color='k', fmt='o', mec='k', capsize=0, \
    markersize=2, ecolor='0.8')
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 70)
# pl.loglog()
pl.savefig('p_vs_a_paper')

# # add clusters
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

logg = np.concatenate((logg, np.ones_like(data1[0])*4.5))

a_errp = age_err
a_errm = age_err
p_err = period_err
pl.clf()
pl.errorbar(age, period, xerr = (a_errp, a_errm), yerr = p_err, \
            color = ocols[1], fmt = 'o', mec = ocols[1], capsize = 0, \
    markersize = ms, ecolor = '0.8')
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 70)
pl.savefig('cool_stars52')

# add sun
data2 = np.genfromtxt('/Users/angusr/Python/Gyro/data/sun.txt',skip_header=2)
bv_sun = data2[0]
bv_sun_err = data2[1]
period_sun = data2[2]
period_sun_err = data2[3]
age_sun = data2[4]
age_sun_err = data2[5]
teff_sun = 5778.
teff_sun_err = 10.

pl.clf()
pl.errorbar(age, period, xerr = (a_errp, a_errm), yerr = p_err, \
            color = ocols[1], fmt = 'o', mec = ocols[1], capsize = 0, \
    markersize = ms, ecolor = '0.8')
pl.errorbar(age_sun, period_sun, xerr = age_sun_err, yerr = period_sun_err, \
            color = 'k', fmt = 'o', mec = 'k', capsize = 0, \
            markersize = ms, ecolor = '0.8')
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 70)
pl.savefig('cool_stars6')

# hot stars
l1 = bv < kbreak
pl.clf()
pl.errorbar(age, period, xerr = (a_errp, a_errm), yerr = p_err, \
            color = ocols[1], fmt = 'o', mec = ocols[1], capsize = 0, \
    markersize = ms, ecolor = '0.8', label="$\mathrm{Cool~stars}$")
pl.errorbar(age_sun, period_sun, xerr = age_sun_err, yerr = period_sun_err, \
            color = 'k', fmt = 'o', mec = 'k', capsize = 0, \
            markersize = ms, ecolor = '0.8')
pl.errorbar(age[l1], period[l1], xerr = (a_errp[l1], a_errm[l1]), yerr = p_err[l1], \
            color = ocols[0], fmt = 'o', mec = ocols[0], capsize = 0, \
    markersize = ms, ecolor = '0.8', label='$\mathrm{Hot~stars}$', zorder = 3)
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 70)
pl.legend(loc="upper right")
pl.savefig('cool_stars7')

# subgiants
l2 = logg < subgiants
pl.clf()
pl.errorbar(age, period, xerr = (a_errp, a_errm), yerr = p_err, \
            color = ocols[1], fmt = 'o', mec = ocols[1], capsize = 0, \
    markersize = ms, ecolor = '0.8', label="$\mathrm{Cool~dwarfs}$")
pl.errorbar(age[l1], period[l1], xerr = (a_errp[l1], a_errm[l1]), yerr = p_err[l1], \
            color = ocols[0], fmt = 'o', mec = ocols[0], capsize = 0, \
    markersize = ms, ecolor = '0.8', label='$\mathrm{Hot~dwarfs}$', zorder = 3)
pl.errorbar(age[l2], period[l2], xerr = (a_errp[l2], a_errm[l2]), yerr = p_err[l2], \
            color = ocols[2], fmt = 'o', mec = ocols[2], capsize = 0, \
            markersize = ms, ecolor = '0.8', label='$\mathrm{Subgiants}$', zorder = 4)
pl.errorbar(age_sun, period_sun, xerr = age_sun_err, yerr = period_sun_err, \
            color = 'k', fmt = 'o', mec = 'k', capsize = 0, \
            markersize = ms, ecolor = '0.8')
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 120)
pl.legend(loc="upper right")
pl.savefig('cool_stars8')

# remove
l3 = (logg > subgiants) * (bv > .45)
pl.clf()
pl.errorbar(age[l3], period[l3], xerr = (a_errp[l3], a_errm[l3]), yerr = p_err[l3], \
            color = '.2', fmt = 'o', mec = '.2', capsize = 0, \
    markersize = ms, ecolor = '0.8', label="$\mathrm{Cool~dwarfs}$")
pl.ylabel("$P_{rot}~\mathrm{(days)}$")
pl.xlabel("$\mathrm{Age~(Gyr)}$")
pl.xlim(0.1, 15)
pl.ylim(0, 80)
pl.legend(loc="upper right")
pl.savefig('cool_stars9')
