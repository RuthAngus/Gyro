# This calculates the typical uncertainity you would expect on the gyro age of a star
import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv_err

c = .45
nsamp = 1000

def period_model(par, age, bv, c):
    return par[0] * (age*1e3)**par[1] * (bv-c)**par[2]

def model(a_samp, n_samp, b_samp, age, bv, c):
    return a_samp * (age*1e3)**n_samp * (bv-c)**b_samp

def age_model(par, period, bv, c):
    return (period/par[0]/(bv-c)**par[2])**(1./par[1])/1000.

def amodel(a_samp, n_samp, b_samp, period, bv, c):
    return (period/a_samp/(bv-c)**b_samp)**(1./n_samp)/1000.

# calculate uncertainites from just the formal errors in the gyro relation
def MC_errors_simple(par, errp, errm, ap, bv, c, nsamp, xmodel, model):
    if bv<c: return 0., 0.

    a, n, b = par[:3]
    a_errp, n_errp, b_errp = errp[:3]
    a_errm, n_errm, b_errm = errp[:3]
    a_err, n_err, b_err = .5*(a_errp+a_errm), .5*(n_errp+n_errm), .5*(b_errp+b_errm)

    a_samp = a+a_err*np.random.randn(nsamp)
    n_samp = n+n_err*np.random.randn(nsamp)
    b_samp = b+b_err*np.random.randn(nsamp)

    mean = xmodel(par, ap, bv, c)
    distribution = model(a_samp, n_samp, b_samp, ap, bv, c)
    err = np.std(distribution)

    return mean, err

# calculate errors for each observation
def MC_errors(par, errp, errm, ap, ap_err, bv, bv_err, c, nsamp, xmodel, model):
    if bv<c: return 0., 0.

    a, n, b = par[:3]
    a_errp, n_errp, b_errp = errp[:3]
    a_errm, n_errm, b_errm = errp[:3]
    a_err, n_err, b_err = .5*(a_errp+a_errm), .5*(n_errp+n_errm), .5*(b_errp+b_errm)

    a_samp = a+a_err*np.random.randn(nsamp)
    n_samp = n+n_err*np.random.randn(nsamp)
    b_samp = b+b_err*np.random.randn(nsamp)
    ap_samp = ap+ap_err*np.random.randn(nsamp)
    bv_samp = bv+bv_err*np.random.randn(nsamp)
    l = (bv_samp>c) * (ap_samp>0)

    mean = xmodel(par, ap, bv, c)
    distribution = model(a_samp[l], n_samp[l], b_samp[l], ap_samp[l], bv_samp[l], c)
    err = np.std(distribution)

    return mean, err

# result = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters_45.txt').T
result = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters_45_2acf.txt').T
pars = result[0]
errp = result[1]
errm = result[2]

# load data to find typical observational uncertainties
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/garcia_all_astero.txt')
t = data[1]
t_err = data[2]
a = data[3]
a_errp = data[4]
a_errm = data[5]
a_err = .5*(a_errp+a_errm)
p = data[6]
p_err = data[7]
logg = data[8]
logg_errp = data[9]
logg_errm = data[10]
logg_err = .5*(logg_errp+logg_errm)
feh = data[11]
feh_err = data[12]
bv, bv_err = teff2bv_err(t, logg, feh, t_err, logg_err, feh_err)

periods1 = np.zeros_like(p)
period_errs1 = np.zeros_like(p)
periods2 = np.zeros_like(p)
period_errs2 = np.zeros_like(p)
for i, ps in enumerate(p):
    periods1[i], period_errs1[i] = MC_errors(pars, errp, errm, a[i], a_err[i], bv[i], bv_err[i], c, nsamp, period_model, model)
    periods2[i], period_errs2[i] = MC_errors_simple(pars, errp, errm, a[i], bv[i], c, nsamp, period_model, model)

l1 = (period_errs1>0) * (periods1>0)
l2 = (period_errs2>0) * (periods2>0)
print 'Typical uncertainty on rotation period, given typical observational uncertainties:'
print np.median(period_errs1[l1]/periods1[l1])*100, '%'

print 'Typical intrinsic uncertainty on rotation period:'
print np.median(period_errs2[l2]/periods2[l2])*100, '%'

ages1 = np.zeros_like(p)
age_errs1 = np.zeros_like(p)
ages2 = np.zeros_like(p)
age_errs2 = np.zeros_like(p)
for i, ass in enumerate(a):
    ages1[i], age_errs1[i] = MC_errors(pars, errp, errm, p[i], p_err[i], bv[i], bv_err[i], c, nsamp, age_model, model)
    ages2[i], age_errs2[i] = MC_errors_simple(pars, errp, errm, p[i], bv[i], c, nsamp, age_model, model)

l1 = (age_errs1>0) * (ages1>0)
l2 = (age_errs2>0) * (ages2>0)
print 'Typical uncertainty on age, given typical observational uncertainties:'
print np.median(age_errs1[l1]/ages1[l1])*100, '%'

print 'Typical intrinsic uncertainty on age:'
print np.median(age_errs2[l2]/ages2[l2])*100, '%'

pl.clf()
pl.hist(age_errs2[l2]/ages2[l2])
pl.show()
