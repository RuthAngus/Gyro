import matplotlib.pyplot as pl
import numpy as np
import emcee
import triangle
from all_plotting import load_dat
import h5py
import datetime
from gyro_like import lnlike, period_model
from sampling import asym_sample

def lnprior(m):
    if -10. < m[0] < 10. and .0< m[1] < 1. and 0. < m[2] < 1. \
            and 0 < m[3] < 30. and 0 < m[4] < 100.\
            and 0 < m[5] < 30. and 0 < m[6] < 100.\
            and 0 < m[7] < 30. and 0 < m[8] < 100.\
            and 0. < m[9] < 1.:
        return 0.0
    return -np.inf

def lnprob(m, age_samp, bv_samp, period_samp, logg_samp, age_obs, age_err, \
        bv_obs, bv_err, period_obs, period_err, logg_obs, \
        logg_err, c):
    lp = lnprior(m)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(m, age_samp, bv_samp, \
            period_samp, logg_samp, age_obs, age_err, bv_obs, bv_err, period_obs, \
            period_err, logg_obs, logg_err, c)

def MCMC(fname, c):

    # load MAP values
    try:
        params = np.genfromtxt('%s/results/parameters%s.txt' % (DIR, fname)).T
        par_true = params[0]
        print 'initialising with MAP values'
        print par_true
    except:
#         par_true = [0.7725, 0.5189, .601, 5., 10., \ # Barnes
        par_true = [0.6, 0.5189, .601, 5., 10., \
                8., 30., 9., 5., .67] # better initialisation
        print par_true

    # load real data
    age_obs, age_err, age_errp, age_errm, period_obs, period_err, bv_obs, bv_err, \
            logg_obs, logg_err, logg_errp, logg_errm, flag = load_dat(fname, False, False)

    pl.clf()
    pl.errorbar(age_obs, period_obs, xerr=age_err, yerr=period_err, fmt='k.',
                 capsize=0, ecolor='.8')
    pl.show()

    # Now generate samples
    age_samp, bv_samp, logg_samp, period_samp = \
            asym_sample(age_obs, age_errp, age_errm, bv_obs, bv_err, period_obs,
                        period_err, logg_obs, logg_errp, logg_errm, 100, 12)

    print 'initial likelihood = ', lnlike(par_true, age_samp, bv_samp, \
            period_samp, logg_samp, age_obs, age_err, bv_obs, bv_err, \
            period_obs, period_err, logg_obs, logg_err, c)

    nwalkers, ndim = 32, len(par_true)
    p0 = [par_true+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    args = (age_samp, bv_samp, period_samp, logg_samp, age_obs, age_err, bv_obs, \
            bv_err, period_obs, period_err, logg_obs, logg_err, c)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = args)

    print("Burn-in")
#     p0, lp, state = sampler.run_mcmc(p0, 5000)
    p0, lp, state = sampler.run_mcmc(p0, 5)
    sampler.reset()
    print("Production run")
#     nstep = 20000
#     nruns = 2000.
    nstep = 2
    nruns = 2

    for j in range(int(nstep/nruns)):

        print fname
        print datetime.datetime.now()
        print 'run', j
        p0, lp, state = sampler.run_mcmc(p0, nruns)

        flat = sampler.chain[:, 50:, :].reshape((-1, ndim))
        mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                          zip(*np.percentile(flat, [16, 50, 84], axis=0)))
        mres = np.array(mcmc_result)[:, 0]
        print 'mcmc_result = ', mres

        print "saving samples"
        f = h5py.File("%s/results/samples_%s" % (DIR, fname), "w")
        data = f.create_dataset("samples", np.shape(sampler.chain))
        data[:,:] = np.array(sampler.chain)
        f.close()

        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

        likelihood = lnlike(mres, age_samp, bv_samp, period_samp,
                            logg_samp, age_obs, age_err, bv_obs, bv_err,
                            period_obs, period_err, logg_obs, logg_err, c)
        print 'likelihood = ', likelihood

    # save parameters and print to screen
    print 'initial values', par_true
    np.savetxt("%s/results/parameters%s.txt" % (DIR, fname), np.array(mcmc_result))
    mcmc_result = np.array(mcmc_result)[:, 0]
    mcmc_result = [mcmc_result[0], mcmc_result[1], mcmc_result[2], c, \
            mcmc_result[3], mcmc_result[4]]

    # save samples
    f = h5py.File("%s/results/samples_%s" % (DIR, fname), "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:,:] = np.array(sampler.chain)
    f.close()

if __name__ == "__main__":

    DIR = '/Users/angusr/Python/Gyro'

    # proper runs
    fname = 'asym_ACHF45'

    print fname, "fname"
    MCMC(fname, .45)