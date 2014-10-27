import numpy as np
import h5py
import matplotlib.pyplot as plt
import random

def age_model(par, p, bv):
    c = 0.45
    return 10**((np.log10(p) - np.log10(par[0]) - \
            par[2]*np.log10(bv - c)) / par[1]) / 1000.

def ndim_age_model(par, p, bv, c):
    return 10**((np.log10(p) - np.log10(par[:, 0]) - \
            par[:, 2]*np.log10(bv - c)) / par[:, 1]) / 1000.

def posterior(mu, sig, n):
    return sig * np.random.randn(n) + mu

def age_calc(age, age_err, bv, bv_err, p, p_err, flatchain, c):
    nsamp = len(flatchain)
    a_samp = posterior(age, age_err, nsamp)
    bv_samp = posterior(bv, bv_err, nsamp)
    p_samp = posterior(p, p_err, nsamp)
    age_samples = ndim_age_model(flatchain[:, :3], p_samp, bv_samp, c)
    p = np.percentile(age_samples, [16, 50, 84])
    return np.array([p[1], p[2]-p[1], p[1]-p[0]]), age_samples

if __name__ == "__main__":

#     fname = 'ACHF45irfm2'
    fnames = ['loo_5ACHF45', 'loo_4ACHF45', 'loo_3ACHF45', 'loo_2ACHF45', 'loo_1ACHF45']

    # alpha cen, 18 Sco, 16 Cyg B, sun
    age = [6.0, 6.0, 3.66, 6.4, 4.568]
    age_err = [1.0, 1.0, 0.2, 0.4, 0.001]
    bv = [0.69, 0.90, 0.64, 0.66, 0.65]
    bv_err = [0.01, 0.01, 0.01, 0.01, 0.01]
    p = [28.8, 38.7, 22.7, 31.5, 26.9]
    p_err = [2.5, 5.0, 0.5, 6.5, 0.1]

    for i, fname in enumerate(fnames):
        with h5py.File("samples_%s" %fname, "r") as f:
            samples = f["samples"][:, 50:, :]
        nwalkers, n, ndim = samples.shape
        flatchain = samples.reshape((-1, ndim))
        mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                          zip(*np.percentile(flatchain, [16, 50, 84], axis=0)))
        mres = np.array(mcmc_result)[:, 0]

        # Barnes
        bpar = [.7725, .5189, .601, .4]
        bpar_err = [.011, .007, .024, 0.]
        nsamp, ndim = np.shape(flatchain)
        bflatchain = np.empty_like(flatchain)
        for j in range(3):
            random.seed(123)
            bflatchain[:, j] = bpar_err[j]*np.random.randn(nsamp)+bpar[j]

        # M&H
        mpar = [.407, .566, .325, .495]
        mpar_err = [.021, .008, .024, .01]
    #     bpar_err = [.0070, .011, .024, 0.]
        mflatchain = np.empty_like(flatchain)
        for j in range(4):
            random.seed(123)
            mflatchain[:, j] = mpar_err[j]*np.random.randn(nsamp)+mpar[j]

        # mine without using posterior samples
        par2 = mres[:3]
        par_err2 = np.mean(np.array(mcmc_result)[:3, 1:], axis=1)
        flatchain2 = np.empty_like(flatchain)
        for j in range(3):
            random.seed(123)
            flatchain2[:, j] = par_err2[j]*np.random.randn(nsamp)+par2[j]

#         for i in range(0, len(age)):
        print i
        print age[i]
        gyro_age, samples = age_calc(age[i], age_err[i], bv[i], bv_err[i],
                       p[i], p_err[i], flatchain, 0.45)
        print 'my age:', gyro_age
        age2, samples2 = age_calc(age[i], age_err[i], bv[i], bv_err[i],
                       p[i], p_err[i], flatchain2, 0.45)
        print 'my age2:', age2
        b_age, bsamples = age_calc(age[i], age_err[i], bv[i], bv_err[i],
                       p[i], p_err[i], bflatchain, 0.4)
        print 'Barnes age:', b_age
        m_age, msamples = age_calc(age[i], age_err[i], bv[i], bv_err[i],
                p[i], p_err[i], mflatchain, mflatchain[:, 3])
        print 'M&H age:', m_age

        plt.clf()
        plt.hist(age_err[i]*np.random.randn(len(flatchain))+age[i],
                 50, normed=1, alpha=.5, color='w', edgecolor='b',
                 histtype='stepfilled') # true
        plt.hist(samples, 50, normed=1, alpha=.5, color='w', edgecolor='r',
                 histtype='stepfilled') # mine
        plt.hist(bsamples, 50, normed=1, alpha=.5, color='w', edgecolor='k',
                 histtype='stepfilled') # barnes
        plt.hist(msamples, 50, normed=1, alpha=.3, color='w', edgecolor='k',
                 histtype='stepfilled') # m&h
        plt.axvline(age[i], color='b')
        plt.axvline(gyro_age[0], color='r')
        plt.axvline(b_age[0], color='k')
        plt.axvline(m_age[0], color='k', alpha=.5)
        print '\n', 'me ', abs(age[i]-gyro_age[0]), \
                'barnes ', abs(age[i]-b_age[0]), \
                'm&h ', abs(age[i]-m_age[0]), '\n'
        plt.ylim(0, 0.7)
        plt.xlim(0, 10)

    # calculate minimum uncertainty
    gyro_age, samples = age_calc(age[-1], 1e-16, bv[-1], 1e-16,
                                 p[-1], 1e-16, flatchain, 0.45)
    print 'minimum uncertainty', gyro_age
