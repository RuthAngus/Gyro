import numpy as np
import matplotlib.pyplot as pl

def simple_sample(age_obs, age_err, bv_obs, bv_err, period_obs, period_err,
                  logg_obs, logg_err, nsamp, s):
    np.random.seed(s)
    age_samp = np.vstack([x0+xe*np.random.randn(nsamp)
                         for x0, xe in zip(age_obs, age_err)])
    np.random.seed(s)
    bv_samp = np.vstack([x0+xe*np.random.randn(nsamp)
                        for x0, xe in zip(bv_obs, bv_err)])
    np.random.seed(s)
    logg_samp = np.vstack([x0+xe*np.random.randn(nsamp)
                          for x0, xe in zip(logg_obs, logg_err)])
    np.random.seed(s)
    period_samp = np.vstack([x0+xe*np.random.randn(nsamp)
                            for x0, xe in zip(period_obs, period_err)])
    return age_samp, bv_samp, logg_samp, period_samp

def asym_sample(age_obs, age_errp, age_errm, bv_obs, bv_err, period_obs, period_err,
                logg_obs, logg_errp, logg_errm, nsamp, s):

    nrough = 1000
    age_samp = np.zeros((len(age_obs), nsamp))
    logg_samp = np.zeros((len(age_obs), nsamp))
    for i, obs in enumerate(age_obs):

        np.random.seed(s)
        age_sampp = obs+age_errp[i]*np.random.randn(nrough)
        np.random.seed(s)
        age_sampm = obs+age_errm[i]*np.random.randn(nrough)
        varp, varm = np.std(age_sampp), np.std(age_sampm)
        r = varm/varp
        age_samp_tempp = obs+age_errp[i]*np.random.randn(nrough)
        age_samp_tempm = obs+age_errm[i]*np.random.randn(nrough*r)
        lp, lm = age_samp_tempp > obs, age_samp_tempm < obs
        age_samp_temp = np.concatenate((age_samp_tempp[lp], age_samp_tempm[lm]))
        ra = np.random.randint(0, len(age_samp_temp), nsamp)
        age_samp[i, :] = age_samp_temp[ra]

        np.random.seed(s)
        logg_sampp = logg_obs[i]+logg_errp[i]*np.random.randn(nrough)
        np.random.seed(s)
        logg_sampm = logg_obs[i]+logg_errm[i]*np.random.randn(nrough)
        varp, varm = np.std(logg_sampp), np.std(logg_sampm)
        r = varm/varp
        logg_samp_tempp = obs+logg_errp[i]*np.random.randn(nrough)
        logg_samp_tempm = obs+logg_errm[i]*np.random.randn(nrough*r)
        lp, lm = logg_samp_tempp > obs, logg_samp_tempm < obs
        logg_samp_temp = np.concatenate((logg_samp_tempp[lp], logg_samp_tempm[lm]))
        ra = np.random.randint(0, len(logg_samp_temp), nsamp)
        logg_samp[i, :] = logg_samp_temp[ra]

#         pl.clf()
#         pl.subplot(2,1,1)
#         pl.hist(age_sampp, 50, normed=True, color='r', alpha=.3)
#         pl.hist(age_sampm, 50, normed=True, color='b', alpha=.3)
#         pl.subplot(2,1,2)
#         pl.hist(age_samp[i, :], 100, normed=True, alpha=.3)
#         pl.xlim(-2, 14)
#         pl.show()
#
#         pl.clf()
#         pl.subplot(2,1,1)
#         pl.hist(logg_sampp, 50, normed=True, color='r', alpha=.3)
#         pl.hist(logg_sampm, 50, normed=True, color='b', alpha=.3)
#         pl.subplot(2,1,2)
#         pl.hist(logg_samp[i, :], 100, normed=True, alpha=.3)
#         pl.show()
#         raw_input('enter')

    np.random.seed(s)
    bv_samp = np.vstack([x0+xe*np.random.randn(nsamp) for x0, xe in zip(bv_obs, bv_err)])
    np.random.seed(s)
    period_samp = np.vstack([x0+xe*np.random.randn(nsamp) for x0, xe in zip(period_obs, period_err)])
    return age_samp, bv_samp, period_samp, logg_samp
