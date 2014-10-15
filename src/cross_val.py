import numpy as np
import matplotlib.pyplot as pl
from sklearn import cross_validation
from sklearn import datasets
from sklearn import svm
from all_plotting import load_dat
from gyro_like import period_model
from p_vs_bv import log_age_model

def period_model(par, age, bv, c):
    return par[0] * (age*1e3)**par[1] * (bv-.45)**par[2]

def scoring(fname, n, test):
    # load test data
    a_test, age_err, age_errp, a_errm, p_test, p_err, bv_test, bv_err, g, g_err, \
            g_errp, g_errm, flag = load_dat(fname, test, cv=True)

    # load parameters from training data
    data = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters%s%s.txt'%(n,fname)).T
    pars = data[0][:3]

    # calculate score
    p_pred = period_model(pars, a_test, bv_test, .45)
    l = np.isfinite(p_pred)
    n = len(p_pred[l])

    # root mean squared error
    return np.sqrt((np.sum((p_test[l] - p_pred[l])**2))/n)

def loo():

#     load field star data
    data = np.genfromtxt('/Users/angusr/Python/Gyro/data/clusters.txt').T
    print np.shape(data)
    bv, bv_err, p, p_err, a, a_err, g, g_err, flag = data
    l = flag==8
    bv, bv_err, p, p_err, a, a_err, g, g_err, flag = bv[l], \
            bv_err[l], p[l], p_err[l], a[l], a_err[l], g[l], \
            g_err[l], flag[l]

    for i, age in enumerate(a):
        print i, age
        raw_input('enter')
        data = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters_loo_%sACHF45.txt'                              % (i+1)).T
        pars = np.zeros(4)
        pars[:3] = data[0][:3]
        pars[3] = 0.45
        print age, log_age_model(pars, np.log10(p[i]), bv[i])

if __name__ == "__main__":
    loo()
    raw_input('enter')

    trains = ['CF45', 'HF45', 'PF45', 'PF5']
    tests = ['AHP', 'ACP', 'ACH', 'ACH']
    RMS = []

    for i in range(len(tests)):
        # load test data
        a_test, age_err, age_errp, a_errm, p_test, p_err, bv_test, bv_err, g, g_err, \
                g_errp, g_errm, flag = load_dat(tests[i], False, False)
        print len(p_test)

        # load parameters from training data
        data = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters%s.txt'%trains[i]).T
        pars = data[0][:3]

        # calculate score
        p_pred = period_model(pars, a_test, bv_test, .45)

        l = np.isfinite(p_pred)
        n = len(p_pred[l])

        # root mean squared error
        RMS.append(np.sqrt((np.sum(((p_test[l] - p_pred[l])**2)))/n))

#         # reduced chi squared
#         print 'chi', sum(((p_test[l] - p_pred[l])**2)/p_err[l]**2)/n
#
#         # log likelihood
        L = (.5*np.sum(((p_pred[l] - p_test[l]) / p_err[l])**2)) #- np.log(2*np.pi*n)
#         print p_err[l]
#         raw_input('enter')
#         print 'likelihood = ', L
#
        print trains[i], RMS[i], L
