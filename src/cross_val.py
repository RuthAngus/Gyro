import numpy as np
import matplotlib.pyplot as pl
from sklearn import cross_validation
from sklearn import datasets
from sklearn import svm
from all_plotting import load_dat
from gyro_like import period_model

def period_model(par, age, bv, c):
    return par[0] * (age*1e3)**par[1] * (bv-.45)**par[2]

def scoring(fname, test):
    # load test data
    a_test, age_err, age_errp, a_errm, p_test, p_err, bv_test, bv_err, g, g_err, \
            g_errp, g_errm, flag = load_dat(fname, test, cv=True)

    # load parameters from training data
    data = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters%s.txt'%fname).T
    pars = data[0][:3]

    # calculate score
    p_pred = period_model(pars, a_test, bv_test, .45)
    l = np.isfinite(p_pred)
    n = len(p_pred[l])

    # root mean squared error
    return np.sqrt((np.sum((p_test[l] - p_pred[l])**2))/n)

if __name__ == "__main__":
    trains = ['p_PF45', 'CF45', 'NF45', 'HF45']
    tests = ['CANH', 'APNH', 'CAPH', 'CAPN']
    RMS = []

    for i in range(len(tests)):
        # load test data
        a_test, age_err, age_errp, a_errm, p_test, p_err, bv_test, bv_err, g, g_err, \
                g_errp, g_errm, flag = load_dat(tests[i])

        # load parameters from training data
        data = np.genfromtxt('/Users/angusr/Python/noisy-plane/parameters%s.txt'%trains[i]).T
        pars = data[0][:3]

        # calculate score
        p_pred = period_model(pars, a_test, bv_test, .45)

        l = np.isfinite(p_pred)
        n = len(p_pred[l])

        # root mean squared error
        RMS.append(np.sqrt((np.sum((p_test[l] - p_pred[l])**2))/n))

        print trains[i], RMS[i]
