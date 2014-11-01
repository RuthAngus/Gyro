import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from teff_bv import teff2bv_orig, teff2bv_err

def load_dat(fname, tn, cv):

    data = np.genfromtxt('/Users/angusr/Python/Gyro/data/garcia_irfm.txt')
    KID = data[0]
    t = data[1]
    p = data[6]
    g = data[8]
    p_err = data[7]

    # remove periods <= 0 and teff < 100
    l = (p > 0.)*(t > 100.)*(g > 0.) * (p_err > 0.)

    KID = data[0][l]
    p = p[l]
    p_err = data[7][l]
    t = t[l]
    t_err = data[2][l]
    g = g[l]
    g_errp = data[9][l]
    g_errm = data[10][l]
    a = data[3][l]
    a_errp = data[4][l]
    a_errm = data[5][l]
    feh = data[11][l]
    feh_err = data[12][l]
    flag = data[13][l]

    # convert temps to bvs
    bv_obs, bv_err = teff2bv_err(t, g, feh, t_err, .5*(g_errp+g_errm), feh_err)

    # add clusters
    data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt", skip_header=1).T
    bv_obs = np.concatenate((bv_obs, data[0]))
    bv_err = np.concatenate((bv_err, data[1]))
    p = np.concatenate((p, data[2]))
    p_err = np.concatenate((p_err, data[3]))
    a = np.concatenate((a, data[4]))
    a_errp = np.concatenate((a_errp, data[5]))
    a_errm = np.concatenate((a_errm, data[5]))
    g = np.concatenate((g, data[6]))
    g_errp = np.concatenate((g_errp, data[7]))
    g_errm = np.concatenate((g_errm, data[7]))
    flag = np.concatenate((flag, data[8]))

    # obviously comment these lines out if you want to use temps
    t = bv_obs
    t_err = bv_err

    # select star group
    flag[flag==2] = 9
    flag -= 3
    flag[flag<0] = 0
    fnames = ['A', 'H', 'P', 'N', 'C', 'F', 'V']
    flist = []
    for i in range(len(fnames)):
        if fname.find(fnames[i]) >= 0:
            print fnames[i], i
            flist.append(i)
    l = (np.sum([flag == i for i in flist], axis=0)) == 1
    t = t[l]; t_err = t_err[l]
    p = p[l]; p_err = p_err[l]
    a = a[l]; a_errp = a_errp[l]; a_errm = a_errm[l]
    g = g[l]; g_errp = g_errp[l]; g_errm = g_errm[l]
    flag = flag[l]

    # reduce errorbars if they go below zero
    # (only necessary if you don't use asym)
    g_err = .5*(g_errp+g_errm)
    a_err = .5*(a_errp + a_errm)
    diff = a - a_err
    l = diff < 0
    a_err[l] = a_err[l] + diff[l] - np.finfo(float).eps

    # LOO
    if cv:
        print 'cv', cv
        tn = np.ones_like(a)
        select = cv
        print select
        tn[select] = 0
        tn = tn==1
        print len(p), len(p[tn])

    if cv:
        return a[tn], a_err[tn], a_errp[tn], a_errm[tn], p[tn], p_err[tn], \
                t[tn], t_err[tn], g[tn], g_err[tn], g_errp[tn], g_errm[tn], flag[tn]
    print len(p)

    return a, a_err, a_errp, a_errm, p, p_err, t, t_err, g, g_err, g_errp, g_errm, flag
