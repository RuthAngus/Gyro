# This function takes in a list of KIDs and an ndarray of KIDs and variables
# (the KIDs must be the first column)
# It will map the variables to the KIDs and output an ndarray of the reassigned
# variables in the same order as the original KID list.

import numpy as np

def match(KID, X):
    n = len(KID)
    m, l = X.shape
    KID_obs = X[0,:]
    Xm = np.zeros((m, n))
    for i in range(n):
        l = np.where(KID_obs == KID[i])[0]
        if len(l) > 0:
            Xm[:,i] = X[:,l].reshape(Xm[:,i].shape)
    return Xm

def assemble():

    # load period data
    pdata = np.genfromtxt('/Users/angusr/Python/Gyro/data/all_data.txt').T
    KID = pdata[0]
    p = pdata[1]
    p_err = pdata[2]
    print KID


    # Load Astero data from table 1 - KIDs, teffs and feh
    # Columns: KID, nu, nu_err, dnu, dnu_err, SDSS_teff, st_err,
    # IRFM_teff, It_err, feh, feh_err
    data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJS91604R2tables.txt', \
            skiprows=30, skip_footer=1343, invalid_raise=False, usecols=(0,5,6,9,10)).T
    table1 = match(KID, data1)

    # load astero data from table 5 - KID, mass, logg, age
    data2 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJS91604R2tables.txt', \
            skiprows=699, skip_footer=675, invalid_raise=False, \
            usecols=(0,1,2,3,11,12,13,14,15)).T
    table5 = match(KID, data2)

    # Assemble data in the following order:
    # KID, period, p_err, teff, teff_err, feh, feh_err, mass, mass_errp,
    # mass_errm, logg, logg_errp, logg_errm, age, age_errp, age_errm
#     data = np.ndarray((len(KID), 16))
    data = np.zeros((len(KID), 16))
    data[:,0] = KID
    data[:,1] = p
    data[:,2] = p_err
    data[:,3:7] = table1[1:,:].T
    data[:,8:16] = table5[1:,:].T

    return data

data = assemble()
np.savetxt("/Users/angusr/Python/Gyro/data/data.txt", data)
