# This function takes in a list of KIDs and an ndarray of KIDs and variables
# (the KIDs must be the first column)
# It will map the variables to the KIDs and output an ndarray of the reassigned
# variables in the same order as the original KID list.

import numpy as np
from find_extra_amy import extra

# match takes in the list of KIDs that you want to map the values of the matrix X to.
# the first column of X should also be a list of KIDs
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

def assemble(KID, p, p_err):

#     # find and load extra amy data
#     KID2, p2, p_err2 = extra(KID)
#     KID = np.concatenate((KID, KID2))
#     p = np.concatenate((p, p2))
#     p_err = np.concatenate((p_err, p2))

    # Load Astero data from table 1 - KIDs, teffs and feh
    # Columns: KID, nu, nu_err, dnu, dnu_err, SDSS_teff, st_err,
    # IRFM_teff, It_err, feh, feh_err
#     data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJS91604R2tables.txt', \
    data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/filled_in.txt', \
            skiprows=30, skip_footer=1343, invalid_raise=False, usecols=(0,5,6,9,10)).T
    table1 = match(KID, data1)
    dnu_data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJS91604R2tables.txt', \
            skiprows=30, skip_footer=1343, invalid_raise=False, usecols=(0,3,4)).T
    dnu_table1 = match(KID, dnu_data1)

    # load astero data from table 5 - KID, mass, logg, age
    data2 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJS91604R2tables.txt', \
            skiprows=699, skip_footer=675, invalid_raise=False, \
            usecols=(0,1,2,3,10,11,12,13,14,15)).T
    table5 = match(KID, data2)

    # Assemble data in the following order:
    # [0]KID, [1]period, [2]p_err, [3]teff, [4]teff_err, [5]feh, [6]feh_err, [7]mass, [8]mass_errp,
    # [9]mass_errm, [10]logg, [11]logg_errp, [12]logg_errm, [13]age, [14]age_errp, [15]age_errm, [17]dnu, [18]dnu_err
#     data = np.ndarray((len(KID), 16))
    data = np.zeros((len(KID), 19))
    data[:,0] = KID
    data[:,1] = p
    data[:,2] = p_err
    data[:,3:7] = table1[1:,:].T
    data[:,7:16] = table5[1:,:].T
    data[:,17] = dnu_table1[1]
    data[:,18] = dnu_table1[2]

    return data

if __name__ == "__main__":

#     # load period data
#     pdata = np.genfromtxt('/Users/angusr/Python/Gyro/data/all_data.txt').T
# #     pdata = np.genfromtxt('/Users/angusr/angusr/ACF2/periods.txt').T
#     data = np.genfromtxt('/Users/angusr/Python/Gyro/data/matched_data.txt').T
#     data = np.genfromtxt('/Users/angusr/Python/Gyro/data/extra_amy.txt').T
    data = np.genfromtxt('/Users/angusr/Python/Gyro/data/Garcia_et_al_2014_Table_1.txt').T
    KID = data[0]
    p = data[1]
    p_err = data[2]
    data = assemble(KID, p, p_err)

#     np.savetxt("/Users/angusr/Python/Gyro/data/new_data.txt", data)
#     np.savetxt("/Users/angusr/Python/Gyro/data/old_data.txt", data) # should be matched_data in the proper format
#     np.savetxt("/Users/angusr/Python/Gyro/data/recovered.txt", data) # this was supposed to be
#     np.savetxt("/Users/angusr/Python/Gyro/data/extra_amy_matched.txt", data)
#     np.savetxt("/Users/angusr/Python/Gyro/data/new_matched.txt", data)
    np.savetxt("/Users/angusr/Python/Gyro/data/garcia.txt", data)
