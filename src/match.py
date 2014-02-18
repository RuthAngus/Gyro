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

# # Assemble data from various places and concatenate
# data = np.genfromtxt('/Users/angusr/Python/Gyro/all_data.txt').T
# l = data[0] != 10124866.0 # Remove missing data point (now only have 42)
# KID = data[0,l]
# period = data[1,l]
# period_err = data[2,l]
#
# # Load Astero data
# data = np.genfromtxt('/Users/angusr/Desktop/astero_ages.txt').T
