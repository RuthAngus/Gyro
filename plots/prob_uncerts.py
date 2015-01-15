import numpy as np
import matplotlib.pyplot as plt
import h5py

# takes filename and produces line parameters
def sigmas(fname):

    # load samples
    DIR = '/Users/angusr/Python/Gyro'
    with h5py.File("%s/code/samples_%s" % (DIR, fname), "r") as f:
        samples = f["samples"][:, 5000:, :]

    # take all "a" samples and reshape
    a_temp = samples[:, :, 0]
    a = np.reshape(a_temp, np.shape(a_temp)[0]*np.shape(a_temp)[1])

    # take all "b" samples and reshape
    b_temp = samples[:, :, 2]
    b = np.reshape(b_temp, np.shape(b_temp)[0]*np.shape(b_temp)[1])

    # take all "n" samples and reshape
    n_temp = samples[:, :, 1]
    n = np.reshape(n_temp, np.shape(n_temp)[0]*np.shape(n_temp)[1])

    # sample from a, b and n chains
    np.random.seed(123)
    r = np.random.randint(0, len(a), 100)
    samps = np.zeros((3, len(r)))
    for i in range(len(r)):
        samps[0, i] = a[r[i]]
        samps[2, i] = b[r[i]]
        samps[1, i] = n[r[i]]

    # output 2d array with parameter, samples
    return samps

if __name__ == "__main__":

    sigmas("ACHF")
