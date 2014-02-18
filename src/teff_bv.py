# Convert Teffs to B-V colours, or vice-versa
# using the sekiguchi and fukugita conversion.

import numpy as np
import matplotlib.pyplot as pl

def bv2teff(bv, logg, feh):
    # best fit parameters
    c = [3.929883, -0.360726, 0.168806, -0.048300]
    f = [0.025761, 0.004537]
    g1 = 0.007367
    h1 = -0.01069
    return c[0] + c[1]*bv + c[2]*bv**2 + c[3]*bv**3 + \
            f[0]*feh + f[1]*feh**2 + g1*logg + h1*bv*logg

def teff2bv(teff, logg, feh):
    # best fit parameters
    t = [-813.3175, 684.4585, -189.923, 17.40875]
    f = [1.2136, 0.0209]
    d1 = -0.294
    g1 = -1.166
    e1 = 0.3125
    return t[0] + t[1]*np.log10(teff) + t[2]*(np.log10(teff))**2 + \
            t[3]*(np.log10(teff))**3 + f[0]*feh + f[1]*feh**2 \
            + d1*feh*np.log10(teff) + g1*logg + e1*logg*np.log10(teff)

# load data
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/data.txt').T
teff = data[3]
logg = data[10]
feh = data[5]

bv = teff2bv(teff, logg, feh)

saveas = np.ndarray((len(data[0]), 2)).T
saveas[0,:] = data[0]
saveas[1,:] = bv

np.savetxt("/Users/angusr/Python/Gyro/data/colours.txt", saveas)

pl.clf()
pl.plot(teff, bv, 'k.')
pl.savefig("/Users/angusr/Python/Gyro/plots/teff_bv")
