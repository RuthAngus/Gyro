import numpy as np
import matplotlib.pyplot as pl
from teff_bv import teff2bv_orig, teff2bv_err
import pretty5
import sys
from colours import plot_colours
ocols = plot_colours()

plotpar = {'axes.labelsize': 18, 'text.fontsize': 10,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)

# load data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/garcia_irfm.txt")
subgiant = 4.2
g = data[8] > subgiant
p1 = data[6][g]
p_err1 = data[7][g]
t1 = data[1][g]
t_err1 = data[2][g]
a1 = data[3][g]
a_errp1 = data[4][g]
a_errm1 = data[5][g]
a_err1 = .5 * (a_errp1 + a_errm1)
logg1 = data[8][g]
logg_err1 = .5*(data[9][g]+data[10][g])
feh1 = data[11][g]
feh_err1 = data[12][g]
flag1 = data[13][g]
bv1, bv_err1 = teff2bv_err(t1, logg1, feh1, t_err1, logg_err1, feh_err1)
# add clusters
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/clusters.txt").T
l = data[4] > 0
bv2 = data[0][l]; bv_err2 = data[1][l]
p2 = data[2][l]; p_err2 = data[3][l]
a2 = data[4][l]; a_err2 = data[5][l]; a_errp2 = data[5][l]; a_errm2 = data[5][l]
logg2 = data[6][l]; logg_err2 = data[7][l]
flag2 = data[8][l]
# combine astero and cluster
bv = np.concatenate((bv1, bv2))
bv_err = np.concatenate((bv_err1, bv_err2))
p = np.concatenate((p1, p2))
p_err = np.concatenate((p_err1, p_err2))
a = np.concatenate((a1, a2))
a_err = np.concatenate((a_err1, a_err2))
a_errp = np.concatenate((a_errp1, a_errp2))
a_errm = np.concatenate((a_errm1, a_errm2))
logg = np.concatenate((logg1, logg2))
logg_err = np.concatenate((logg_err1, logg_err2))
flag = np.concatenate((flag1, flag2))
# turn astero target flags to '3' and take 3 away
flag[flag==2] = 9
flag -= 3
flag[flag<0] = 0

age = 4.568
fnames = ['A', 'H', 'P', 'N', 'C', 'F', 'V']
fname = "HCNPF45"

# select star group
flist = []
for i in range(len(fnames)):
    if fname.find(fnames[i]) >= 0:
        print fnames[i], 'yes'
        flist.append(i)
l = (np.sum([flag == i for i in flist], axis=0)) == 1

pl.clf()
bv -= 0.45
hl = a[l]==.625
cl = a[l]==.5
prl = a[l]==.588
nl = a[l]==1.1
ms = 6

pl.errorbar(bv[l][prl], p[l][prl], xerr=bv_err[l][prl], yerr=p_err[l][prl],
        fmt='^', color=ocols.orange, capsize=0, ecolor='.7', markersize=ms,
        label="$\mathrm{Praesepe}$", mec=ocols.orange)
pl.errorbar(bv[l][hl], p[l][hl], xerr=bv_err[l][hl], yerr=p_err[l][hl],
        fmt='o', color=ocols.blue, capsize=0, ecolor='.7', markersize=ms,
        label="$\mathrm{Hyades}$", mec=ocols.blue)
pl.errorbar(bv[l][cl], p[l][cl], xerr=bv_err[l][cl], yerr=p_err[l][cl],
        fmt='o', color=ocols.maroon, capsize=0, ecolor='.7', markersize=ms,
        label="$\mathrm{Coma~Ber}$", mec=ocols.maroon)
pl.errorbar(bv[l][nl], p[l][nl], xerr=bv_err[l][nl], yerr=p_err[l][nl],
        fmt='^', color=ocols.green, capsize=0, ecolor='.7', markersize=ms,
        label="$\mathrm{NGC~6811}$", mec=ocols.green)
pl.xlabel('$\mathrm{B-V-}~c$')
pl.ylabel('$\mathrm{Period~(days)}$')
pl.subplots_adjust(hspace=.3)
pl.xlim(10**-1.6, 10**0.1)
pl.ylim(10**.6, 10**1.3)
pl.loglog()
pl.legend(loc="best")
pl.savefig("/Users/angusr/Python/Gyro/gyro_paper/show%s.pdf"%fname)
