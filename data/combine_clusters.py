import numpy as np
import matplotlib.pyplot as pl
from colour_conversion import gr2bv

c_err = .01

# add hyades
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/hyades.txt", skip_header=2).T
bv = data[0]
bv_err = data[1]
p = data[2]
p_err = data[3]
a = data[4]
a_err = data[5]
a_errp = data[5]
a_errm = data[5]

# add praesepe
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/praesepe.txt').T
bv = np.concatenate((bv, (data[5]-data[6])))
bv_err = np.concatenate((bv_err, np.ones_like(data[5])*c_err))
p = np.concatenate((p, 1./data[3]))
p_err = np.concatenate((p_err, (1./data[3])*(data[4]/data[3])))
# a = np.concatenate((a, np.ones_like(data[5])*.59))
a = np.concatenate((a, np.ones_like(data[5])*.588))
a_err = np.concatenate((a_err, np.ones_like(data[5])*.137))
a_errp = np.concatenate((a_errp, np.ones_like(data[5])*.137))
a_errm = np.concatenate((a_errm, np.ones_like(data[5])*.137))

# add NGC6811
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/NGC6811", skip_header=2).T
p = np.concatenate((p, data[3]))
p_err = np.concatenate((p_err, data[4]))
mbv, mbv_err = gr2bv(data[1], data[2])
bv = np.concatenate((bv, mbv))
bv_err = np.concatenate((bv_err, mbv_err))
a = np.concatenate((a, np.ones_like(data[3])*1.1))
a_err = np.concatenate((a_err, np.ones_like(data[3])*.1))
# a_err = np.concatenate((a_err, np.ones_like(data[3])*.2))
pl.clf()
pl.plot(mbv, data[3], 'ko')
pl.savefig('NGC6811')

logg = np.ones_like(bv)*4.5
logg_err = np.ones_like(bv)*.01

data = np.zeros((8, len(bv)+1))
data[0,:-1] = bv
data[1,:-1] = bv_err
data[2,:-1] = p
data[3,:-1] = p_err
data[4,:-1] = a
data[5,:-1] = a_err
data[6,:-1] = logg
data[7,:-1] = logg_err

# add the sun
datasun = list(np.genfromtxt("/Users/angusr/Python/Gyro/data/sun.txt", skip_header=2).T)
datasun.append(4.5)
datasun.append(logg_err[-1])
datasun = np.array(datasun)
data[:,-1] = datasun

np.savetxt("clusters.txt", np.transpose(data))
