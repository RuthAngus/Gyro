import numpy as np
import matplotlib.pyplot as pl

# load Garcia data
data = np.genfromtxt("garcia.txt").T
print np.shape(data)
KID, p, p_err, t, t_err, feh, feh_err, e7, e8, e9, logg, \
        logg_errp, logg_errm, a, a_errp, a_errm, flag, e17, e18 = data
flag+=3
print len(KID)

# load my data
data = np.genfromtxt("all_astero.txt").T
print np.shape(data)
myKID, myt, myt_err, mya, mya_errp, mya_errm, myp, myp_err, mylogg, mylogg_errp, \
        mylogg_errm, myfeh, myfeh_err, myflag = data
print len(myKID)

# add extras
eKID=[]; et=[]; et_err=[]; ea=[]; ea_errp=[]; ea_errm=[]; ep=[]; ep_err=[]
elogg=[]; elogg_errp=[]; elogg_errm=[]; efeh=[]; efeh_err=[]; eflag=[]
for i, kid in enumerate(myKID):
    l = KID==kid
    if sum(l) == 0:
        eKID.append(myKID[i])
        et.append(myt[i])
        et_err.append(myt_err[i])
        ea.append(mya[i])
        ea_errp.append(mya_errp[i])
        ea_errm.append(mya_errm[i])
        ep.append(myp[i])
        ep_err.append(myp_err[i])
        elogg.append(mylogg[i])
        elogg_errp.append(mylogg_errp[i])
        elogg_errm.append(mylogg_errm[i])
        efeh.append(myfeh[i])
        efeh_err.append(myfeh_err[i])
        eflag.append(myflag[i])
KID = np.concatenate((KID, np.array(eKID)))
t = np.concatenate((t, np.array(et)))
t_err = np.concatenate((t_err, np.array(et_err)))
a = np.concatenate((a, np.array(ea)))
a_errp = np.concatenate((a_errp, np.array(ea_errp)))
a_errm = np.concatenate((a_errm, np.array(ea_errm)))
p = np.concatenate((p, np.array(ep)))
p_err = np.concatenate((p_err, np.array(ep_err)))
logg = np.concatenate((logg, np.array(elogg)))
logg_errp = np.concatenate((logg_errp, np.array(elogg_errp)))
logg_errm = np.concatenate((logg_errm, np.array(elogg_errm)))
feh = np.concatenate((feh, np.array(efeh)))
feh_err = np.concatenate((feh_err, np.array(efeh_err)))
flag = np.concatenate((flag, np.ones_like(eKID)))
print len(KID)

# I checked for duplicates and amy data already!

# load victor and travis data and replace
data = np.genfromtxt("vandt.txt", skip_header=1).T
n=0
for i, kid in enumerate(data[0]):
    l = KID==kid
    print KID[l]
    t[l] = data[1][i]
    t_err[l] = data[2][i]
    a[l] = data[3][i]
    a_errp[l] = data[4][i]
    a_errm[l] = data[5][i]
    p[l] = data[6][i]
    p_err[l] = data[7][i]
    logg[l] = data[8][i]
    logg_errp[l] = data[9][i]
    logg_errm[l] = data[10][i]
    feh[l] = data[11][i]
    feh_err[l] = data[12][i]
    flag[l] = 2
    n+=1
print n, 'n'

l = (a!=0) * (t!=0)
data = np.zeros((14, len(KID[l])))
data = KID[l], t[l], t_err[l], a[l], a_errp[l], a_errm[l], p[l], p_err[l], logg[l], \
        logg_errp[l], logg_errm[l], feh[l], feh_err[l], flag[l]
np.savetxt("garcia_all_astero.txt", data)

pl.clf()
pl.errorbar(a[l], p[l], xerr=(a_errp[l], a_errm[l]), yerr=p_err[l], capsize=0, fmt='k.')
pl.savefig("garcia_p_vs_a")

pl.clf()
pl.errorbar(t[l], p[l], xerr=t_err[l], yerr=p_err[l], capsize=0, fmt='k.')
pl.savefig("garcia_p_vs_t")
