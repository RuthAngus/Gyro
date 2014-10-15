import numpy as np
import matplotlib.pyplot as pl

# load Garcia data
data = np.genfromtxt("garcia.txt").T
print np.shape(data)
KID, p, p_err, t, t_err, feh, feh_err, e7, e8, e9, logg, \
        logg_errp, logg_errm, a, a_errp, a_errm, flag, e17, e18 = data
flag+=3

data = np.genfromtxt('/Users/angusr/Python/Gyro/data/Garcia_et_al_2014_Table_1.txt').T
KID = data[0]
p = data[1]
p_err = data[2]

sdss, sdss_err = 5, 6
irfm, irfm_err = 7, 8
teff, teff_err = irfm, irfm_err
# Load Astero data from table 1 - KIDs, teffs and feh
data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=30, skip_footer=1335, invalid_raise=False, usecols=(0,teff,teff_err,7,8,9,10)).T
KID1 = data1[0]
IRFM = data1[3]
# Load Astero data from table 1 - KIDs, teffs and feh
data1 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=30, skip_footer=1335, invalid_raise=False, usecols=(0,teff,teff_err,9,10)).T
table = np.zeros((11, len(KID)))
# Load Astero data from table 2 - KIDs, teffs and feh 5,6,7,8
data12 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=576, skip_footer=1222, invalid_raise=False, usecols=(0,teff,teff_err,7,8)).T

sf = 671
# load astero data from table 4 - KID, logg, age with SDSS teffs
data2 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=698, skip_footer=sf, invalid_raise=False, \
        usecols=(0,10,11,12,13,14,15)).T
KID2 = data2[0]
# load astero data from table 4 - KID, logg, age with SDSS teffs
data2 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=698, skip_footer=sf, invalid_raise=False, \
        usecols=(10,11,12,13,14,15)).T
# load astero data from table 5 - KID, logg, age
data3 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=1251, skip_footer=122, invalid_raise=False, \
        usecols=(0, 10,11,12,13,14,15)).T
KID3 = data3[0]
# load astero data from table 5 - KID, logg, age
data3 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=1251, skip_footer=122, invalid_raise=False, \
        usecols=(10,11,12,13,14,15)).T

# load astero data from table 6 - KID, logg, age
data4 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=1804, invalid_raise=False, \
        usecols=(0, 10,11,12,13,14,15)).T
KID4 = data4[0]
# load astero data from table 6 - KID, logg, age
data4 = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt', \
        skiprows=1804, invalid_raise=False, \
        usecols=(10,11,12,13,14,15)).T

# assemble data from 2 tables
for i, kid in enumerate(KID):
    l = KID1==kid
    if sum(l):
        table[:5].T[i] = data1.T[:][l]
#     else:
#         print kid
    ll = KID2==kid
    if sum(ll):
        table[5:].T[i] = data2.T[:][ll][0]

# replace missing SDSS temps with IRFM temps
for i in range(len(table.T)):
    if table[1][i]==0 and IRFM[i]!=0:
        table[1][i] = IRFM[i]

# replace values with Bruntt
for i, kid in enumerate(KID4):
    l = table[0]==kid
    if sum(l):
        table[1][l] = data12[1][i]
        table[2][l] = data12[2][i]
        table[3][l] = data12[3][i]
        table[4][l] = data12[4][i]

KID1, t, t_err, feh, feh_err, logg, logg_errp, logg_errm, a, a_errp, a_errm = table
print KID1[KID1==0]
print a[a==0]

# # load my data
# data = np.genfromtxt("all_astero.txt").T
# print np.shape(data)
# myKID, myt, myt_err, mya, mya_errp, mya_errm, myp, myp_err, mylogg, mylogg_errp, \
#         mylogg_errm, myfeh, myfeh_err, myflag = data
# print len(myKID)
#
# # add extras
# eKID=[]; et=[]; et_err=[]; ea=[]; ea_errp=[]; ea_errm=[]; ep=[]; ep_err=[]
# elogg=[]; elogg_errp=[]; elogg_errm=[]; efeh=[]; efeh_err=[]; eflag=[]
# for i, kid in enumerate(myKID):
#     l = KID==kid
#     if sum(l) == 0:
#         eKID.append(myKID[i])
#         et.append(myt[i])
#         et_err.append(myt_err[i])
#         ea.append(mya[i])
#         ea_errp.append(mya_errp[i])
#         ea_errm.append(mya_errm[i])
#         ep.append(myp[i])
#         ep_err.append(myp_err[i])
#         elogg.append(mylogg[i])
#         elogg_errp.append(mylogg_errp[i])
#         elogg_errm.append(mylogg_errm[i])
#         efeh.append(myfeh[i])
#         efeh_err.append(myfeh_err[i])
#         eflag.append(myflag[i])
# KID = np.concatenate((KID, np.array(eKID)))
# t = np.concatenate((t, np.array(et)))
# t_err = np.concatenate((t_err, np.array(et_err)))
# a = np.concatenate((a, np.array(ea)))
# a_errp = np.concatenate((a_errp, np.array(ea_errp)))
# a_errm = np.concatenate((a_errm, np.array(ea_errm)))
# p = np.concatenate((p, np.array(ep)))
# p_err = np.concatenate((p_err, np.array(ep_err)))
# logg = np.concatenate((logg, np.array(elogg)))
# logg_errp = np.concatenate((logg_errp, np.array(elogg_errp)))
# logg_errm = np.concatenate((logg_errm, np.array(elogg_errm)))
# feh = np.concatenate((feh, np.array(efeh)))
# feh_err = np.concatenate((feh_err, np.array(efeh_err)))
# flag = np.concatenate((flag, np.ones_like(eKID)))
# print len(KID)

# I checked for duplicates and amy data already!

# load victor and travis data and replace
data = np.genfromtxt("vandt.txt", skip_header=1).T
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

f = KID==0
print KID[f]
print len(KID)
l = (a!=0) * (t!=0)
ll = a==0
# print KID[ll], a[ll]
data = np.zeros((14, len(KID[l])))
data = KID[l], t[l], t_err[l], a[l], a_errp[l], a_errm[l], p[l], p_err[l], logg[l], \
        logg_errp[l], logg_errm[l], feh[l], feh_err[l], flag[l]
print len(KID[l])
# print data[ll]
# np.savetxt("garcia_all_astero.txt", data)
np.savetxt("garcia_irfm.txt", data)

pl.clf()
pl.errorbar(a[l], p[l], xerr=(a_errp[l], a_errm[l]), yerr=p_err[l], capsize=0, fmt='k.')
pl.savefig("garcia_p_vs_a")

pl.clf()
pl.errorbar(t[l], p[l], xerr=t_err[l], yerr=p_err[l], capsize=0, fmt='k.')
pl.savefig("garcia_p_vs_t")
