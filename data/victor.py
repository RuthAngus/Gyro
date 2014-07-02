import numpy as np
from match import match
import subprocess

ldata = np.genfromtxt("/Users/angusr/Python/Gyro/data/Victor_params.txt", skip_header=1).T
# pdata = np.genfromtxt("/Users/angusr/Python/Gyro/data/victor_periods.txt").T
# data = np.zeros((len(ldata[0]), 19))
# data[:,0] = ldata[0] # KID
# data[:,1] = pdata[1] # period
# data[:,2] = pdata[2] # p_err
# data[:,3] = ldata[3] # teff
# data[:,4] = ldata[4] # t_err
# data[:,5] = ldata[1] # feh
# data[:,6] = ldata[2] # feh_err
# data[:,7] = ldata[8] # mass
# data[:,8] = ldata[9] # mass_errp
# data[:,9] = ldata[10] # mass_errm
# data[:,10] = ldata[11] # logg
# data[:,11] = ldata[12] # logg_errp
# data[:,12] = ldata[13] # logg_errm
# data[:,13] = ldata[14] # age
# data[:,14] = ldata[15] # age_errp
# data[:,15] = ldata[16] # age_errm
#
# KID = ldata[0]
# period = pdata[1]
# period_err = pdata[2]
#
# # load ACF results
# data2 = np.genfromtxt("/Users/angusr/angusr/ACF2/results.txt").T
# KID2 = data2[0]
# period2 = data2[1]
# period_err2 = data2[2]
#
# data3 = match(KID, data2)
# KID3 = data3[0]
# period3 = data3[1]
# period_err3 = data3[2]
#
# for i in range(len(KID)):
#     if period_err[i] > period_err3[i]:
#         print KID[i], KID3[i]
#         print period[i], period3[i]
#         print period_err[i], period_err3[i], '\n'
#         raw_input('enter')
#     else:
#         data[2][i] = period_err3[i]
#
# np.savetxt("/Users/angusr/Python/Gyro/data/victor_p_errs.txt", data)
#
# print data

KID = ldata[0]
for k, i in enumerate(KID):
    print int(i)
#     subprocess.call("open ~/angusr/ACF/PDCQ3_output/plots_acf/%s_full.png"\
#             %(int(i)), shell=True)
#     subprocess.call("open ~/angusr/ACF/PDCQ4_output/plots_acf/%s_full.png"\
#             %(int(i)), shell=True)
#     subprocess.call("open ~/angusr/ACF/PDCQ5_output/plots_acf/%s_full.png"\
#             %(int(i)), shell=True)
#     subprocess.call("open ~/angusr/ACF/PDCQ6_output/plots_acf/%s_full.png"\
#             %(int(i)), shell=True)
    data = np.genfromtxt('/Users/angusr/angusr/ACF2/%s_all_result.txt'%int(KID[k])).T
    print data
    data = np.genfromtxt('/Users/angusr/Python/Gyro/data/victor_periods.txt').T
    print data[1][k], data[2][k]
    raw_input('enter')

# data = np.genfromtxt('/Users/angusr/Python/Gyro/data/travis_data.txt', skip_header=1).T
# KID = data[0]
#
# periods = np.empty((len(KID), 3))
# for i in range(len(KID)):
#     periods[i][0] = KID[i]
#     data = np.genfromtxt('/Users/angusr/angusr/ACF2/%s_all_result.txt'%int(KID[i])).T
#     print data[1]
#     periods[i][1:] = data
#
# np.savetxt('/Users/angusr/Python/Gyro/data/travis_actual_periods.txt', periods)
