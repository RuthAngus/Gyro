# Fit data to model in 3 dimensions where there are uncertainties in all.
# The way to do this is to draw n samples from a gaussian with mean x_i and stdv sigma_i.
#(basically pretend you have many more measurements than you actually do).
# Try it in 2-d first!

import numpy as np
import gyro_calibration
import random

# ''' Load data '''
# data = np.genfromtxt('/Users/angusr/Documents/rotation/final_catalogue.txt').T
# KIDs = data[0]
# measured_p = data[1]
# errors = data[2]
# B_Vs = data[3]
# # B_Vs = 0.5 #FIXME
# p0 = [0.5189, 0.7725, 0.601]

def totls2d():

    data = np.genfromtxt('/Users/angusr/Documents/rotation/final_catalogue.txt').T
    B_Vs = 0.5 #FIXME
    p0 = [0.5189, 0.7725, 0.601]
    
    x = gyro_calibration.Period_measurements.x_values()
    x_err = 0.1*np.ones(len(x))

    y = data[1]
    yerr = data[2]

    # for each x, draw from gausssian
    ndraws = 100
    for i in range(len(x)):
        for j in ndraws:
            new_x = random.gauss(x[i], x_err[i])
    
