# Fit data to model in 3 dimensions where there are uncertainties in all.
# The way to do this is to draw n samples from a gaussian with mean x_i and stdv sigma_i.
#(basically pretend you have many more measurements than you actually do).
# Try it in 2-d first!

import numpy as np
import gyro_calibration as gc
import random
import pylab as p
from scipy.optimize import fmin


''' Load data '''
data = np.genfromtxt('/Users/angusr/Documents/rotation/final_catalogue.txt').T
KIDs = data[0]
measured_p = data[1]
errors = data[2]
B_Vs = data[3]
B_Vs = 0.5 #FIXME
p0 = [0.5189, 0.7725, 0.601]


class totls2d(object):

    def __init__(self, x, x_err, y, y_err, p0):
        self.x = x
        self.x_err = x_err
        self.y = y
        self.y_err = y_err
        self.p0 = p0

    def minimise(self):

        # Optimise parameters
        new_p = fmin(gyro_cal.residuals, p0, args = (self.y, self.x))
        print new_p

        # Plotting
        x_plot = np.arange(0, max(self.x), 10)
        B_Vs = 0.5 # FIXME        
        y_plot = 10**(self.p0[0]*(np.log10(x_plot)) + np.log10(self.p0[1]) + \
                                 self.p0[2]*(np.log10(B_Vs - 0.4)))
        new_params_y_plot = 10**(new_p[0]*(np.log10(x_plot)) + np.log10(new_p[1]) + \
                                 new_p[2]*(np.log10(B_Vs - 0.4)))
        
        p.close(1)
        p.figure(1)
        # p.subplot(2,1,1)
        p.plot(x_plot, y_plot, 'r-')
        p.plot(x_plot, new_params_y_plot, 'b-')
        p.errorbar(self.x, self.y, yerr = self.y_err, fmt = 'k.')
        # p.subplot(2,1,2)
        # p.plot(x_plot, y_plot, 'r-')
        # p.plot(x_plot, new_params_y_plot, 'b-')
        # p.errorbar(new_x, new_y, yerr = new_yerr, fmt = 'k.')
        p.ylabel('Period (days)')
        p.xlabel('Age (Myr)')
        p.show()
    
        def residuals(self, p, y, x):

            # for each x, draw from gausssian and assign the same y value to each 
            ndraws = 100
            tot_x = np.zeros(ndraws*len(x))
            tot_y = np.zeros(ndraws*len(x))
            tot_yerr = np.zeros(ndraws*len(x))
            n = 0
            for i in range(len(x)):
                for j in range(ndraws):
                    tot_x[n] = random.gauss(x[i], x_err[i]) 
                    tot_y[n] = y[i]
                    tot_yerr[n] = yerr[i]
                    n += 1
                    
            B_Vs = 0.5 # FIXME
            model_periods = 10**(p[0]*(np.log10(tot_x)) + np.log10(p[1]) + \
                                 p[2]*(np.log10(B_Vs - 0.4)))
            err = sum((tot_y - model_periods(x))/tot_yerr)
            return err


# Define x
B_Vs = 0.5 # FIXME
gyro_cal = gc.Period_measurements(KIDs, measured_p, errors, B_Vs, p0)
x = gyro_cal.x_values()
x_err = 0.35*np.ones(len(x)) # make up uncertainties

find_params = totls2d(x, x_err, measured_p, errors, p0)
find_params.minimise()
