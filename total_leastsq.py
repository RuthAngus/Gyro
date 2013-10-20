# Fit data to model in 3 dimensions where there are uncertainties in all.
# The way to do this is to draw n samples from a gaussian with mean x_i and stdv sigma_i.
#(basically pretend you have many more measurements than you actually do).
# Try it in 2-d first!

import numpy as np
import gyro_calibration as gc
import random
import pylab as p
from scipy.optimize import fmin

# Load data
data = np.genfromtxt('/Users/angusr/Documents/rotation/final_catalogue.txt').T
KIDs = data[0]
measured_p = data[1]
errors = data[2]
B_Vs = data[3]
B_Vs = 0.5 #FIXME

# Gyro parameters: n, a, b
p0 = np.array([0.5189, 0.7725, 0.601])


class totls2d(object):

    def __init__(self, x, x_err, y, y_err, p0):
        self.x = x
        self.x_err = x_err
        self.y = y
        self.y_err = y_err
        self.p0 = p0

    def minimise(self, func, p0 = None):
        if p0 == None:
            p0 = self.p0
        else:
            self.p0 = p0            

        # Optimise parameters
        new_p = fmin(func, p0)
        final_resid = self.residuals(new_p)
        print new_p
        self.params = new_p

        # Plotting
        x_plot = np.r_[x.min():x.max():101j]
        B_Vs = 0.5 # FIXME        
        y_plot = model_periods(p0, x_plot, B_Vs)
        new_params_y_plot = model_periods(new_p, x_plot, B_Vs)

        p.close(1)
        p.figure(1)
        # p.subplot(2,1,1)
        p.plot(x_plot, y_plot, 'r-')
        p.plot(x_plot, new_params_y_plot, 'b-')
        p.errorbar(self.x, self.y, yerr = self.y_err, xerr = self.x_err, fmt = 'k.')
        # p.subplot(2,1,2)
        # p.plot(x_plot, y_plot, 'r-')
        # p.plot(x_plot, new_params_y_plot, 'b-')
        # p.errorbar(new_x, new_y, yerr = new_yerr, fmt = 'k.')
        p.ylabel('Period (days)')
        p.xlabel('Age (Myr)')
        p.show()
    
    def residuals(self, params):
        x, y, yerr = self.x, self.y, self.y_err
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
        err = sum(((tot_y - model_periods(params, tot_x, B_Vs))/tot_yerr)**2)
        # print 'RESIDUALS...'
        # print err
        # #raw_input('Press return to continue')
        # print '...RESIDUALS'        
        return err

def model_periods(params, x, B_Vs):
    a = np.log10(x)
    b = np.log10(B_Vs - 0.4)
    c = params[0] * a + np.log10(params[1]) + params[2] * b
    # print 'MODEL_PERIODS...'
    # print params
    # print a.min(), a.max()
    # print b.min(), b.max()
    # print c.min(), c.max()
    # print '...MODEL_PERIODS'
    return 10.0**c

def measured_ages(params, x, B_Vs):
    ages = 10** ((1./params[0]) * (np.log10(x) - np.log10(params[1]) - \
params[2]*np.log10(B_Vs - 0.4)))
    
def model_ages():
    data = np.genfromtxt('/Desktop/aster_ages.txt')
    kicic = data[0]
    ages = data[-1]
    print ages
    return ages

# Define x
B_Vs = 0.5 # FIXME
gyro_cal = gc.Period_measurements(KIDs, measured_p, errors, B_Vs, p0)
x = gyro_cal.x_values()
x_err = 0.05*x # make up uncertainties
p_true = np.array([0.56, 0.72, 0.53])

F = totls2d(x, x_err, measured_p, errors, p0)
p.clf()
p.errorbar(F.x, F.y, xerr = F.x_err, yerr = F.y_err, fmt = 'ko')
p.plot(F.x, model_periods(p0, F.x, B_Vs), 'r+')
p.plot(F.x, model_periods(p_true, F.x, B_Vs), 'cx')
#print F.residuals(p0)
F.minimise(F.residuals)
                          
                          
                          
                      
                          
                          
