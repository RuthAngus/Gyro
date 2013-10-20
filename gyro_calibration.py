
import numpy as np
import pylab as p
from scipy.optimize import fmin
import mpfit


''' Load data '''
data = np.genfromtxt('/Users/angusr/Documents/rotation/final_catalogue.txt').T
KIDs = data[0]
measured_p = data[1]
errors = data[2]
B_Vs = data[3]
# B_Vs = 0.5 #FIXME
p0 = [0.5189, 0.7725, 0.601]

class Period_measurements(object):

     def __init__(self, KIDs, measured_p, errors, B_Vs, p0):
        self.KIDs = KIDs
        self.measured_p = measured_p
        self.errors = errors
        self.B_Vs = B_Vs
        self.p0 = p0

     def minimise(self, doplot = True):

        nf = 0.56; af = 0.72; bf = 0.53 # Don't worry about these - just here for fake ages
        model_ages = self.period_age(nf, af, bf)        

        ''' scipy.optimize.leastsq '''

        y = self.measured_p
        x = model_ages
        result = fmin(self.residuals, self.p0, args = (y, x))
        print result
            
     def residuals(self, p, y, x):
        # y = self.measured_p, x = model_ages
        #n, a, b = p
        #model_p = self.model_periods(p0)
        err = sum((y - self.model_periods(x))/self.errors)
        return err

     def peval(self, x, p):
        return self.model_periods(x)

     def model_periods(self, x):
        return self.age_period(x)

     def period_age(self, n, a, b):
        age = 10**((1./n)*(np.log10(self.measured_p)-np.log10(a)-\
                                    b*(np.log10(self.B_Vs-0.4))))
        x = self.B_Vs - 0.4
        errors = 0.02*(np.sqrt(3.0 + 0.5*np.log(age)**2.0+2.0*self.measured_p**0.6+\
                (0.6/x)**2+(2.4*np.log(x))**2)) # FIXME
        return age
 
     def age_period(self, age):
        return 10**(self.p0[0]*(np.log10(age)) + np.log10(self.p0[1]) + \
                    self.p0[2]*(np.log10(self.B_Vs - 0.4)))

     def x_values(self):
          # nf = 0.5189; af = 0.7725; bf = 0.601
          nf = 0.56; af = 0.72; bf = 0.53 # Don't worry about these - just here for fake ages
          return self.period_age(nf, af, bf)   

find_params = Period_measurements(KIDs, measured_p, errors, B_Vs, p0)
find_params.minimise()


def model_ages():
    data = np.genfromtxt('/Desktop/aster_ages.txt')
    kicic = data[0]
    ages = data[-1]
    print ages
    return ages
