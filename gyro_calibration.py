# This code is for dumb-ass least-squares fitting of my rotation periods to the # asteroseismic
# ages and masses!

import numpy as np
import pylab as p
from scipy.optimize import leastsq
import mpfit

''' Load data '''
data = np.genfromtxt('/Users/angusr/Documents/rotation/final_catalogue.txt').T
KIDs = data[0]
measured_p = data[1]
errors = data[2]
B_Vs = data[3]
# B_Vs[-5] = 0.4, this is producing a divide by zero

class Period_measurements(object):

    def __init__(self, KIDs, measured_p, errors, B_Vs):
        self.KIDs = KIDs
        self.measured_p = measured_p
        self.errors = errors
        self.B_Vs = B_Vs

    def residuals(self, p, model_p):
        #n, a, b = p
        #model_p = self.model_periods(p0)
        err = self.measured_p - model_p
        return err

    def peval(self, p):
        return self.model_periods(p[0], p[1], p[2])

    def model_ages(self):
        nf = 0.5189; af = 0.7725; bf = 0.601 # Don't worry about these - just here for fake ages
        return 1.2*(self.convert_Barnes(nf, af, bf) + 0.2) # FIXME: actually want to change the params
        
    def model_periods(self, n, a, b):
        return self.age_period(self.model_ages(), n, a, b)

    def minimise(self):
        
        ''' scipy.optimize.leastsq '''

        p0 = [0.5189, 0.7725, 0.601]

        # This should be a scalar value
        print type(self.residuals(p0, self.model_periods(p0[0], p0[1], p0[2])))
        print type(p0)
        print type(self.measured_p)
        print type(self.model_ages())

        plsq = leastsq(self.residuals(p0, self.model_periods(p0[0], p0[1], p0[2])), p0, \
                       args = (self.measured_p, self.model_ages()))
        print plsq[0]

        # x0 = np.array[0.5189, 0.7725, 0.601]
        # res = scipy.optimize.minimize(convert_Barnes(

        ''' mpfit '''

        # # A total of 3 parameters, with starting values of 0.5189, 0.7725 and 0.601 are given.
        # # all constrained to be between 0 and 1
        
        # parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}]*3
        # parinfo[0]['limited'][0] = 1; parinfo[0]['limited'][1] = 1
        # parinfo[1]['limited'][0] = 1; parinfo[1]['limited'][1] = 1
        # parinfo[2]['limited'][0] = 1; parinfo[2]['limited'][1] = 1
        # parinfo[0]['limits'][0]  = 0.; parinfo[0]['limits'][1]  = 1.
        # parinfo[1]['limits'][0]  = 0.; parinfo[1]['limits'][1]  = 1.  # FIXME: do I need to specify B_Vs as a param
        # parinfo[2]['limits'][0]  = 0.; parinfo[2]['limits'][1]  = 1.  # to be held fixed?
        # values = [0.5189, 0.7725, 0.601]
        # for i in range(3): parinfo[i]['value']=values[i]

        
        # x = self.measured_p
        # p0 = [0.5189, 0.7725, 0.601]
        # y = 10**((1./p0[0])*(np.log10(measured_p)-np.log10(p0[1])- p0[2]*(np.log10(self.B_Vs-0.4))))
        # fa = {'x':x, 'y':y, 'err':self.errors}
        # m = mpfit.mpfit('myfunct', p0, functkw=fa)
        # print 'status = ', m.status
        # if (m.status <= 0): print 'error message = ', m.errmsg
        # print 'parameters = ', m.params
 

    def calibrate(self):

        n = 0.5189; a = 0.7725; b = 0.601
        
        # Calculate ages from measured periods with Barnes 2007
        ages = Period_measurements.convert_Barnes(self, n, a, b, doplot = False)

        model_ages = 1.2*(ages + 0.2)

        # Calculate the periods these fake ages would give you with Barnes 2007
        # Imagine that the self.B_Vs colours are the asteroseismic colours (v. small uncertainties)
        n = 0.5189; a = 0.7725; b = 0.601
        model_p = Period_measurements.age_period(self, model_ages, n, a, b)

        p.close(2)
        p.figure(2)
        p.errorbar(self.B_Vs, self.measured_p, yerr = errors, fmt = 'k.')
        p.plot(self.B_Vs, model_p, 'r.')

        dof = len(data) - 3 - 1
        # chi_sq = Period_measurements.chi_squared(self, measured_p, model_p, errors, dof)
        # print 'Chi squared =  ', chi_sq

    def convert_Barnes(self, n, a, b, doplot = False):
        
        
        # Calculating Barnes ages
        ages = np.zeros(len(self.KIDs))
        ages_Gyr = np.zeros(len(self.KIDs))
        errors = np.zeros(len(self.KIDs))

        # print self.B_Vs
        # print self.B_Vs-0.4
        log_age = (1./n)*(np.log10(measured_p)-np.log10(a)- b*(np.log10(self.B_Vs-0.4)))
        ages = 10**log_age
        ages_Gyr = ages/10.**3
        x = self.B_Vs - 0.4
        errors = 0.02*(np.sqrt(3.0 + 0.5*np.log(ages)**2.0+2.0*measured_p**0.6+\
                (0.6/x)**2+(2.4*np.log(x))**2)) # FIXME
            

        if doplot == True:
            #Calculating isochrones
            age_labels = ['100 Myr', '500 Myr', '1 Gyr', '2 Gyr', '3 Gyr','5 Gyr']
            cols = ['#6666ff','#00cc99','#66cc00','#ff9900','#ff0000', '#ff0099']
            p.close(1)
            p.figure(1)
            #ages.remove(ages[-5])
            iso_bv = np.linspace(min(self.B_Vs),max(self.B_Vs)+0.01, num=1000)
            iso_ages = [100, 500, 1000, 2000, 3000, 5000]
            iso_periods = np.zeros(len(iso_bv))
            for i in range(len(iso_ages)):
                for j in range(len(iso_periods)):
                    iso_periods[j] = 0.7725 * ((iso_bv[j] - 0.4)**0.601) \
                        * (iso_ages[i]**0.5189)
                p.plot(iso_bv, iso_periods, cols[i])
                p.text(0.76, iso_periods[800]-1.5, age_labels[i], color = cols[i])
            p.errorbar(self.B_Vs, self.measured_p, yerr = errors, fmt = 'k.')
            p.xlabel('B-V')
            p.ylabel('Rotation Period (days)')
            p.xlim(min(self.B_Vs)-0.01, max(self.B_Vs)+0.01)
            p.savefig('/Users/angusr/Documents/rotation/Gyro_results')
            p.savefig('/Users/angusr/Documents/rotation/Gyro_results.ps')
            
        return ages

    def age_period(self, ages, n, a, b):
        log_periods = n*(np.log10(ages)) + np.log10(a) + b*(np.log10(self.B_Vs - 0.4)) 
        periods = 10**log_periods
        return periods

    def chi_squared(self, data, model, uncertainties, dof):

        
        chi_sq = (1./float(dof))* sum( (data - model)**2 / uncertainties**2 )
        return chi_sq

    # # params is an array containing parameter values
    # def minimise_chi_sq(self, data, model, uncertainties, params):

        

    #     dof = len(data) - len(params) - 1
    #     chi_sq = (1./float(dof))* sum( (data - model)**2 / uncertainties**2 )

    #     # function
        
        
        
        
        

dataset = Period_measurements(KIDs, measured_p, errors, B_Vs)
dataset.calibrate()
