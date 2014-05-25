# This version plots different groups of coeval stars

import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as sci
import teff_bv

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 15,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#99CC99', '#9933FF', '#CC0000']

from matplotlib import rc
rc("font", size=20, family="serif", serif="Computer Sans")
rc("text", usetex=True)

# class plotting(object):
class plotting(object):

    def __init__(self, stars):
        self.KID = stars[0]
        self.period = stars[1]
        self.p_err = stars[2]
        self.teff = stars[3]
        self.t_err = stars[4]
        self.feh = stars[5]
        self.feh_err = stars[6]
        self.logg = stars[10]
        self.l_errp = stars[11]
        self.l_errm = stars[12]
        self.age = stars[13]
        self.a_errp = stars[14]
        self.a_errm = stars[15]

    # note - matched_data.txt has
    # [0]KID, [1]period, [2]period_err, [3]age, [4]age_errp, [5]age_errm, [6]mass, [7]mass_errp,
    # [8]mass_errm, [9]logg, [10]logg_errp, [11]logg_errm, [12]teff, [13]teff_err, [14]empty

#         self.KID = stars[0]
#         self.period = stars[1]
#         self.p_err = stars[2]
#         self.age = stars[3]
#         self.a_errp = stars[4]
#         self.a_errm = stars[5]
#         self.logg = stars[9]
#         self.l_errp = stars[10]
#         self.l_errm = stars[11]
#         self.teff = stars[12]
#         self.t_err = stars[13]
#         self.feh = stars[14]-.2

    def log_period_model(self, par, log_age, teff, logg, col):
        if col==False:
            return par[0] + par[1] * log_age + par[2] * np.log10(6250 - teff) # temp
#         bv = teff_bv.teff2bv(teff, logg, np.ones_like(teff)*-.2, error=False)
        return par[0] + par[1] * log_age + par[2] * np.log10(teff - .4) # colour

    def iso_calc(self, ag, pars, age, model, col):

        if type(ag) != np.ndarray: a_ref = ag
        else: a_ref = 10.0**(np.median(np.log10(age[ag])))

        bv_or_t = np.linspace(1.3, .4, 100)
        if col==False: bv_or_t = np.linspace(4000, 6249, 100)

        logg = np.ones_like(bv_or_t)*4.5
        period = 10**model(pars, np.log10(a_ref*1000), bv_or_t, logg, col)
        return bv_or_t, period, a_ref

    # remove kraft break and subgiant stars
    def kraft_sub(self, kraft):
        k = self.logg>4.0
#         if kraft==True: k = (self.mass < 1.4)*(self.logg>4.0)
        if kraft==True: k = (self.teff < 6250)*(self.logg>4.0)
        return self.period[k], self.p_err[k], self.age[k], self.teff[k], self.t_err[k], \
                self.a_errp[k], self.a_errm[k], self.logg[k], self.l_errp[k], \
                self.l_errm[k]

    def p_vs_t(self, pars, model, col):

        # convert teff data to colours
        if col==True:
            self.teff = teff_bv.teff2bv(self.teff, self.logg, self.feh)
            self.t_err = np.ones_like(self.teff)*0.01

        # remove subgiants and kraft break with kraft toggle
        period, p_err, age, teff, t_err, a_errp, a_errm, logg, logg_errp,\
                logg_errm = self.kraft_sub(kraft=False)

        # Remove zeros
        zero = teff > 0

        # age masks
        a_lim = [0, 1, 2, 3, 4, 5, 8, max(age)]

        # bin according to age
        for i in range(len(a_lim)-1):
            a = (a_lim[i] < age)*(age < a_lim[i+1])*zero

            # Plot data
            pl.clf()
            pl.errorbar(teff[a], period[a], xerr=t_err[a], \
                    yerr=p_err[a], color='k', fmt='o', mec='k', \
                    capsize=0, markersize=5, ecolor='0.8', zorder = 4, \
                    label='$%s < \mathrm{Age} <%s \mathrm{Gyr}$'%(a_lim[i], a_lim[i+1]))

            # Add Isochrones
            xs, ys, t_ref = plotting.iso_calc(self, a, pars, age, model, col)
            pl.plot(xs, ys, color = ocols[i], linestyle='-', linewidth = 2, \
                    label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ \
                    $\mathrm{(M\&H~2008)}$' % t_ref, zorder = 1)
            a = a_lim[i]
            xs, ys, t_ref = plotting.iso_calc(self, a, pars, age, model, col)
            pl.plot(xs, ys, color = ocols[i], linestyle = '--', linewidth = 2, zorder = 2)
            a = a_lim[i+1]
            xs, ys, t_ref = plotting.iso_calc(self, a, pars, age, model, col)
            pl.plot(xs, ys, color = ocols[i], linestyle = '--', linewidth = 2, zorder = 3)

            pl.xlabel("$\mathrm{B-V}$")
            if col==False: pl.xlabel("$\mathrm{T_{eff (K)}}$")
            pl.ylabel("$\mathrm{P_{rot} (days)}$")
    #             pl.xlim(pl.gca().get_xlim()[::-1])

            pl.xlim(.2, .8)
            if col==False: pl.xlim(7000, 5000)
            pl.ylim(0, 70)
            pl.legend(loc='upper left')
            pl.savefig("p_vs_t%s"%i)

if __name__ == "__main__":

    data = np.genfromtxt('/Users/angusr/Python/Gyro/data/recovered.txt').T
#     data = np.genfromtxt('/Users/angusr/Python/Gyro/data/matched_data.txt').T

    # Load data
    plots = plotting(data)
#     pars = [-.6, 0.5189, 0.2] # fitting by eye
#     pars = [0.14510016, 0.59600838, 0.32905815] # results of working
    pars = [np.log10(.7725), 0.5189, 0.601] # Barnes

    plots.p_vs_t(pars, plots.log_period_model, col=True)
#     plots.p_vs_t(pars, plots.log_period_model, col=False)
