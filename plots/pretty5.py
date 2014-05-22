# This version plots different groups of coeval stars

import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as sci

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

class plotting(object):

    def __init__(self):
        data = np.genfromtxt('/Users/angusr/Python/Gyro/data/matched_data.txt').T
        self.KID = data[0]
        self.period = data[1]
        self.p_err = data[2]
        self.age = data[3]
        self.a_errp = data[4]
        self.a_errm = data[5]
        self.mass = data[6]
        self.m_errp = data[7]
        self.m_errm = data[8]
        self.logg = data[9]
        self.l_errp = data[10]
        self.l_errm = data[11]
        self.teff = data[12]
        self.t_err = data[13]
        self.j_k = data[14]

    def log_period_model(self, par, log_age, temp):
        return par[0] + par[1] * log_age + par[2] * np.log10(par[3] - temp)

    def log_period_colour_model(self, par, log_age, bv):
        return par[0] + par[1] * log_age + par[2] * np.log10(bv - par[3])

    def iso_calc(self, ag, pars, age):

        # calculate isochrones
        tt = np.array([4000, 4250, 4500, 5000, 5500, 6000, 6500])
        bv = np.array([1.393, 1.306, 1.237, 1.063, 0.815, 0.656, 0.498])
        # http://www.astro.yale.edu/demarque/green.tbl

        if type(ag) != np.ndarray: t_ref = ag
#         else: t_ref = 10.0**(np.log10(age[ag]).mean())
        else: t_ref = 10.0**(np.median(np.log10(age[ag])))
#         model = a*(bv - c)**b*(t_ref*1e3)**n
        teff = np.linspace(4000, pars[-1]-10, 100)
        model = 10**self.log_period_model(pars, np.log10(t_ref*1000), teff)
#         spl = sci.UnivariateSpline(tt, model)
#         xs = np.linspace(min(tt), max(tt), 1000)
#         ys = spl(xs)
        xs = teff
        ys = model
        return xs, ys, t_ref

    # remove kraft break and subgiant stars
    def kraft_sub(self):
        k = (self.mass < 1.4)*(self.logg>4.0)
        return self.period[k], self.p_err[k], self.age[k], self.teff[k], self.t_err[k]

    # remove subgiants only
    def sub(self):
        k = self.logg>4.0
        return self.period[k], self.p_err[k], self.age[k], self.teff[k], self.t_err[k]

    def p_vs_t(self, pars):

        # remove k break and subgiants
        period, p_err, age, teff, t_err = plotting.kraft_sub(self)

        # Remove zeros
        zero = teff > 0
        # age masks
        a_lim = [0, 1, 2, 3, 4, 5, 8, max(age)]

        # plot
        for i in range(len(a_lim)-1):
            a = (a_lim[i] < age)*(age < a_lim[i+1])*zero

            # Plot data
            pl.clf()
            pl.errorbar(teff[a], period[a], xerr=t_err[a], \
                    yerr=p_err[a], color='k', fmt='o', mec='k', \
                    capsize=0, markersize=5, ecolor='0.8', zorder = 4, \
                    label='$%s < \mathrm{Age} <%s \mathrm{Gyr}$'%(a_lim[i], a_lim[i+1]))

            # Add Isochrones
            xs, ys, t_ref = plotting.iso_calc(self, a, pars, age)
            pl.plot(xs, ys, color = ocols[i], linestyle='-', linewidth = 2, \
                    label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ \
                    $\mathrm{(M\&H~2008)}$' % t_ref, zorder = 1)
            a = a_lim[i]
            xs, ys, t_ref = plotting.iso_calc(self, a, pars, age)
            pl.plot(xs, ys, color = ocols[i], linestyle = '--', linewidth = 2, zorder = 2)
            a = a_lim[i+1]
            xs, ys, t_ref = plotting.iso_calc(self, a, pars, age)
            pl.plot(xs, ys, color = ocols[i], linestyle = '--', linewidth = 2, zorder = 3)

            pl.xlabel("$\mathrm{T_{eff (K)}}$")
            pl.ylabel("$\mathrm{P_{rot} (days)}$")
    #             pl.xlim(pl.gca().get_xlim()[::-1])
            pl.xlim(7000, 5000)
            pl.ylim(0, 70)
            pl.legend(loc='upper left')
            pl.savefig("p_vs_t%s"%i)

if __name__ == "__main__":
    # Load data
    plots = plotting()
#     pars  = [0.407, 0.325, 0.495, 0.566]
#     pars = [np.log10(0.7725), 0.5189, .2, 6300.]
#     pars = [1.51613663, .185575537, -.245929036, 9.04129937e+03]
    pars  = [0.407, 0.325, 0.495, 0.566]
    plots.p_vs_t(pars)
