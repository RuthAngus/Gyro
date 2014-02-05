# This version plots different groups of coeval stars


import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as sci

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 10,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
pl.rcParams.update(plotpar)
ocols = ['#FF9933','#66CCCC' , '#FF33CC', '#3399FF', '#CC0066', '#99CC99', '#9933FF', '#CC0000']

from matplotlib import rc
rc("font", size=20, family="serif", serif="Computer Sans")
rc("text", usetex=True)

# Load data
data = np.genfromtxt('/Users/angusr/Python/Gyro/data/matched_data.txt').T

class plotting(object):

    def __init__(self, data):
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

    def iso_calc(self, ag, iso_method, age):
        if iso_method == 'mh':
            # Mamajek and Hillenbrand
            a = 0.407 #pm 0.021
            b = 0.325 #pm 0.024
            c = 0.495 #pm 0.010
            n = 0.566 #pm 0.008

        elif iso_method == 'b':
            # Barnes
            a = 0.7725
            b = 0.601
            c = 0.4
            n = 0.5189

        # calculate isochrones
        tt = np.array([4873, 5331, 5670, 6140, 6560, 6840, 6980, 7260])#, 7650]) # from http://www.isthe.com/chongo
        # /tech/astro/HR-temp-mass-table-bymass.html
        tt = np.array([4000, 4250, 4500, 5000, 5500, 6000, 6500])
        bv = np.array([1.38,1.25,1.05,0.87,0.79,0.69,0.59,0.50])#,0.4])
        bv = np.array([1.393, 1.306, 1.237, 1.063, 0.815, 0.656, 0.498])
        # http://www.astro.yale.edu/demarque/green.tbl
        pp = a * (bv - c)**b
        t_ref = 10.0**(np.log10(age[ag]).mean())
        g = (t_ref*1e3)**n
        spl = sci.UnivariateSpline(tt, pp*g)
        xs = np.linspace(min(tt), max(tt), 1000)
        ys = spl(xs)
        return xs, ys, t_ref

    # remove kraft break and subgiant stars
    def kraft_sub(self):
        k = (self.mass < 1.4)*(self.logg>4.0)
        return self.period[k], self.p_err[k], self.age[k], self.teff[k], self.t_err[k]

    # remove subgiants only
    def sub(self):
        k = self.logg>4.0
        return self.period[k], self.p_err[k], self.age[k], self.teff[k], self.t_err[k]

    def p_vs_t(self):

        # remove k break and subgiants
        period, p_err, age, teff, t_err = plotting.kraft_sub(self)

        # Remove zeros
        zero = teff > 0
        # age masks
        a_lim = [0, 1, 2, 3, 4, 5, max(age)]

        # plot
        for i in range(len(a_lim)-1):
            a = (a_lim[i] < age)*(age < a_lim[i+1])*zero

            # Plot data
            pl.clf()
            pl.errorbar(teff[a], period[a], xerr=t_err[a], \
                    yerr=p_err[a], color='k', fmt='o', mec='k', \
                    capsize=0, markersize=5, ecolor='0.8',\
                    label='$%s < \mathrm{Age} <%s \mathrm{Gyr}$'%(a_lim[i], a_lim[i+1]))

            # Add Isochrones
            xs, ys, t_ref = plotting.iso_calc(self, a, 'b', age)
            pl.plot(xs, ys, color = ocols[i], linestyle = '--', linewidth = 2, \
                    label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ $\mathrm{(Barnes~2007)}$' % t_ref)
            xs, ys, t_ref = plotting.iso_calc(self, a, 'mh', age)
            pl.plot(xs, ys, color = ocols[i], linestyle = '-.', linewidth = 2, \
                    label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ $\mathrm{(Mamajaek~\&~Hillenbrand~2008)}$' % t_ref)

            pl.xlabel("$\mathrm{T_{eff}}$")
            pl.ylabel("$\mathrm{P_{rot}}$")
    #             pl.xlim(pl.gca().get_xlim()[::-1])
            pl.xlim(7000, 5000)
            pl.ylim(0, 70)
            pl.legend(loc='upper left')
            pl.savefig("p_vs_t%s"%i)

        pl.clf()
        # Make subplots
        for i in range(len(a_lim)-1):
            a = (a_lim[i] < age)*(age < a_lim[i+1])*zero

            pl.subplot(len(a_lim)-1, 1, i+1)
            # Plot data
            pl.errorbar(teff[a], period[a], xerr=t_err[a], \
                    yerr=p_err[a], color='k', fmt='o', mec='k', \
                    capsize=0, markersize=5, ecolor='0.8',\
                    label='$%s < \mathrm{Age} <%s \mathrm{Gyr}$'%(a_lim[i], a_lim[i+1]))

            # Add Isochrones
            xs, ys, t_ref = plotting.iso_calc(self, a, 'b', age)
            pl.plot(xs, ys, color = ocols[i], linestyle = '--', linewidth = 2, \
                    label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ $\mathrm{(Barnes~2007)}$' % t_ref)
            xs, ys, t_ref = plotting.iso_calc(self, a, 'mh', age)
            pl.plot(xs, ys, color = ocols[i], linestyle = '-.', linewidth = 2, \
                    label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ $\mathrm{(Mamajaek~\&~Hillenbrand~2008)}$' % t_ref)

#             pl.legend(loc='upper left')
            pl.xlim(7000, 5000)
            pl.ylim(0, 70)
            pl.subplots_adjust(hspace=None)

        pl.xlabel("$\mathrm{T_{eff} (K)}$")
        pl.ylabel("$\mathrm{P_{rot} (days)}$")
        pl.savefig("p_vs_t_sub")

        acols = ['#8856a7','#9ebcda']
        # Make figure with all ages (these will still have subgiants in!)
        pl.clf()
        pl.errorbar(self.teff, self.period, xerr=self.t_err, \
                yerr=self.p_err, color=ocols[2], fmt='o', mec=ocols[2], \
                capsize=0, markersize=5, ecolor='0.8')
        a = self.age > 0
        xs, ys, t_ref = plotting.iso_calc(self, a, 'mh', self.age)
        pl.plot(xs, ys, color = '0.5', linestyle = '-.', linewidth = 2, \
                label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ $\mathrm{(Barnes~2007)}$' % t_ref)
        pl.xlim(7000, 5000)
        pl.ylim(0, 70)
        pl.xlabel("$\mathrm{T_{eff} (K)}$")
        pl.ylabel("$\mathrm{P_{rot} (days)}$")
        pl.legend(loc='upper left')
        pl.savefig("p_vs_t_all")


        # Make figure with all ages and a different colour across the kraft break
        pl.clf()
        b = self.teff > 6500
        pl.errorbar(self.teff, self.period, xerr=self.t_err, \
                yerr=self.p_err, color=ocols[2], fmt='o', mec=ocols[2], \
                capsize=0, markersize=5, ecolor='0.8')
        pl.errorbar(self.teff[b], self.period[b], xerr=self.t_err[b], \
                yerr=self.p_err[b], color=ocols[1], fmt='o', mec=ocols[1], \
                capsize=0, markersize=5, ecolor='0.8')
        a = self.age > 0
        xs, ys, t_ref = plotting.iso_calc(self, a, 'mh', self.age)
        pl.plot(xs, ys, color = '0.5', linestyle = '-.', linewidth = 2, \
                label = '$\mathrm{Age} = %.1f$\,$\mathrm{Gyr}$ $\mathrm{(Barnes~2007)}$' % t_ref)
        pl.xlim(7000, 5000)
        pl.ylim(0, 70)
        pl.xlabel("$\mathrm{T_{eff} (K)}$")
        pl.ylabel("$\mathrm{P_{rot} (days)}$")
        pl.legend(loc='upper left')
        pl.savefig("p_vs_t_bw")

        # remaking period vs mass as p_vs_t
        # remove subgiants
        period, p_err, age, teff, t_err = plotting.sub(self)
#         period, p_err, age, teff, t_err = self.period, self.p_err, self.age, self.teff, self.t_err

        pl.clf()
        b = age > 3.
        c = (0 < age)*(age < 3.)

        pl.errorbar(teff[b], period[b], xerr=t_err[b], \
                yerr=p_err[b], color=ocols[2], fmt='o', mec=ocols[2], \
                capsize=0, markersize=5, ecolor='0.8', label = '$>3\mathrm{Gyr}$')
        pl.errorbar(teff[c], period[c], xerr=t_err[c], \
                yerr=p_err[c], color=ocols[1], fmt='o', mec=ocols[1], \
                capsize=0, markersize=5, ecolor='0.8', label = '$<3\mathrm{Gyr}$')

        xs, ys, t_ref = plotting.iso_calc(self, b, 'mh', self.age)
        pl.plot(xs, ys, color=ocols[2], linestyle = '-.', linewidth = 2, \
                label = '$%.1f$\,$\mathrm{Gyr}$ $\mathrm{(M\&H~2008)}$' % t_ref)
        xs, ys, t_ref = plotting.iso_calc(self, c, 'mh', self.age)
        pl.plot(xs, ys, color=ocols[1], linestyle = '-.', linewidth = 2, \
                label = '$%.1f$\,$\mathrm{Gyr}$ $\mathrm{(M\&H~2008)}$' % t_ref)

        pl.xlim(7000, 5000)
        pl.ylim(0, 70)
        pl.xlabel("$\mathrm{T_{eff} (K)}$")
        pl.ylabel("$\mathrm{P_{rot} (days)}$")
        pl.legend(loc='upper left')
        pl.savefig("p_vs_t_orig2")


plots = plotting(data)
plots.p_vs_t()
