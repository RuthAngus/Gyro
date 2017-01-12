import numpy as np
import matplotlib.pyplot as plt
from colours import plot_colours
cols = plot_colours()

plotpar = {'axes.labelsize': 18,
           'font.size': 18,
           'legend.fontsize': 14,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)

# load data
bv, bv_err, p, p_err, a, a_err, logg, logg_err, f = \
    np.genfromtxt("../data/clusters.txt").T
has, hps = a[f == 4], p[f == 4]  # hyades
hbvs = bv[f == 4]
pas, pps = a[f == 5], p[f == 5]  # praesepe
pbvs = bv[f == 5]
nas, nps = a[f == 6], p[f == 6]  # NGC6811
nbvs = bv[f == 6]
cas, cps = a[f == 7], p[f == 7]  # Coma ber
cbvs = bv[f == 7]

# load 16Cyg A
_, _, p16, perr16, a16, aerr16, _, _, _, _ = \
    np.genfromtxt("../data/16CygA.txt").T

# load NGC6819
bv19s, p19s = np.genfromtxt("../data/NGC6819.txt", skip_header=1).T
a19s = np.ones_like(p19s)*2.5

# select based on colour
select = lambda min_bv, max_bv, bvs: (min_bv < bvs) * (bvs < max_bv)
min_bv, max_bv = .6, .7
m = select(min_bv, max_bv, hbvs)
has2, hps2 = has[m], hps[m]
m = select(min_bv, max_bv, pbvs)
pas2, pps2 = pas[m], pps[m]
m = select(min_bv, max_bv, nbvs)
nas2, nps2 = nas[m], nps[m]
m = select(min_bv, max_bv, cbvs)
cas2, cps2 = cas[m], cps[m]
m = select(min_bv, max_bv, bv[-5:])
fas, fps = a[-5:][m], p[-5:][m]
fa_err, fp_err = a_err[-5:][m], p_err[-5:][m]
fas, fps = a[-5:], p[-5:]
fa_err, fp_err = a_err[-5:], p_err[-5:]
m = select(min_bv, max_bv, bv19s)
a192, p192 = a19s[m], p19s[m]

# find medians of clusters
ca, cp = np.median(cas), np.median(cps2)
ha, hp = np.median(has), np.median(hps2)
pa, pp = np.median(pas), np.median(pps2)
na, npp = np.median(nas), np.median(nps2)
a19, p19 = np.median(a192), np.median(p192)


# plot model
def log_period_model(par, log_a, bv):
    return np.log10(par[0]) + par[1] * log_a + par[2] * np.log10(bv - par[3])
pars = [.407, .566, .325, .495]  # MH
pars_err = [0.008, 0.021, 0.024, 0.010]
xs = np.linspace(0, 20, 100)
plt.clf()
plt.plot(xs, 10**log_period_model(pars, np.log10(xs*1000), .65), '.5',
         ls="--")

# plot clusters + field stars
ms = 8
plt.errorbar(a16, p16, xerr=aerr16, yerr=perr16, fmt="o",
             color="k", mec="k", capsize=0, ms=6)
plt.errorbar(fas, fps, xerr=fa_err, yerr=fp_err, fmt="o",
             color="k", mec="k", capsize=0, ms=6)
# plt.errorbar(ha, hp, xerr=([.1], [.1]), yerr=(np.std(hps2)), fmt="o",
#              color=cols.lightblue, mec=cols.lightblue, capsize=0, ms=ms)
# plt.errorbar(pa, pp, xerr=([.1], [.1]), yerr=(np.std(pps2)), fmt="o",
#              color=cols.orange, mec=cols.orange, capsize=0, ms=ms)
# plt.errorbar(na, npp, xerr=([.1], [.1],), yerr=(np.std(nps2)), fmt="o",
#              color=cols.pink, mec=cols.pink, capsize=0, ms=ms)
# plt.errorbar(ca, cp, xerr=([.1], [.1]), yerr=(np.std(cps2)), fmt="o",
#              color=cols.purple, mec=cols.purple, capsize=0, ms=ms)
# plt.errorbar(a19, p19, xerr=([.1], [.1]), yerr=(np.std(p19)), fmt="o",
#              color=cols.green, mec=cols.green, capsize=0, ms=ms)

plt.errorbar(ca, cp, xerr=([.1], [.1]),
             yerr=([cp-min(cps2)], [max(cps2)]-cp), fmt="o",
             color=cols.purple, mec=cols.purple, capsize=0, ms=ms,
             label="$\mathrm{Coma~Berenices}$")
plt.errorbar(pa, pp, xerr=([.1], [.1]),
             yerr=([pp-min(pps2)], [max(pps2)-pp]), fmt="o",
             color=cols.orange, mec=cols.orange, capsize=0, ms=ms,
             label="$\mathrm{Praesepe}$")
plt.errorbar(ha, hp, xerr=([.1], [.1]),
             yerr=([hp-min(hps2)], [max(hps2)-hp]),
             fmt="o", color=cols.lightblue, mec=cols.lightblue, capsize=0,
             ms=ms, label="$\mathrm{Hyades}$")
plt.errorbar(na, npp, xerr=([.1], [.1],),
             yerr=([npp-min(nps2[1:])], [max(nps2)-npp]), fmt="o",
             color=cols.pink, mec=cols.pink, capsize=0, ms=ms,
             label="$\mathrm{NGC~6811}$")
plt.errorbar(a19, p19, xerr=([.1], [.1]),
             yerr=([p19-min(p192)], [max(p192)-p19]), fmt="o",
             color=cols.green, mec=cols.green, capsize=0, ms=ms,
             label="$\mathrm{NGC~6819}$")

plt.loglog()
plt.xlim(np.log(.01), np.log(100000000))
plt.ylim(4, 10e1)
plt.xlabel("$\mathrm{Age~(Gyr)}$")
plt.ylabel("$\mathrm{Rotation~period~(days)}$")

# plot astero
data = np.genfromtxt("../data/garcia_irfm.txt")
subgiant = 4.2
g = data[8] > subgiant
p1 = data[6][g]
p_err1 = data[7][g]
t1 = data[1][g]
t_err1 = data[2][g]
a1 = data[3][g]
a_errp1 = data[4][g]
a_errm1 = data[5][g]
a_err1 = .5 * (a_errp1 + a_errm1)
logg1 = data[8][g]
logg_err1 = .5*(data[9][g]+data[10][g])
feh1 = data[11][g]
feh_err1 = data[12][g]
flag1 = data[13][g]

tmin, tmax = 5300, 6000
m = select(tmin, tmax, t1)
plt.errorbar(fas[-1], fps[-1], xerr=fa_err[-1], yerr=fp_err[-1], fmt="o",
             color=cols.blue, mec=cols.blue, capsize=0, ms=ms,
             label="$\mathrm{The~Sun}$")
plt.errorbar(a1[m], p1[m], xerr=a_errp1[m], yerr=p_err1[m], fmt="o",
             color="k", mec="k", capsize=0, ms=6,
             label="$\mathrm{Asteroseismic~targets}$")
plt.legend(loc="best")
plt.savefig("AAS_talk")

plt.clf()
plt.plot(xs, 10**log_period_model(pars, np.log10(xs*1000), .65), '.5',
         ls="--")

# plot clusters + field stars
ms = 8
# plt.errorbar(a16, p16, xerr=aerr16, yerr=perr16, fmt="o",
#              color="k", mec="k", capsize=0, ms=6)
# plt.errorbar(fas, fps, xerr=fa_err, yerr=fp_err, fmt="o",
#              color="k", mec="k", capsize=0, ms=6)

# plt.errorbar(ha, hp, xerr=([.1], [.1]), yerr=(np.std(hps2)), fmt="o",
#              color=cols.lightblue, mec=cols.lightblue, capsize=0, ms=ms)
# plt.errorbar(pa, pp, xerr=([.1], [.1]), yerr=(np.std(pps2)), fmt="o",
#              color=cols.orange, mec=cols.orange, capsize=0, ms=ms)
# plt.errorbar(na, npp, xerr=([.1], [.1],), yerr=(np.std(nps2)), fmt="o",
#              color=cols.pink, mec=cols.pink, capsize=0, ms=ms)
# plt.errorbar(ca, cp, xerr=([.1], [.1]), yerr=(np.std(cps2)), fmt="o",
#              color=cols.purple, mec=cols.purple, capsize=0, ms=ms)
# plt.errorbar(a19, p19, xerr=([.1], [.1]), yerr=(np.std(p19)), fmt="o",
#              color=cols.green, mec=cols.green, capsize=0, ms=ms)

plt.errorbar(ca, cp, xerr=([.1], [.1]),
             yerr=([cp-min(cps2)], [max(cps2)]-cp), fmt="o",
             color=cols.purple, mec=cols.purple, capsize=0, ms=ms,
             label="$\mathrm{Coma~Berenices}$")
plt.errorbar(pa, pp, xerr=([.1], [.1]),
             yerr=([pp-min(pps2)], [max(pps2)-pp]), fmt="o",
             color=cols.orange, mec=cols.orange, capsize=0, ms=ms,
             label="$\mathrm{Praesepe}$")
plt.errorbar(ha, hp, xerr=([.1], [.1]),
             yerr=([hp-min(hps2)], [max(hps2)-hp]),
             fmt="o", color=cols.lightblue, mec=cols.lightblue, capsize=0,
             ms=ms, label="$\mathrm{Hyades}$")
plt.errorbar(na, npp, xerr=([.1], [.1],),
             yerr=([npp-min(nps2[1:])], [max(nps2)-npp]), fmt="o",
             color=cols.pink, mec=cols.pink, capsize=0, ms=ms,
             label="$\mathrm{NGC~6811}$")
plt.errorbar(a19, p19, xerr=([.1], [.1]),
             yerr=([p19-min(p192)], [max(p192)-p19]), fmt="o",
             color=cols.green, mec=cols.green, capsize=0, ms=ms,
             label="$\mathrm{NGC~6819}$")

plt.loglog()
plt.xlim(np.log(.01), np.log(100000000))
plt.ylim(4, 10e1)
plt.xlabel("$\mathrm{Age~(Gyr)}$")
plt.ylabel("$\mathrm{Rotation~period~(days)}$")

# plot astero
data = np.genfromtxt("../data/garcia_irfm.txt")
subgiant = 4.2
g = data[8] > subgiant
p1 = data[6][g]
p_err1 = data[7][g]
t1 = data[1][g]
t_err1 = data[2][g]
a1 = data[3][g]
a_errp1 = data[4][g]
a_errm1 = data[5][g]
a_err1 = .5 * (a_errp1 + a_errm1)
logg1 = data[8][g]
logg_err1 = .5*(data[9][g]+data[10][g])
feh1 = data[11][g]
feh_err1 = data[12][g]
flag1 = data[13][g]

tmin, tmax = 5300, 6000
m = select(tmin, tmax, t1)

plt.errorbar(fas[-1], fps[-1], xerr=fa_err[-1], yerr=fp_err[-1], fmt="o",
             color=cols.blue, mec=cols.blue, capsize=0, ms=ms,
             label="$\mathrm{The~Sun}$")
plt.errorbar(a1[m], p1[m], xerr=a_errp1[m], yerr=p_err1[m], fmt="o",
             color="k", mec="k", capsize=0, ms=6,
             label="$\mathrm{Asteroseismic~targets}$")

# add M67 (barnes
epic, bv, v, p, perr = np.genfromtxt("barnes_M67.txt", skip_header=1).T
m = (.6 < bv) * (bv < .7)
mean_p = np.mean(p[m])
rms = np.mean(p[m]**2)**.5
print(rms)
# plt.errorbar(4, mean_p, xerr=1, yerr=rms, fmt="o",
#              color="Crimson", mec="Crimson", capsize=0, ms=ms,
#              label="$\mathrm{M67~(Barnes~2016)}$")
print(len(p[m]))

# add M67 (rebecca)
epic, p, bv, bv_0, _ = np.genfromtxt("acf_ls_results.txt", skip_header=1,
                                     delimiter=",").T
mean_p = np.mean(p[m])
print(len(p[m]))
rms = np.median(((p[m]-np.median(p[m]))**2)**.5)
print("rms = ", rms)
plt.errorbar(4, mean_p, xerr=1, yerr=rms, fmt="o",
             color="tomato", mec="tomato", capsize=0, ms=ms,
             label="$\mathrm{M67~(Esselstein,~in~prep)}$")
# plt.plot(np.ones_like(p)*4, p, "o", color="Orchid", mec="Orchid",
#          ms=ms, label="$\mathrm{M67~(Esselstein,~in~prep)}$")

plt.legend(loc="best")
plt.savefig("AAS_talk_M67_astero.pdf")

plt.clf()
plt.hist(p)
plt.savefig("Test")
