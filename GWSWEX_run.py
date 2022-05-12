#%%
import numpy as np
from scipy.integrate import quad
import GWSWEX
import os
import matplotlib.pyplot as plt

def vanGI(d):
    d = d/100
    def theta(h_c):
        theta_r = vanG_pars[0]
        theta_s = vanG_pars[1]
        alpha = vanG_pars[2]
        n = vanG_pars[3]
        m = (1-(1/n))
        return theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
    return np.float64(quad(theta,d,0)[0])*100

#%% in mm and s
elems = int(100)
nts = int(1000)
dt = int(600)
gok = np.random.default_rng().uniform(-3, 3, elems) + 1000
bot = gok - 800
n = np.full(elems, 0.4)
k = np.full(elems, 333e-5)
vanG_pars = np.array([0.1, 0.4, 0.7, 1.5])
chd = np.full(elems, False, dtype=bool)
p = np.full((elems,nts+1), 515e-5)
p[:,-500:] = 490e-5
et = np.full((elems,nts+1), 500e-5)
gws_ini = bot + 400
sws_ini = np.random.default_rng().uniform(0, 1e-1, elems)
epv_ini = (gok-gws_ini)*n
for x in range(elems):
    sm_ini = vanGI(bot[x]-gws_ini[x])

(gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff) = GWSWEX.run(gws_ini, sws_ini, sm_ini, epv_ini, nts, elems, n, dt, k, bot, chd, p, et, gok, vanG_pars)
gws = np.array(gws, dtype=np.double, order="c")
sws = np.array(sws, dtype=np.double, order="c")
sm = np.array(sm, dtype=np.double, order="c")
epv = np.array(epv, dtype=np.double, order="c")
gw_dis = np.array(gw_dis, dtype=np.double, order="c")
sw_dis = np.array(sw_dis, dtype=np.double, order="c")
sm_dis = np.array(sm_dis, dtype=np.double, order="c")
Qin = np.array(Qin, dtype=np.double, order="c")
Qout = np.array(Qout, dtype=np.double, order="c")
Qdiff = np.array(Qdiff, dtype=np.double, order="c")

#%%
fig_path = "output/figs/"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)
elem = 0
plotWlev = False
plotPrec = False
plotDis = False
plotBal = False
savefig = True
dDPI = 90
pDPI = 900
alpha_scatter = 0.7
scatter_size = 3
format = "jpg" #svg, png, jpg
pal = ["#E74C3C", "#2ECC71", "#5EFCA1", "#E7AC3C", "#2980B9", "#1A3BFF", "#FF6600"] #[gw, sm, epv, sv, sw, p, et]

def disPlot(elem):
    plt.figure(dpi=dDPI)
    plt.xlabel("Time Steps")
    plt.ylabel("Discharges in Storage")
    plt.scatter(range(0,nts), gw_dis[elem,:], label="GW_dis", color=pal[0],\
    alpha=alpha_scatter, s=scatter_size)
    plt.scatter(range(0,nts), sm_dis[elem,:], label="SM_dis", color=pal[1],\
    alpha=alpha_scatter, s=scatter_size)
    plt.scatter(range(0,nts), sw_dis[elem,:], label="SW_dis", color=pal[4],\
    alpha=alpha_scatter, s=scatter_size)
    plt.legend(loc="best", fontsize="small")

def wlevPlot(elem,gws,sws,sm):
    plt.figure(dpi=dDPI)
    plt.xlabel("Time Steps")
    plt.ylabel("Water Levels")
    plt.ylim([bot[elem]-10, np.nanmax(sws[elem,:])+25+gok[elem]])
    plt.stackplot(range(0,nts-1), gws[elem,1:], sm[elem,1:],\
    epv[elem,1:]-sm[elem,1:], (np.full(nts-1,gok[elem])-gws[elem,1:])*(1-n[elem]),\
    sws[elem,1:], labels=["Groundwater","Soil Moisture", "Effective Pore Volume", "Soil Volume", "Surface Water"], colors=pal)
    plt.plot(range(0,nts+1), np.full(nts+1,gok[elem]), color="brown", linewidth=0.5, label="Ground Level")
    plt.plot(range(0,nts+1), np.full(nts+1,bot[elem]), color="black", linewidth=0.75, label="Bottom")
    plt.legend(loc=1, fontsize=3)

def precPlot():
    plt.figure(dpi=dDPI)
    plt.xlabel("Time Steps")
    plt.ylabel("Precipitation")
    plt.scatter(range(0,nts+1), p, s=0.1, color=pal[6])

def balPlot():
    plt.figure(dpi=dDPI)
    ind = np.random.choice(Qdiff.shape[0], 1, replace=False)[0]
    plt.xlabel("Time Steps")
    plt.ylabel("Mass Balance Error (Total = {:.2g})".format(Qdiff[ind].sum()))
    plt.plot(Qdiff[ind], "r")

if plotDis:
    disPlot(0)
    if savefig:
        plt.savefig(fig_path+"discharges."+format, format=format, dpi=pDPI)

if plotWlev:
    wlevPlot(0,gws,sws,sm)
    if savefig:
        plt.savefig(fig_path+"water_levels."+format, format=format, dpi=pDPI)

if plotPrec:
    precPlot()
    if savefig:
        plt.savefig(fig_path+"prec."+format, format=format, dpi=pDPI)

if plotBal:
    balPlot()
    if savefig:
        plt.savefig(fig_path+"mBal."+format, format=format, dpi=pDPI)
