#%%
from fortrapper import gwswex
import numpy as np
from scipy.integrate import quad
import os
import matplotlib.pyplot as plt
#%%
vanG_pars = np.array([0.1, 0.4, 0.42, 2.5], dtype=np.float64,order='F')

def vanGI(d):
    def theta(h_c):
        theta_r = vanG_pars[0]
        theta_s = vanG_pars[1]
        alpha = vanG_pars[2]
        n = vanG_pars[3]
        m = (1-(1/n))
        return theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
    return np.float64(quad(theta,d,0)[0])

def kSM(s):
    #vanG-Mualem (YT)
    ks = 3e-2
    theta_r = vanG_pars[0]
    theta_s = vanG_pars[1]
    n = vanG_pars[3]
    m = (1-(1/n))
    sat = ((s-theta_r)/(theta_s-theta_r))
    k = ks*sat*((1-(1-(sat)**(1/m))**m)**2)
    return np.float64(k)

#%%
elems = int(1)
nts = int(100+1)
dt = int(600)
gok = np.random.default_rng().uniform(-3, 3, elems) + 100
bot = gok - 20
n = np.full(elems, 0.4)
k = np.full(elems, 60e-3)
gwswex.build(elems, nts, dt, gok, bot, n, k, vanG_pars)

chd = np.full(elems, False, dtype=bool)
p = np.full((elems,nts), 550e-5)
p[:,-50:] = 45e-4
et = np.full((elems,nts), 500e-5)
gwswex.init(chd, p, et)

gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff = np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F')
gws[:,0] = np.random.default_rng().uniform(-1, 1, elems) + 40
sws[:,0] = np.random.default_rng().uniform(0, 1e-1, elems)
epv[:,0] = (gok-gws[:,0])*n
sm[:,0] = epv[:,0]*0.5
gwswex.run(vanGI, gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)

#%%
fig_path = "output/figs/"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)
elem = 0
plotWlev = True
plotPrec = True
plotDis = True
plotBal = True
savefig = True
dDPI = 90
pDPI = 1600
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
    gws = gws[elem,1:]
    plt.ylim([gws.min()-50, sws[elem,:].max()+50+gok[elem]])
    plt.stackplot(range(0,nts-1), gws, sm[elem,1:],\
    epv[elem,1:]-sm[elem,1:], (np.full(nts-1,gok[elem])-gws)*(1-n),\
    sws[elem,1:], labels=["Groundwater","Soil Moisture", "Effective Pore Volume", "Soil Volume", "Surface Water"], colors=pal)
    plt.plot(range(0,nts), np.full(nts,gok[elem]), "k", linewidth=0.5, label="Ground Level")
    plt.legend(loc="best", fontsize="small")

def precPlot():
    plt.figure(dpi=dDPI)
    plt.xlabel("Time Steps")
    plt.ylabel("Precipitation")
    plt.scatter(range(0,nts), p, s=0.1, color=pal[6])

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