#%%
from fortrapper import gwswex
import numpy as np
from scipy.integrate import quad
import os
import matplotlib.pyplot as plt

def vanGI(d):
    def theta(h_c):
        theta_s = 0.4
        theta_r = 0.1
        alpha = 0.4
        n = 2.5
        return theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**(1-(1/n)))
    return np.float64(quad(theta,d,0)[0])

def kSM(sm):
    return np.float64(1e-2)

elems = int(1)
nts = int(100)
dt = int(600)
gok = np.random.default_rng().uniform(-3, 3, elems) + 100
n = 0.3
gwswex.build(elems, nts, dt, gok, n)

chd = np.full(elems, False, dtype=bool)
p = np.full((elems,nts), 55e-4)
p[:,-50:] = 45e-4
et = np.full((elems,nts), 50e-4)
gwswex.init(chd, p, et)

gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff = np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F')
gws[:,0] = np.random.default_rng().uniform(-1, 1, elems) + 40
sws[:,0] = np.random.default_rng().uniform(0, 1e-1, elems)
epv[:,0] = (gok-gws[:,0])*n
sm[:,0] = epv[:,0]*0.5
gwswex.run(vanGI, kSM, gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)

#%%
fig_path = "output/figs/"
if not os.path.exists(fig_path):
    os.mkdir(fig_path)
elem = 0
plotWlev = False
plotPrec = False
plotDis = False
plotBal = False
savefig = False
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

def wlevPlot(elem):
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
    disPlot()
    if savefig:
        plt.savefig(fig_path+"discharges."+format, format=format, dpi=pDPI)

if plotWlev:
    wlevPlot()
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