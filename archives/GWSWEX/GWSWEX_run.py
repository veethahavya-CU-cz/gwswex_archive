#%%
from fortrapper import gwswex
import numpy as np
import os
import matplotlib.pyplot as plt

elems = int(50000)
nts = int(150)
dt = int(6)
gok = np.random.default_rng().uniform(-3, 3, elems) + 100
gwswex.build(elems, nts, dt, gok)

n = 0.3
m = 0.2
n_gw = 0.3
alpha = 0.1
beta = 0.85
sw_th = 1e-9
k = np.full(elems, 1e-2) 
gwswex.params(n, m, n_gw, beta, alpha, sw_th, k)

chd = np.full(elems, False, dtype=bool)
p = np.full(nts, 1e-1)
p[-60:] = 5e-2
et = np.full(nts, 75e-3)
gws_ini = np.random.default_rng().uniform(-1, 1, elems) + 40
sws_ini = np.random.default_rng().uniform(0, 1e-1, elems)
fc_ini = (gok-gws_ini)*n
sm_ini = fc_ini*0.5
gwswex.init(chd, p, et, gws_ini, sws_ini, fc_ini, sm_ini)

gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff = np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F'), np.zeros((elems,nts),dtype=np.float64,order='F')
gwswex.run(gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)

#%%
fig_path = "output/figs/"
if not os.path.exists(fig_path):
    os.mkdir(fig_path)
elem = 0
plotPrec = False
plotDis = False
plotBal = False
savefig = False
dDPI = 90
pDPI = 1600
alpha_scatter = 0.7
scatter_size = 3
format = "svg" #svg, png, jpg
pal = ["#E74C3C", "#2ECC71", "#5EFCA1", "#E7AC3C", "#2980B9", "#1A3BFF", "#FF6600"] #[gw, sm, epv, sv, sw, p, et]

if plotDis:
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
    if savefig:
        plt.savefig(fig_path+"discharges."+format, format=format, dpi=pDPI)

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
if savefig:
    plt.savefig(fig_path+"water_levels."+format, format=format, dpi=pDPI)

if plotPrec:
    plt.figure(dpi=dDPI)
    plt.xlabel("Time Steps")
    plt.ylabel("Precipitation")
    plt.scatter(range(0,nts), p, s=0.1, color=pal[6])
    if savefig:
        plt.savefig(fig_path+"prec."+format, format=format, dpi=pDPI)

if plotBal:
    plt.figure(dpi=dDPI)
    ind = np.random.choice(Qdiff.shape[0], 1, replace=False)[0]
    plt.xlabel("Time Steps")
    plt.ylabel("Mass Balance Error (Total = {:.2g})".format(Qdiff[ind].sum()))
    plt.plot(Qdiff[ind], "r")
    if savefig:
        plt.savefig(fig_path+"mBal."+format, format=format, dpi=pDPI)
