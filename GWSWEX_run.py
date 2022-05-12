#%%
import numpy as np
from scipy.io import FortranFile
from scipy.integrate import quad
import os, shutil
import subprocess as sp
import matplotlib.pyplot as plt

vanG_pars = np.array([0.1, 0.4, 0.7, 1.5], dtype=np.float64, order='F')

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

def fwrite(fname, val):
    ip_path = 'exe/fort/input/'
    Ffile = FortranFile(os.path.join(ip_path,fname), 'w')
    Ffile.write_record(val)
    Ffile.close()

def fread(fname):
    shape = (elems, nts)
    op_path = 'exe/fort/output/'
    Ffile = FortranFile(os.path.join(op_path,fname), 'r')
    val = Ffile.read_reals().reshape(shape, order='F')
    Ffile.close()
    return val

#%% in mm and s
if not os.path.exists('exe/fort/output/'):
    os.mkdir('exe/fort/output')
if os.path.exists('exe/fort/input/'):
    shutil.rmtree('exe/fort/input/')
os.mkdir('exe/fort/input/')

elems = int(1)
nts = int(1000)
dt = int(600)
np.savetxt('exe/fort/input/build.dat', np.array([elems, nts, dt], dtype=np.int32), fmt='%d')

gok = np.random.default_rng().uniform(-3, 3, elems)+1000
fwrite('gok.ip', np.array(gok, dtype=np.float64, order='F'))
bot = gok - 800
fwrite('bot.ip', np.array(bot, dtype=np.float64, order='F'))
n = np.full(elems, 0.4)
fwrite('n.ip', np.array(n, dtype=np.float64, order='F'))
k = np.full(elems, 333e-5)
fwrite('k.ip', np.array(k, dtype=np.float64, order='F'))
chd = np.full(elems, False, dtype=bool)
fwrite('chd.ip', np.array(chd, dtype=np.float64, order='F'))
p = np.full((elems,nts+1), 515e-5)
p[:,-500:] = 490e-5
fwrite('p.ip', np.array(p, dtype=np.float64, order='F'))
et = np.full((elems,nts+1), 500e-5)
fwrite('et.ip', np.array(et, dtype=np.float64, order='F'))
fwrite('vanG_pars.ip', np.array(vanG_pars, dtype=np.float64, order='F'))

sm_ini = []
gws_ini = bot + 400
sws_ini = np.random.default_rng().uniform(0, 1e-1, elems)
epv_ini = (gok-gws_ini)*n
for x in range(elems):
    sm_ini.append(vanGI(bot[x]-gws_ini[x]))
fwrite('gws_ini.ip', np.array(gws_ini, dtype=np.float64, order='F'))
fwrite('sws_ini.ip', np.array(sws_ini, dtype=np.float64, order='F'))
fwrite('epv_ini.ip', np.array(epv_ini, dtype=np.float64, order='F'))
fwrite('sm_ini.ip', np.array(sm_ini, dtype=np.float64, order='F'))

wd = os.getcwd()
os.chdir('exe/fort/')
fort_run = sp.Popen('./GWSWEX', shell=True, stdout = sp.PIPE)
fort_run.communicate()
os.chdir(wd)

gws = fread('gws.op')
sws = fread('sws.op')
sm = fread('sm.op')
epv = fread('epv.op')
gw_dis = fread('gw_dis.op')
sw_dis = fread('sw_dis.op')
sm_dis = fread('sm_dis.op')
Qin = fread('Qin.op')
Qout = fread('Qout.op')
Qdiff = fread('Qdiff.op')


#%%
fig_path = "output/figs/"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)
elem = 0
plotWlev = True
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
    plt.ylim([bot[elem]-10, sws[elem,:].max()+25+gok[elem]])
    plt.stackplot(range(0,nts-1), gws[elem,1:], sm[elem,1:],\
    epv[elem,1:]-sm[elem,1:], (np.full(nts-1,gok[elem])-gws[elem,1:])*(1-n[elem]),\
    sws[elem,1:], labels=["Groundwater","Soil Moisture", "Effective Pore Volume", "Soil Volume", "Surface Water"], colors=pal)
    plt.plot(range(0,nts+1), np.full(nts+1,gok[elem]), color="brown", linewidth=0.5, label="Ground Level")
    plt.plot(range(0,nts+1), np.full(nts+1,bot[elem]), color="black", linewidth=0.5, label="Bottom")
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