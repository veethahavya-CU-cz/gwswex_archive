#%%
import numpy as np
from scipy.io import FortranFile
from scipy.integrate import quad
import os, shutil, psutil
import subprocess as sp
import matplotlib.pyplot as plt

os.environ['OMP_NUM_THREADS'] = str(psutil.cpu_count(logical = False))

vanG_pars = np.array([0.02, 0.42, 0.35, 1.25], dtype=np.float64, order='F')

def vanGI(d):
	def theta(h_c):
		theta_r = vanG_pars[0]
		theta_s = vanG_pars[1]
		alpha = vanG_pars[2]
		n = vanG_pars[3]
		m = (1-(1/n))
		return theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
	return np.float64(quad(theta,-d,0)[0])

def fwrite(fname, val):
	ip_path = 'exe/fort/input/'
	Ffile = FortranFile(os.path.join(ip_path,fname), 'w')
	Ffile.write_record(val.T)
	Ffile.close()

def fread(fname):
	shape = (elems, nts)
	op_path = 'exe/fort/output/'
	Ffile = FortranFile(os.path.join(op_path,fname), 'r')
	val = Ffile.read_reals().reshape(shape, order='F')
	Ffile.close()
	return val

def plot(elem, plotWlev=True, plotPrec=True, plotDis=True, plotBal=True, savefig=True, dDPI=90, pDPI=1600, alpha_scatter=0.7, scatter_size=3, format='jpg'):
	#formats = jpg, svg, png, jpg
	fig_path = os.path.join(op_path, 'figs/')
	if not os.path.exists(fig_path):
		os.makedirs(fig_path)
	pal = ['#E74C3C', '#2ECC71', '#5EFCA1', '#E7AC3C', '#2980B9', '#1A3BFF', '#FF6600'] #[gw, sm, epv, sv, sw, p, et]

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
		plt.legend(loc='best', fontsize='small')
		plt.tight_layout()
		plt.xticks(range(0,nts,24*30))

	def wlevPlot(elem,gws,sws,sm):
		plt.figure(dpi=dDPI)
		plt.xlabel("Time Steps (h)")
		plt.ylabel("Water Levels (m.a.s.l.)")
		plt.ylim([bot[elem]-0.5, sws[elem,:].max()+0.5+gok[elem]])
		plt.stackplot(range(0,nts-1), gws[elem,1:], sm[elem,1:],\
		epv[elem,1:]-sm[elem,1:], (np.full(nts-1,gok[elem])-gws[elem,1:])*(1-n[elem]),\
		sws[elem,1:], labels=["Groundwater","Soil Moisture", "Effective Pore Volume", "Soil Volume", "Surface Water"], colors=pal)
		if plotPrec:
			p_dom, et_dom = [], []
			ht = (sws[elem,:].max()+0.5+gok[elem]) + (bot[elem]-0.5)
			for ts in range(nts-1):
				if p[elem,ts] > et[elem,ts]:
					p_dom.append(ht)
					et_dom.append(0)
				else:
					et_dom.append(ht)
					p_dom.append(0)
			plt.stackplot(range(0,nts-1), p_dom, labels=["Precipitation Dominant", ], colors=['#A8EAED'], alpha=0.21)
			plt.stackplot(range(0,nts-1), et_dom, labels=["Evapotranspiration Dominant", ], colors=['#E8A78B'], alpha=0.21)
		plt.plot(range(0,nts+1), np.full(nts+1,gok[elem]), color='#502D16', linewidth=0.5, label="Ground Level")
		plt.plot(range(0,nts+1), np.full(nts+1,bot[elem]), color='black', linewidth=0.5, label="Bottom")
		plt.legend(loc='lower right', fontsize=3)
		plt.tight_layout()
		plt.xticks(range(0,nts,24*30))

	def balPlot():
		plt.figure(dpi=dDPI)
		ind = np.random.choice(Qdiff.shape[0], 1, replace=False)[0]
		plt.xlabel("Time Steps")
		plt.ylabel("Mass Balance Error (Total = {:.2g})".format(Qdiff[ind].sum()))
		plt.plot(Qdiff[ind], "r")
		plt.tight_layout()
		plt.xticks(range(0,nts,24*30))

	if plotDis:
		disPlot(0)
		if savefig:
			plt.savefig(os.path.join(fig_path,"discharges."+format), format=format, dpi=pDPI)

	if plotWlev:
		wlevPlot(0,gws,sws,sm)
		if savefig:
			plt.savefig(os.path.join(fig_path,"water_levels."+format), format=format, dpi=pDPI)

	if plotBal:
		balPlot()
		if savefig:
			plt.savefig(os.path.join(fig_path,"mBal."+format), format=format, dpi=pDPI)

#%% all units in SI
if not os.path.exists('exe/fort/output/'):
	os.mkdir('exe/fort/output')
if os.path.exists('exe/fort/input/'):
	shutil.rmtree('exe/fort/input/')
os.mkdir('exe/fort/input/')

logger_switch = 1 #0: no logger, 1: logger
elems = int(1)
nts = int(24*30*6) #a total of 4320 time steps representing ~6 months
dt = int(60*60*1) #at 1h intervals
np.savetxt('exe/fort/input/build.dat', np.array([elems, nts, dt, logger_switch], dtype=np.int32), fmt='%d')

gok = np.random.default_rng().uniform(-3, 3, elems)+150
fwrite('gok.ip', np.array(gok, dtype=np.float64, order='F'))
bot = gok - 30
fwrite('bot.ip', np.array(bot, dtype=np.float64, order='F'))
n = np.full(elems, vanG_pars[1])
fwrite('n.ip', np.array(n, dtype=np.float64, order='F'))
k = np.full(elems, 50e-5)
fwrite('k.ip', np.array(k, dtype=np.float64, order='F'))
chd = np.full(elems, False, dtype=bool)
fwrite('chd.ip', np.array(chd, dtype=np.float64, order='F'))
gw_sm_interconnectivity = np.full(elems, vanG_pars[0], dtype=np.float64)
fwrite('gw_sm_interconnectivity.ip', np.array(gw_sm_interconnectivity, dtype=np.float64, order='F'))
macropore_inf_degree = np.full(elems, 0, dtype=np.float64)
fwrite('macropore_inf_degree.ip', np.array(macropore_inf_degree, dtype=np.float64, order='F'))
p = np.full((elems,nts+1), 2.5*(1e-3/3600))
#p[:,int(-nts/2):] = 0*(1e-3/3600)
p[:,0:500] = 3.5*(1e-3/3600)
p[:,500:750] = 0*(1e-3/3600)
p[:,1000:1250] = 0*(1e-3/3600)
p[:,1000:1250] = 0*(1e-3/3600)
p[:,1750:2000] = 0*(1e-3/3600)
p[:,2100:nts] = 0*(1e-3/3600)
fwrite('p.ip', np.array(p, dtype=np.float64, order='F'))
et = np.full((elems,nts+1), 0.33*(1e-3/3600))
fwrite('et.ip', np.array(et, dtype=np.float64, order='F'))
fwrite('vanG_pars.ip', np.array(vanG_pars, dtype=np.float64, order='F'))

sm_ini = []
gws_ini = bot + 5
sws_ini = np.random.default_rng().uniform(0, 1e-2, elems)
epv_ini = (gok-gws_ini)*n
for x in range(elems):
	sm_ini.append(vanGI(gok[x]-gws_ini[x]))
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

op_path = 'output/'
if not os.path.exists(op_path):
	os.makedirs(op_path)
np.savez(os.path.join(op_path, 'GWSWEX.npz'), gws=gws, sws=sws, sm=sm, epv=epv, gw_dis=gw_dis, sw_dis=sw_dis, sm_dis=sm_dis, Qin=Qin, Qout=Qout, Qdiff=Qdiff)

plot(0, plotWlev=True, plotPrec=True, plotDis=True, plotBal=True, savefig=True) #True False

#%%
influx = (p[0].sum()-p[0,0])*dt - (et[0].sum()-et[0,0])*dt
delta_storages = sm[0,-1]-sm[0,0] + (gws[0,-1]-gws[0,0])*vanG_pars[1] + sws[0,-1]-sws[0,0]
print("mbal err: {:.2e}".format(influx-delta_storages))
print("mbal % err: {:.2F}".format(((influx-delta_storages)*100/influx)))
