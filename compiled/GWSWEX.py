import time
import numpy as np
import pandas as pd
import sys, os
import subprocess as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.io import FortranFile
from datetime import datetime, timedelta
import netCDF4 as nc
import flopy as fp
import fileinput
import warnings
import logging
import shutil

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
pd.options.mode.chained_assignment = None

if os.path.exists("../output/GWSWEX.log"):
	with open("../output/GWSWEX.log", "w"):
		pass
else:
    os.mkdir("../output/")
    with open("../output/GWSWEX.log", "w"):
        pass
logging.basicConfig(filename="../output/GWSWEX.log", format="%(asctime)s; %(levelname)s: %(message)s", level=logging.DEBUG)

#%%
class bundleIT(object):
	def __init__(self, **argd):
		self.__dict__.update(argd)
	def keys(self):
		for property, value in vars(self).items():
			print(property)
	def __main__(self):
		self.keys()


#%%
class timing():
	def __init__(self, start, end, dt_exchange, ts_exchange):
		self.str = {}
		self.unix = {}
		self.rel = {}
		self.stmp = {}
		self.dt = {}
		self.nTS = {}
		self.ts = {}
		self.run = {}

		self.run["PET.prep"] = []
		self.run["PET.get"] = []
		self.run["Fort.build"] = []
		self.run["Fort.Run"] = []
		self.run["Fort.update"] = []
		self.run["Fort.load"] = []
		self.run["Delft.readNC"] = []
		self.run["Delft.Run"] = []
		self.run["Delft.update"] = []
		self.run["Delft.load"] = []
		self.run["MF.build"] = []
		self.run["MF.Run"] = []
		self.run["MF.update"] = []
		self.run["MF.load"] = []
		self.run["start"] = time.time()
		self.runtimes =[]

		self.str["global_start"] = start
		self.unix["global_start"] = time.mktime(datetime.strptime(self.str["global_start"], "%Y-%m-%d %H:%M:%S").timetuple())
		self.rel["global_start"] = 0.0
		self.stmp["global_start"] = datetime.strptime(self.str["global_start"], "%Y-%m-%d %H:%M:%S")
		self.str["global_end"] = end
		self.unix["global_end"] = time.mktime(datetime.strptime(self.str["global_end"], "%Y-%m-%d %H:%M:%S").timetuple())
		self.rel["global_end"] = self.unix["global_end"] - self.unix["global_start"]
		self.stmp["global_end"] = datetime.strptime(self.str["global_end"], "%Y-%m-%d %H:%M:%S")

		self.nTS["exchange"] = ts_exchange
		self.dt["exchange"] = dt_exchange
		self.dt["FortRun"] = self.dt["exchange"]/self.nTS["exchange"]
		self.nTS["ran"] = 0
		self.nTS["run_num"] = 0
		self.nTS["max"] = int(self.rel["global_end"]/self.dt["exchange"] - 1)

		self.str["local_start"] = start
		self.unix["local_start"] = time.mktime(datetime.strptime(self.str["local_start"], "%Y-%m-%d %H:%M:%S").timetuple())
		self.rel["local_start"] = self.unix["local_start"] - self.unix["global_start"]
		self.stmp["local_start"] = datetime.strptime(self.str["local_start"], "%Y-%m-%d %H:%M:%S")
		self.unix["local_end"] = self.unix["local_start"] + self.dt["exchange"]
		self.str["local_end"] = datetime.fromtimestamp(self.unix["local_end"]).strftime("%Y-%m-%d %H:%M:%S")
		self.rel["local_end"] = self.unix["local_end"] - self.unix["global_start"]
		self.stmp["local_end"] = datetime.strptime(self.str["local_end"], "%Y-%m-%d %H:%M:%S")

		self.ts["local"] = 0
		self.ts["global"] = self.ts["local"]

		start = self.stmp["local_start"].hour*3600 +self.stmp["local_start"].minute*60 + self.stmp["local_start"].second
		stop = self.stmp["local_end"].hour*3600 + self.stmp["local_end"].minute*60 + self.stmp["local_end"].second
		self.rel["delft"] = np.array([start, stop])
		self.rel["delft_tim"] =self.rel["delft"]/60

	def update(self, update_by=None, internal=False):
		self.run["start"] = time.time()
		if update_by is None:
			dt_update = self.dt["exchange"]
		else:
			dt_update = update_by
		self.unix["local_start"] = self.unix["local_start"] + dt_update
		self.str["local_start"] = datetime.fromtimestamp(self.unix["local_start"]).strftime("%Y-%m-%d %H:%M:%S")
		self.rel["local_start"] = self.unix["local_start"] - self.unix["global_start"]
		self.stmp["local_start"] = datetime.strptime(self.str["local_start"], "%Y-%m-%d %H:%M:%S")
		self.unix["local_end"] = self.unix["local_end"] + dt_update
		self.str["local_end"] = datetime.fromtimestamp(self.unix["local_end"]).strftime("%Y-%m-%d %H:%M:%S")
		self.rel["local_end"] = self.unix["local_end"] - self.unix["global_start"]
		self.stmp["local_end"] = datetime.strptime(self.str["local_end"], "%Y-%m-%d %H:%M:%S")
		self.ts["local"] = np.arange(self.rel["local_start"], self.rel["local_end"]+1, self.dt["FortRun"])
		self.ts["global"] = np.append(self.ts["global"], self.ts["local"][1:])
		self.rel["delft"] = self.rel["delft"] + dt_update
		self.rel["delft_tim"] = self.rel["delft_tim"] + dt_update/60
		if internal:
			logging.info("TIMING: Updated times to suit nTS change")
		else:
			logging.info("*TIMING: Updated. Start: {0}   End: {1}".format(self.str["local_start"], self.str["local_end"]))

	def getRuntime(self, n):
		max_n = max(len(self.run["Fort.Run"]), len(self.run["Delft.Run"]), len(self.run["MF.Run"]))
		if n > max_n:
			raise ValueError("Run number is out of bounds")
		runtime = []
		for v in self.run.values():
			if v:
				runtime.append(v[0])
		runtime = sum(runtime)
		self.runtimes.append(runtime)
		logging.info("*RUNTIME{}: {}".format(self.nTS["run_num"], self.runtimes[int(self.nTS["run_num"])]))


#%%
class PET:
	def __init__(self, data_path, times, Fort, p_name="p.dat", et_name="et.dat"):
		self.data_path = data_path
		self.p_path = os.path.join(data_path, p_name)
		self.et_path = os.path.join(data_path, et_name)
		self.times = times
		self.Fort = Fort
		self.pDF = None
		self.etDF = None
		self.p = None
		self.et = None
		self.petDF_sel = None

	def prep(self):
		strt_time = time.time()
		pDF = pd.read_table(self.p_path)
		pDF[pDF.columns[0]] = pd.to_datetime(pDF[pDF.columns[0]])
		pDF.index = pDF[pDF.columns[0]]
		pDF = pDF.drop(pDF.columns[0], 1)
		etDF = pd.read_table(self.et_path)
		etDF[etDF.columns[0]] = pd.to_datetime(etDF[etDF.columns[0]])
		etDF.index = etDF[etDF.columns[0]]
		etDF = etDF.drop(etDF.columns[0], 1)
		start = datetime(self.times.stmp["global_start"].year, self.times.stmp["global_start"].month, self.times.stmp["global_start"].day,\
		self.times.stmp["global_start"].hour)
		end = datetime(self.times.stmp["global_end"].year, self.times.stmp["global_end"].month, self.times.stmp["global_end"].day,\
		self.times.stmp["global_end"].hour)
		p_mask = (pDF.index >= start) & (pDF.index <= end)
		et_mask = (etDF.index >= start) & (etDF.index <= end)
		pDF = pDF.loc[p_mask]
		etDF = etDF.loc[et_mask]
		self.pDF = pDF
		self.etDF = etDF
		self.times.run["PET.prep"].append(time.time() - strt_time)
		logging.info("PET: Prepped")

	def get(self, throttle=False, new_dt_exchange=None):
		strt_time = time.time()
		start = datetime(self.times.stmp["local_start"].year, self.times.stmp["local_start"].month, self.times.stmp["local_start"].day,\
		self.times.stmp["local_start"].hour)
		end = datetime(self.times.stmp["local_end"].year, self.times.stmp["local_end"].month, self.times.stmp["local_end"].day,\
		self.times.stmp["local_end"].hour)
		p_mask = (self.pDF.index >= start) & (self.pDF.index <= end)
		pDF_sel = self.pDF.loc[p_mask]
		p = pDF_sel[pDF_sel.columns[0]].to_numpy()[0]
		if throttle:
			if not new_dt_exchange:
				new_dt_exchange = self.times.dt["exchange"]
			if p.max() == 0:
				self.times.nTS["exchange"] = 1
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.1:
				self.times.nTS["exchange"] = 2
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.2:
				self.times.nTS["exchange"] = 4
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.3:
				self.times.nTS["exchange"] = 6
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 0.5:
				self.times.nTS["exchange"] = 8
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			elif p.max() <= 1:
				self.times.nTS["exchange"] = 10
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			else:
				self.times.nTS["exchange"] = 15
				old_dt_exch = self.times.dt["exchange"]
				self.times.dt["exchange"] = new_dt_exchange
			self.times.dt["FortRun"] = self.times.dt["exchange"]/self.times.nTS["exchange"]
			self.times.update(update_by=self.times.dt["exchange"]-old_dt_exch, internal=True)
			logging.warning("PET: !Throttling! nTS: {0}. ts_exch: {1}s".format(self.times.nTS["exchange"], self.times.dt["exchange"]))
		p = np.repeat(p, (self.times.nTS["exchange"])+1)
		p = p/3600
		logging.info("PET: max P for this ts is {} mm/s".format(p.max()))
		et_mask = (self.etDF.index >= start) & (self.etDF.index <= end)
		etDF_sel = self.etDF.loc[et_mask]
		et = etDF_sel[etDF_sel.columns[0]].to_numpy()[0]
		et = np.repeat(et, (self.times.nTS["exchange"]/et.size)+1)
		et = et/3600
		logging.info("PET: max ET for this ts is {} mm/s".format(et.max()))
		if pDF_sel.index.equals(etDF_sel.index):
			petDF_sel = pDF_sel
			petDF_sel[etDF_sel.columns[0]] = etDF_sel[etDF_sel.columns[0]]
		else:
			logging.errorr("WARNING: P and ET indexes do not match")
			sys.exit("WARNING: P and ET indexes do not match")
		self.Fort.Ini.p = p
		self.Fort.Ini.et = et
		self.Fort.Ini.ts = p.size
		self.Fort.Ini.ts_size = self.times.dt["FortRun"]
		self.p = p
		self.et = et
		self.petDF_sel = petDF_sel
		self.times.run["PET.get"].append(time.time() - strt_time)
		logging.info("PET: Fetched and Passed P and ET values to FORT")


#%%
class Fort:
	def __init__(self, path, times, exe_path="GWSWEX.exe"):
		self.path = path
		self.exe_path = os.path.join(path, exe_path)
		self.times = times
		self.Ini = bundleIT(elems=None, ts=None, ts_size=None, n=None, m=None, beta=None, alpha=None, gok=None, k=None, sw_th=None,\
		p=None, et=None, gws=None, sws=None, epv=None, sm=None, chd=None)
		self.Res = bundleIT(gws=None, sws=None, sm=None, epv=None, Qin=None, Qout=None, Qdiff=None, gw_dis=None, sw_dis=None, sm_dis=None,\
		max_err_abs=None, max_err_perc=None)

	def build(self, run=False, restart=False, res=None):
		strt_time = time.time()
		if restart:
			fort_file = os.path.join(res.fort_path, "wasenmoos_fort_"+(self.times.stmp["local_start"]-timedelta(seconds=\
			self.times.dt["exchange"])).strftime("%Y%m%d_%H%M%S")+".npz")
			fort_rst_file = np.load(fort_file)
			self.Ini.sm = fort_rst_file["sm"][:,-1]
			self.Ini.epv = fort_rst_file["epv"][:,-1]
		self.Ini.chd = np.zeros(self.Ini.elems, dtype=np.int8)
		for elem in range(self.Ini.elems):
			if elem in self.Ini.chd_cells:
				self.Ini.chd[elem] = 1
			else:
				self.Ini.chd[elem] = 0
		missing = []
		for attr in self.Ini.__dict__.keys():
			if getattr(self.Ini, attr) is None:
				missing.append(attr)
		if len(missing) != 0:
			logging.error("FORT: Incomplete set of FortIni Objects. Missing values: ", missing[:])
			raise ValueError("Incomplete set of FortIni Objects. Missing values: ", missing[:])
		p_file = FortranFile(os.path.join(self.path,"p.ip"), "w")
		p_file.write_record(self.Ini.p)
		p_file.close()
		et_file = FortranFile(os.path.join(self.path,"et.ip"), "w")
		et_file.write_record(self.Ini.et)
		et_file.close()
		args = [self.Ini.elems, self.Ini.ts, self.Ini.ts_size, self.Ini.n, self.Ini.n_gw, self.Ini.m, self.Ini.beta,\
		self.Ini.alpha, self.Ini.sw_th]
		np.savetxt(os.path.join(self.path,"args.ip"), args)
		gws_ini_file = FortranFile(os.path.join(self.path,"gws_ini.ip"), "w")
		gws_ini_file.write_record(self.Ini.gws)
		gws_ini_file.close()
		sws_ini_file = FortranFile(os.path.join(self.path,"sws_ini.ip"), "w")
		sws_ini_file.write_record(self.Ini.sws)
		sws_ini_file.close()
		epv_ini_file = FortranFile(os.path.join(self.path,"epv_ini.ip"), "w")
		epv_ini_file.write_record(self.Ini.epv)
		epv_ini_file.close()
		sm_ini_file = FortranFile(os.path.join(self.path,"sm_ini.ip"), "w")
		sm_ini_file.write_record(self.Ini.sm)
		sm_ini_file.close()
		gok_file = FortranFile(os.path.join(self.path,"gok.ip"), "w")
		gok_file.write_record(self.Ini.gok)
		gok_file.close()
		k_file = FortranFile(os.path.join(self.path,"k.ip"), "w")
		k_file.write_record(self.Ini.k)
		k_file.close()
		chd_file = FortranFile(os.path.join(self.path,"chd.ip"), "w")
		chd_file.write_record(self.Ini.chd)
		chd_file.close()
		self.times.run["Fort.build"].append(time.time() - strt_time)
		logging.info("FORT: Built")
		if run:
			self.Run(load=True)

	def update(self, Delft, MF6, run=False):
		strt_time = time.time()
		self.Ini.sws = Delft.Res.sws[-1,:]*1000
		self.Ini.gws = MF6.Res.h[0,-1,:]*1000
		self.Ini.sm = self.Res.sm[:,-1]
		self.Ini.epv = self.Res.epv[:,-1]
		p_file = FortranFile(os.path.join(self.path,"p.ip"), "w")
		p_file.write_record(self.Ini.p)
		p_file.close()
		et_file = FortranFile(os.path.join(self.path,"et.ip"), "w")
		et_file.write_record(self.Ini.et)
		et_file.close()
		args = [self.Ini.elems, self.Ini.ts, self.Ini.ts_size, self.Ini.n, self.Ini.n_gw, self.Ini.m, self.Ini.beta,\
		self.Ini.alpha, self.Ini.sw_th]
		np.savetxt(os.path.join(self.path,"args.ip"), args)
		gws_ini_file = FortranFile(os.path.join(self.path,"gws_ini.ip"), "w")
		gws_ini_file.write_record(self.Ini.gws)
		gws_ini_file.close()
		sws_ini_file = FortranFile(os.path.join(self.path,"sws_ini.ip"), "w")
		sws_ini_file.write_record(self.Ini.sws)
		sws_ini_file.close()
		epv_ini_file = FortranFile(os.path.join(self.path,"epv_ini.ip"), "w")
		epv_ini_file.write_record(self.Ini.epv)
		epv_ini_file.close()
		sm_ini_file = FortranFile(os.path.join(self.path,"sm_ini.ip"), "w")
		sm_ini_file.write_record(self.Ini.sm)
		sm_ini_file.close()
		self.Ini.chd = np.zeros(self.Ini.elems, dtype=np.int8)
		for elem in range(self.Ini.elems):
			if elem in self.Ini.chd_cells:
				self.Ini.chd[elem] = 1
			else:
				self.Ini.chd[elem] = 0
		chd_file = FortranFile(os.path.join(self.path,"chd.ip"), "w")
		chd_file.write_record(self.Ini.chd)
		chd_file.close()
		self.times.run["Fort.update"].append(time.time() - strt_time)
		logging.info("FORT: Updated")
		if run:
			self.Run(load=True)

	def Run(self, load=False):
		strt_time = time.time()
		logging.info("FORT: Run initialized")
		owd = os.getcwd()
		if os.path.exists(self.exe_path):
			dir_path = os.path.split(self.exe_path)[0]
			exe = os.path.split(self.exe_path)[1]
			os.chdir(dir_path)
			fort_run = sp.Popen(exe, shell=True, stdout = sp.PIPE)
			stdout_fort, stderr_fort = fort_run.communicate()
			if not fort_run.returncode == 0:
				logging.error("FORT: Failed!")
				raise RuntimeError("FortRun failed")
			else:
				Fort.success = True
		else:
			logging.error("FORT: executable not found")
			raise NameError("FortRun executable missing")
		os.chdir(owd)
		self.times.run["Fort.Run"].append(time.time() - strt_time)
		logging.info("FORT: Run complete")
		if load:
			self.load()

	def load(self):
		strt_time = time.time()
		shape = (self.Ini.elems, self.Ini.ts)
		Qin_file = FortranFile(os.path.join(self.path,"Qin.op"), "r")
		self.Res.Qin = (Qin_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		Qin_file.close()
		Qout_file = FortranFile(os.path.join(self.path,"Qout.op"), "r")
		self.Res.Qout = (Qout_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		Qout_file.close()
		Qdiff_file = FortranFile(os.path.join(self.path,"Qdiff.op"), "r")
		self.Res.Qdiff = (Qdiff_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		Qdiff_file.close()
		gws_file = FortranFile(os.path.join(self.path,"gws.op"), "r")
		self.Res.gws = gws_file.read_reals().reshape(shape, order="F")
		gws_file.close()
		sws_file = FortranFile(os.path.join(self.path,"sws.op"), "r")
		self.Res.sws = sws_file.read_reals().reshape(shape, order="F")
		sws_file.close()
		sm_file = FortranFile(os.path.join(self.path,"sm.op"), "r")
		self.Res.sm = sm_file.read_reals().reshape(shape, order="F")
		sm_file.close()
		epv_file = FortranFile(os.path.join(self.path,"epv.op"), "r")
		self.Res.epv = epv_file.read_reals().reshape(shape, order="F")
		epv_file.close()
		gw_dis_file = FortranFile(os.path.join(self.path,"gw_dis.op"), "r")
		self.Res.gw_dis = (gw_dis_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		gw_dis_file.close()
		sw_dis_file = FortranFile(os.path.join(self.path,"sw_dis.op"), "r")
		self.Res.sw_dis = (sw_dis_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		sw_dis_file.close()
		sm_dis_file = FortranFile(os.path.join(self.path,"sm_dis.op"), "r")
		self.Res.sm_dis = (sm_dis_file.read_reals().reshape(shape, order="F"))/(self.times.dt["FortRun"]*1000)
		sm_dis_file.close()
		self.Res.max_err_abs = np.amax(self.Res.Qdiff.sum(1))
		self.Res.max_err_perc = np.amax(self.Res.Qdiff.sum(1))/self.Res.Qin.sum(1)[np.argmax(self.Res.Qdiff.sum(1))]*100
		self.times.run["Fort.load"].append(time.time() - strt_time)
		logging.info("FORT: Loaded results")

	def plot(self, elem, plotPrec=False, plotDis=False, savefig=False, dDPI=90, pDPI=1600, alpha_scatter=0.7, scatter_size=3, fromMF6=None):
		pal = ["#E74C3C", "#2ECC71", "#5EFCA1", "#E7AC3C", "#2980B9", "#1A3BFF", "#FF6600"] #[gw, sm, epv, sv, sw, p, et]

		if plotDis:
			plt.figure(dpi=dDPI)
			plt.xlabel("Time Steps")
			plt.ylabel("Discharges in Storage")
			plt.scatter(range(0,self.Ini.ts), self.Res.gw_dis[elem,:], label="GW_dis", color=pal[0],\
			alpha=alpha_scatter, s=scatter_size)
			plt.scatter(range(0,self.Ini.ts), self.Res.sm_dis[elem,:], label="SM_dis", color=pal[1],\
			alpha=alpha_scatter, s=scatter_size)
			plt.scatter(range(0,self.Ini.ts), self.Res.sw_dis[elem,:], label="SW_dis", color=pal[4],\
			alpha=alpha_scatter, s=scatter_size)
			plt.legend(loc="best", fontsize="small")
			if savefig:
				plt.savefig("prec_dist.jpg", format="jpeg", dpi=pDPI)

		plt.figure(dpi=dDPI)
		plt.xlabel("Time Steps")
		plt.ylabel("Water Levels in mm")
		if fromMF6:
			gws = fromMF6.Res.h[0,:,elem]*1000
		else:
			gws = self.Res.gws[elem,1:]
		plt.ylim([gws.min()-50, self.Res.sws[elem,:].max()+50+self.Ini.gok[elem]])
		plt.stackplot(range(0,self.Ini.ts-1), gws, self.Res.sm[elem,1:],\
		self. Res.epv[elem,1:]-self.Res.sm[elem,1:], (np.full(self.Ini.ts-1,self.Ini.gok[elem])-gws)*(1-self.Ini.n),\
		self.Res.sws[elem,1:], labels=["Groundwater","Soil Moisture", "Field Capacity", "Soil Volume", "Surface Water"], colors=pal)
		plt.plot(range(0,self.Ini.ts), np.full(self.Ini.ts,self.Ini.gok[elem]), "k", linewidth=0.5, label="Ground Level")
		plt.legend(loc="best", fontsize="small")
		if savefig:
			plt.savefig("water_levels.jpg", format="jpeg", dpi=pDPI)

		if plotPrec:
			plt.figure(dpi=dDPI)
			plt.xlabel("Time Steps")
			plt.ylabel("Precipitation in mm")
			plt.scatter(range(0,self.Ini.ts), self.Ini.p, s=0.1, color=pal[6])
			if savefig:
				plt.savefig("prec.jpg", format="jpeg", dpi=pDPI)

	def delIP(self):
		dat_files = [f for f in os.listdir(self.path) if f.endswith(".ip") ]
		for f in dat_files:
			os.remove(os.path.join(self.path, f))
	def delOP(self):
		dat_files = [f for f in os.listdir(self.path) if f.endswith(".op") ]
		for f in dat_files:
			os.remove(os.path.join(self.path, f))


#%%
class Delft:
	def __init__(self, path, times, Fort, OP_path="output/", mdu_file="wasenmoos.mdu", bat_file="run.bat", lat_path="", prepBC=False):
		self.path = path
		self.model_name = os.path.splitext(mdu_file)[0]
		self.OP_path_abs = os.path.join(self.path + OP_path)
		self.OP_path_rel = OP_path
		self.mapNC_path = os.path.join(self.OP_path_abs, self.model_name+"_map.nc")
		self.hisNC_path = os.path.join(self.OP_path_abs, self.model_name+"_his.nc")
		self.mdu_file = os.path.join(self.path, mdu_file)
		self.bat_path = os.path.join(self.path, bat_file)
		self.lat_path =  os.path.join(self.path, lat_path)
		if not os.path.exists(self.lat_path):
			os.mkdir(self.lat_path)
		self.lat_path_rel = lat_path
		self.lat_file_rel = os.path.join(lat_path, "laterals.ext")
		self.lat_file = os.path.join(self.lat_path, "laterals.ext")
		self.BC_dat_path = None
		self.BC_ext_path = None
		self.BC_tim_path = None
		self.prepBC = prepBC
		self.times = times
		self.Fort = Fort
		self.NC = bundleIT(elems=None, nodes=None, edges=None, maxVerts=None)
		self.MDU = bundleIT(IniFieldFile=None, ExtForceFileNew=None, mapint=None)
		self.Res = bundleIT()
		if self.MDU.IniFieldFile is None:
			self.MDU.IniFieldFile = "IniConds_rst.Ini"
		if self.MDU.ExtForceFileNew is None:
			self.MDU.ExtForceFileNew = self.lat_file_rel
		if prepBC:
			if not self.BC_dat_path:
				self.BC_dat_path = "../data/delft_grabenBC.dat"
			if not self.BC_tim_path:
				self.BC_tim_path = os.path.join(self.path + "grabenBC_0001.tim")
			if not self.BC_ext_path:
				self.BC_ext_path = "grabenBC.ext"
			bcDF = pd.read_table(self.BC_dat_path)
			bcDF[bcDF.columns[0]] = pd.to_datetime(bcDF[bcDF.columns[0]])
			bcDF.index = bcDF[bcDF.columns[0]]
			bcDF = bcDF.drop(bcDF.columns[0], 1)
			start = datetime(self.times.stmp["global_start"].year, self.times.stmp["global_start"].month, self.times.stmp["global_start"].day,\
			self.times.stmp["global_start"].hour)
			end = datetime(self.times.stmp["global_end"].year, self.times.stmp["global_end"].month, self.times.stmp["global_end"].day,\
			self.times.stmp["global_end"].hour)
			bc_mask = (bcDF.index >= start) & (bcDF.index <= end)
			bcDF = bcDF.loc[bc_mask]
			self.bcDF = bcDF.drop_duplicates()
			query_time = self.times.stmp["local_start"] + (self.times.stmp["local_end"] - self.times.stmp["local_start"])/2
			idx = self.bcDF.index.get_loc(query_time , method="nearest")
			self.bc_val = self.bcDF[self.bcDF.keys()[0]][idx]
			self.bc_val = np.repeat(self.bc_val,2)
			logging.info("Delft: BC DF Prepped")

	def readNC(self, nc_path=None, v=4, internal=False, lazy=False, lazy_path="../data/NCread.npz"):
		strt_time = time.time()
		if not lazy:
			if nc_path is not None:
				self.mapNC_path = nc_path
			if v == 4:
				logging.info("DELFT: Reading mapNC from " + self.mapNC_path + "(v4)")
				mapNC = nc.Dataset(self.mapNC_path)
				nelems = mapNC.dimensions["mesh2d_nFaces"].size
				elemX = np.array(mapNC["mesh2d_face_x"][:])
				elemY = np.array(mapNC["mesh2d_face_y"][:])
				elemZ = np.array(mapNC["mesh2d_flowelem_bl"][:])
				area = np.array(mapNC["mesh2d_flowelem_ba"][:])
				ma_array = mapNC["mesh2d_face_nodes"][:]
				elems = bundleIT(n=nelems, x=elemX, y=elemY, z=elemZ, area=area, ma_array=ma_array)
				self.NC.elems = elems
				nnodes = mapNC.dimensions["mesh2d_nNodes"].size
				nodeX = np.array(mapNC["mesh2d_node_x"][:])
				nodeY = np.array(mapNC["mesh2d_node_y"][:])
				nodeZ = np.array(mapNC["mesh2d_node_z"][:])
				nodes = bundleIT(n=nnodes, x=nodeX, y=nodeY, z=nodeZ)
				self.NC.nodes = nodes
				self.NC.edges = mapNC.dimensions["mesh2d_nEdges"].size
				self.NC.maxVerts = mapNC.dimensions["mesh2d_nMax_face_nodes"].size
				mapNC.close()
				if not internal:
					self.Fort.Ini.elems = nelems
					self.Fort.Ini.gok = elemZ*1000
				logging.info("DELFT: Read")
			elif v == 3.0:
				logging.info("DELFT: Reading mapNC from " + self.mapNC_path + "(v3)")
				mapNC = nc.Dataset(self.mapNC_path)
				nelems = mapNC.dimensions["nFlowElem"].size
				elemX = np.array(mapNC["FlowElem_xcc"][:])
				elemY = np.array(mapNC["FlowElem_ycc"][:])
				elemZ = np.array(mapNC["FlowElem_bl"][:])
				area = np.array(mapNC["FlowElem_bac"][:])
				ma_array = mapNC["NetElemNode"][:]
				elems = bundleIT(n=nelems, x=elemX, y=elemY, z=elemZ, area=area, ma_array= ma_array)
				self.NC.elems = elems
				nnodes = mapNC.dimensions["nNetNode"].size
				nodeX = np.array(mapNC["NetNode_x"][:])
				nodeY = np.array(mapNC["NetNode_y"][:])
				nodeZ = np.array(mapNC["NetNode_z"][:])
				nodes = bundleIT(n=nnodes, x=nodeX, y=nodeY, z=nodeZ)
				self.NC.nodes = nodes
				self.NC.edges = mapNC.dimensions["nNetLink"].size
				self.NC.maxVerts = mapNC.dimensions["nNetElemMaxNode"].size
				if not internal:
					self.Fort.Ini.elems = nelems
					self.Fort.Ini.gok = elemZ*1000
				mapNC.close()
				logging.info("DELFT: Read")
			elif v == 3.1:
				logging.info("DELFT: Reading mapNC from " + self.mapNC_path + "(v3.1)")
				mapNC = nc.Dataset(self.mapNC_path)
				nelems = mapNC.dimensions["nFlowElem"].size
				elemX = np.array(mapNC["FlowElem_xzw"][:])
				elemY = np.array(mapNC["FlowElem_yzw"][:])
				elemZ = np.array(mapNC["FlowElem_bl"][:])
				area = np.array(mapNC["FlowElem_bac"][:])
				ma_array = mapNC["NetElemNode"][:]
				elems = bundleIT(n=nelems, x=elemX, y=elemY, z=elemZ, area=area, ma_array= ma_array)
				self.NC.elems = elems
				nnodes = mapNC.dimensions["nNetNode"].size
				nodeX = np.array(mapNC["NetNode_x"][:])
				nodeY = np.array(mapNC["NetNode_y"][:])
				nodeZ = np.array(mapNC["NetNode_z"][:])
				nodes = bundleIT(n=nnodes, x=nodeX, y=nodeY, z=nodeZ)
				self.NC.nodes = nodes
				self.NC.edges = mapNC.dimensions["nNetLink"].size
				self.NC.maxVerts = mapNC.dimensions["nNetElemMaxNode"].size
				if not internal:
					self.Fort.Ini.elems = nelems
					self.Fort.Ini.gok = elemZ*1000
				np.savez("../data/NCread", nelems=nelems,elemsX=elemX,elemsY=elemY,elemsZ=elemZ,area=area,\
				nnodes=nnodes,nodeX=nodeX,nodeY=nodeY,nodeZ=nodeZ, edges=self.NC.edges, maxVerts=self.NC.maxVerts)
				ma_array.dump("../data/NCread")
				mapNC.close()
				logging.info("DELFT: Read")
			elif v == 0:
				logging.info("DELFT: Reading netNC from " + self.mapNC_path)
				netNC = nc.Dataset(self.mapNC_path)
				nelems = mapNC.dimensions["mesh2d_nFaces"].size
				elemX = np.array(mapNC["mesh2d_face_x"][:])
				elemY = np.array(mapNC["mesh2d_face_y"][:])
				elemZ = None
				area = None
				elems = bundleIT(n=nelems, x=elemX, y=elemY, z=elemZ, area=area)
				self.NC.elems = elems
				nnodes = mapNC.dimensions["mesh2d_nNodes"].size
				nodeX = np.array(mapNC["mesh2d_node_x"][:])
				nodeY = np.array(mapNC["mesh2d_node_y"][:])
				nodeZ = np.array(mapNC["mesh2d_node_z"][:])
				nodes = bundleIT(n=nnodes, x=nodeX, y=nodeY, z=nodeZ)
				self.NC.nodes = nodes
				self.NC.edges = mapNC.dimensions["mesh2d_nEdges"].size
				self.NC.maxVerts = mapNC.dimensions["mesh2d_nMax_face_nodes"].size
				if not internal:
					self.Fort.Ini.elems = nelems
				netNC.close()
				logging.info("DELFT: Read")
			else:
				logging.error("DELFT: NCread: Invalid version number")
				raise ValueError("Invalid version number")
		else:
			lazy_file = np.load(lazy_path)
			logging.info("DELFT: Reading lazy NC file from " + lazy_path)
			self.NC.elems = bundleIT()
			self.NC.elems.n = int(lazy_file["nelems"])
			self.NC.elems.x = lazy_file["elemsX"]
			self.NC.elems.y = lazy_file["elemsY"]
			self.NC.elems.z = lazy_file["elemsZ"]
			self.NC.elems.area = lazy_file["area"]
			self.NC.nodes = bundleIT()
			self.NC.nodes.n = int(lazy_file["nnodes"])
			self.NC.nodes.x = lazy_file["nodeX"]
			self.NC.nodes.y = lazy_file["nodeY"]
			self.NC.nodes.z = lazy_file["nodeZ"]
			self.NC.edges = lazy_file["edges"]
			self.NC.maxVerts = lazy_file["maxVerts"]
			ma_path = os.path.join(os.path.split(lazy_path)[0]+"/", os.path.splitext(os.path.split(lazy_path)[1])[0])
			self.NC.elems.ma_array = np.load(ma_path, allow_pickle=True)
			if not internal:
				self.Fort.Ini.elems = self.NC.elems.n
				self.Fort.Ini.gok = self.NC.elems.z*1000
			logging.info("DELFT: Read")
		if not internal and v!=3.1:
			self.load(internal=True)
			self.Fort.Ini.sws = self.Res.sws[-1,:]*1000
			logging.info("DELFT: Passed ini SW to FORT")
			self.times.run["Delft.readNC"].append(time.time() - strt_time)

	def update(self, run=False):
		strt_time = time.time()
		if self.prepBC:
			query_time = self.times.stmp["local_start"] + (self.times.stmp["local_end"] - self.times.stmp["local_start"])/2
			idx = self.bcDF.index.get_loc(query_time , method="nearest")
			self.bc_val = self.bcDF[self.bcDF.keys()[0]][idx]
			self.bc_val = np.repeat(self.bc_val,2)
			tim_array = np.stack([self.times.rel["delft_tim"], self.bc_val], axis=1)
			np.savetxt(self.BC_tim_path,tim_array)
			logging.info("DELFT: grabenBC tim file updated")
		self.sw_flux = np.zeros(self.Fort.Res.sw_dis.shape)
		for i in range(self.sw_flux.shape[1]):
			self.sw_flux[:,i] = np.multiply(self.Fort.Res.sw_dis[:,i], self.NC.elems.area.T).T
		sw_dis = self.sw_flux
		fHEADER = "[General]\nfileVersion = 2.01\nfileType = extForce\n\n\n"
		lHEADER = "[Lateral]\n"
		ID = "id = l{0}\n"
		NAME = "name = lat_{0}\nnumCoordinates = 1\n"
		X = "xCoordinates = {0}\n"
		Y = "yCoordinates = {0}\n"
		DIS = "discharge = {0}\n\n"
		ext_file = open(os.path.join(self.path, self.MDU.ExtForceFileNew), "w")
		ext_file.writelines(fHEADER)
		ext_data = ""
		for i in range(self.NC.elems.n):
			if sw_dis[i,:].any() != 0:
				dis_avg = sw_dis[i,1:].mean()
				ext_data += lHEADER + ID.format(i) + NAME.format(i) + \
				X.format(self.NC.elems.x[i]) + Y.format(self.NC.elems.y[i]) + DIS.format(dis_avg)
		ext_file.write(ext_data)
		ext_file.close()
		logging.info("DELFT: Update: Wrote sw_dis from FORT into EXT file")
		self.times.run["Delft.update"].append(time.time() - strt_time)
		if run:
			self.Run(load=True)

	def updateTIM(self, run=False):
		strt_time = time.time()
		if self.prepBC:
			query_time = self.times.stmp["local_start"] + (self.times.stmp["local_end"] - self.times.stmp["local_start"])/2
			idx = self.bcDF.index.get_loc(query_time , method="nearest")
			self.bc_val = self.bcDF[self.bcDF.keys()[0]][idx]
			self.bc_val = np.repeat(self.bc_val,2)
			tim_array = np.stack([self.times.rel["delft_tim"], self.bc_val], axis=1)
			np.savetxt(self.BC_tim_path,tim_array)
			logging.info("DELFT: grabenBC tim file updated")
		if not os.path.exists(self.lat_path):
			os.mkdir(self.lat_path)
		self.sw_flux = np.zeros(self.Fort.Res.sw_dis.shape)
		for i in range(self.sw_flux.shape[1]):
			self.sw_flux[:,i] = np.multiply(self.Fort.Res.sw_dis[:,i], self.NC.elems.area.T).T
		sw_dis = self.sw_flux
		tim_times = np.linspace(self.times.rel["delft_tim"][0], self.times.rel["delft_tim"][1], num=self.times.nTS["exchange"]+1)
		fHEADER = "[General]\nfileVersion = 2.01\nfileType = extForce\n\n"
		lHEADER = "[Lateral]\n"
		ID = "id = l{0}\n"
		NAME = "name = lat_{0}\nnumCoordinates = 1\n"
		X = "xCoordinates = {0}\n"
		Y = "yCoordinates = {0}\n"
		DIS = "discharge = " + self.lat_path_rel + "l{n}.tim\n\n"
		ext_file = open(self.lat_file, "w")
		ext_file.writelines(fHEADER+"\n")
		ext_data = ""
		for i in range(self.NC.elems.n):
			if sw_dis[i,:].any() != 0:
				ext_data += lHEADER + ID.format(i) + NAME.format(i) + \
				X.format(self.NC.elems.x[i]) + Y.format(self.NC.elems.y[i]) + DIS.format(n=i)
				tim_file = open(os.path.join(self.lat_path, "l{0}.tim".format(i)), "w")
				tim_data = ""
				for j in range(tim_times.size):
					tim_data += str(tim_times[j]) + "\t" + str(sw_dis[i,j]) + "\n"
				tim_file.writelines(tim_data)
				tim_file.close()
		ext_file.write(ext_data)
		ext_file.close()
		logging.info("DELFT: Update: Wrote sw_dis from FORT into EXT and TIM files")
		self.times.run["Delft.update"].append(time.time() - strt_time)
		if run:
			self.Run(load=True)

	def Run(self, initialRun=False, firstRestart=False, load=False, update=False, updateTIM=False):
		strt_time = time.time()
		self.err = False
		if update:
			self.update()
		if updateTIM:
			self.updateTIM()
		if not initialRun:
			self.writeRestart(self.times, firstRestart=firstRestart)
		owd = os.getcwd()
		if os.path.exists(self.bat_path):
			dir_path = os.path.split(self.bat_path)[0]
			bat_path = os.path.split(self.bat_path)[1]
			os.chdir(dir_path)
			logging.info("DELFT: Run initialized")
			delft_run = sp.Popen(os.path.abspath(bat_path), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
			stdout, stderr = delft_run.communicate()
			delft_cli_dump_check = "process 0" in stderr.decode(sys.stdout.encoding)
			if delft_run.returncode != 0 or delft_cli_dump_check:
				logging.error("DELFT: Run Failed. Check DIA file\n\n\n")
				raise RuntimeError("DelftRun failed. Check DIA file for more info", stderr)
				self.err = True
			else:
				Delft.success = True
		else:
			logging.error("DELFT: BAT file not found")
			raise ValueError("Incorrect BAT file")
		os.chdir(owd)
		logging.info("DELFT: Run complete")
		self.times.run["Delft.Run"].append(time.time() - strt_time)
		if load:
			self.load()

	def writeRestart(self, Restart_from="rst", firstRestart=False):
		for line in fileinput.input(self.mdu_file, inplace=True):
			rstrt_file = self.OP_path_rel + self.model_name + "_" + self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S") + "_rst.nc"
			date = self.times.stmp["local_start"].strftime("%Y%m%d")
			if self.times.dt["FortRun"] == self.times.dt["exchange"]:
				dtmax = self.times.dt["FortRun"]/3
			else:
				dtmax = self.times.dt["FortRun"]
			if line.startswith("RestartFile"):
				print("RestartFile\t\t\t\t= " + rstrt_file)
			elif line.startswith("MapInterval"):
				print("MapInterval\t\t\t\t=", int(self.times.dt["FortRun"]))
			elif line.startswith("HisInterval"):
				print("HisInterval\t\t\t\t=", int(self.times.dt["FortRun"]))
			elif line.startswith("RstInterval"):
				print("RstInterval\t\t\t\t=", int(self.times.dt["exchange"]))
			elif Restart_from == "map" and line.startswith("RestartDateTime"):
				print("RestartDateTime\t\t\t= " + self.times.stmp["local_start"].strftime("%Y%m%d%H%M%S"))
			elif line.startswith("RefDate") and firstRestart:
				print("RefDate\t\t\t\t\t=", date)
			elif line.startswith("TStart"):
				print("TStart\t\t\t\t\t=", self.times.rel["delft"][0])
			elif line.startswith("TStop"):
				print("TStop\t\t\t\t\t=", self.times.rel["delft"][1])
			elif line.startswith("DtUser"):
				print("DtUser\t\t\t\t\t=", dtmax)
			elif line.startswith("DtMax"):
				print("DtMax\t\t\t\t\t=", dtmax)
			elif line.startswith("Dtinit"):
				print("Dtinit\t\t\t\t\t=", 0.001)
			elif firstRestart and line.startswith("IniFieldFile"):
				print("IniFieldFile\t\t\t=", self.MDU.IniFieldFile)
			elif firstRestart and line.startswith("ExtForceFileNew"):
				print("ExtForceFileNew\t\t\t=", self.MDU.ExtForceFileNew)
			elif self.prepBC and firstRestart and line.startswith("ExtForceFile"):
				print("ExtForceFile\t\t\t=", self.BC_ext_path)
			elif firstRestart and line.startswith("RefDate"):
				print("RefDate\t\t\t\t\t=", self.times.stmp["global_start"].strftime("%Y%m%d"))
			else:
				sys.stdout.write(line)
		fileinput.close()
		logging.info("DELFT: MDU file modified for restart")

	def load(self, internal=False):
		strt_time = time.time()
		map_file = self.path + self.OP_path_abs + self.model_name + "_map.nc"
		mapNC = nc.Dataset(map_file)
		self.Res.sws_abs = np.array(mapNC["mesh2d_s1"][:])
		self.Res.sws = np.array(mapNC["mesh2d_waterdepth"][:])
		self.Res.sws_prev = np.array(mapNC["mesh2d_s0"][:])
		mapNC.close()
		logging.info("DELFT: Loaded results")
		if not internal:
			self.times.run["Delft.load"].append(time.time() - strt_time)

	def plot(self, MF6, ts=-1):
		fig = plt.figure()
		fig.suptitle("Water Level in m")
		pmv = fp.plot.PlotMapView(model=MF6.model)
		pmv.plot_grid(linewidth=0.1)
		hplot = pmv.plot_array(self.Res.sws[ts,:], cmap=mpl.cm.Blues)
		plt.colorbar(hplot)

#%%
class ModFlow():
	def __init__(self, path, times, Delft, Fort, exe_path="../../MF6/res/mf6.exe"):
		self.path = path
		self.exe_path = exe_path
		self.times = times
		self.Delft = Delft
		self.Fort = Fort
		self.name = self.Delft.model_name
		self.head_file = "{}.hds".format(self.name)
		self.Params = bundleIT()
		self.Params.internalTS = 1
		self.Res = bundleIT()
		self.boundary_elems_file = "../data/boundaries.npz"
		boundaries = np.load(self.boundary_elems_file)
		self.Params.ridge_elems = boundaries["ridge_boundary_elems"] - 1
		self.Params.graben_elems = boundaries["graben_boundary_elems"] - 1

	def build(self, run=False, lazy=False, restart=False, res=None):
		start_time = time.time()
		if lazy:
			self.sim = fp.mf6.MFSimulation.load(sim_name=self.name, exe_name=self.exe_path, version="mf6", sim_ws=self.path, verbosity_level=0)
			self.tdis = self.sim.tdis
			self.ims = self.sim.ims
			self.model = self.sim.get_model(self.name)
			self.disv = self.model.disv
			self.npf = self.model.npf
			self.ic = self.model.ic
			if restart:
				head_file = self.name+"_"+(self.times.stmp["local_start"]-timedelta(seconds=self.times.dt["exchange"]))\
				.strftime("%Y%m%d_%H%M%S").strftime("%Y%m%d_%H%M%S")+".hds"
				if res:
					head_path = os.path.join(res.mf_path, head_file)
				else:
					head_path = os.path.join(self.path, head_file)
				hds = fp.utils.binaryfile.HeadFile(head_path)
				start = np.zeros((self.disv.nlay.data, self.Delft.NC.elems.n))
				start = hds.get_data(kstpkper=(self.Params.internalTS-1, len(hds.get_kstpkper())-1)).mean(axis=1)
				self.ic.strt = start
				self.ic.write()
				self.Fort.Ini.gws = self.ic.strt.data[0,:]*1000
				logging.info("MF: Copied initial heads from previous run and passed to FORT as well")
			self.sto = self.model.sto
			self.chd = self.model.chd
			self.rch = self.model.rch
			self.oc = self.model.oc
			self.obs = self.model.obs
			self.Fort.Ini.gws = self.ic.strt.data[0,:]*1000
			self.Fort.Ini.k = self.npf.k.data[0,:]*1000
			self.Params.top = self.disv.top.data
			logging.info("MF: Built lazily. k and gws_ini values passed to FORT")

		else:
			self.Params.nlay = 3
			self.Params.ncells = self.Delft.NC.elems.n
			self.Params.top = self.Delft.NC.elems.z
			self.Params.aq_thickness = np.zeros((self.Params.nlay, self.Params.ncells))
			self.Params.aq_thickness[0] = None
			self.Params.aq_thickness[1] = 1
			self.Params.aq_thickness[2] = 10
			self.Params.k = np.zeros((self.Params.nlay, self.Params.ncells))
			self.Params.k[0] = 5e-4
			self.Params.k[1] = 3e-9
			self.Params.k[2] = 1e-4
			self.Params.ini_h = [None]*3
			self.Params.ini_h[0] = ["ratio", 0.8]
			self.Params.ini_h[1] = ["rel", 0.1]
			self.Params.ini_h[2] = ["rel", 0.1]
			self.Params.sy = np.zeros(self.Params.nlay)
			self.Params.sy[0] = 0.4
			self.Params.sy[1] = 0.15
			self.Params.sy[2] = 0.15
			self.Params.ss = np.zeros(self.Params.nlay)
			self.Params.ss[0] = 6e-3
			self.Params.ss[1] = 1e-3
			self.Params.ss[2] = 2e-5
			#
			def flatten(xs):
				result = []
				if isinstance(xs, (list, tuple)):
					for x in xs:
						result.extend(flatten(x))
				else:
							result.append(xs)
				return result
			#
			def getDisV():
				xor = self.Delft.NC.nodes.x.min()
				yor = self.Delft.NC.nodes.y.min()
				nodes = self.Delft.NC.nodes.n
				node_index = np.array(range(0,nodes))
				nodeX = self.Delft.NC.nodes.x - xor
				nodeY = self.Delft.NC.nodes.y - yor
				verteces = []
				for n in node_index:
					verteces.append([node_index[n], nodeX[n], nodeY[n]])
				elements = self.Params.ncells
				elem_index = np.array(range(elements))
				elemental_nodes_ma = self.Delft.NC.elems.ma_array - 1
				elemental_nodes_ma = np.flip(elemental_nodes_ma, axis=1)
				elemental_nodes = [en.compressed().tolist() for en in elemental_nodes_ma]
				elemX = self.Delft.NC.elems.x - xor
				elemY = self.Delft.NC.elems.y - yor
				cell2d = []
				for e in elem_index:
					temp = [elem_index[e], elemX[e], elemY[e], len(elemental_nodes[e]), elemental_nodes[e]]
					temp = flatten(temp)
					cell2d.append(temp)
				return verteces, cell2d, elements, nodes, elem_index
			#
			self.sim = fp.mf6.MFSimulation(sim_name=self.name, exe_name=self.exe_path, version="mf6", sim_ws=self.path)
			#
			period_data = []
			for i in range(self.times.nTS["exchange"]):
				period_data.append((self.times.dt["FortRun"], self.Params.internalTS, 1))
			self.tdis = fp.mf6.ModflowTdis(self.sim, time_units="seconds", nper=len(period_data), perioddata=period_data)
			#
			self.ims = fp.mf6.ModflowIms(self.sim, complexity="complex", no_ptcrecord="all", linear_acceleration="BICGSTAB",\
			outer_dvclose=1e-7, outer_rclosebnd=1e-7, outer_maximum=10000, relaxation_factor=0.97,\
			under_relaxation="dbd", under_relaxation_gamma=0.25, under_relaxation_theta=0.7, under_relaxation_kappa=0.1,\
			backtracking_number=10, backtracking_tolerance=10e4, backtracking_reduction_factor=0.2, backtracking_residual_limit=100,\
			scaling_method="L2NORM", reordering_method="MD")
			self.sim.register_ims_package(self.ims, [self.name])
			#
			self.model = fp.mf6.ModflowGwf(self.sim, modelname=self.name, model_nam_file="{}.nam".format(self.name),\
			newtonoptions="under_relaxation")
			#
			verteces, cell2d, elements, nodes, self.Params.elem_index = getDisV()
			bottom_elev = np.ones((self.Params.nlay, elements))
			bottom_elev[0] = np.full(elements, self.Delft.NC.elems.z.min()- 1e-7)
			bottom_elev[1] = bottom_elev[0] - self.Params.aq_thickness[1]
			bottom_elev[2] = bottom_elev[1] - self.Params.aq_thickness[2]
			self.disv = fp.mf6.ModflowGwfdisv(self.model, length_units="meters", nlay=self.Params.nlay,\
			top=self.Params.top, botm=bottom_elev, ncpl=elements,\
			nvert=nodes, vertices=verteces, cell2d=cell2d)
			self.Params.aq_volume =  np.zeros((self.Params.nlay,elements))
			self.Params.aq_thickness[0] = bottom_elev[0] - self.Params.top
			self.Params.aq_volume[0] = np.multiply(self.Params.aq_thickness[0], self.Delft.NC.elems.area.T)
			self.Params.aq_volume[1] = np.multiply(self.Params.aq_thickness[1], self.Delft.NC.elems.area.T)
			self.Params.aq_volume[2] = np.multiply(self.Params.aq_thickness[2], self.Delft.NC.elems.area.T)
			#
			ict = np.ones((self.Params.nlay,elements))
			ict[1] = -1
			ict[2] = -1
			k = np.ones((self.Params.nlay,elements))
			k[0] = self.Params.k[0]
			self.Fort.Ini.k = k[0,:]*1000
			k[1] = self.Params.k[1]
			k[2] = self.Params.k[2]
			self.npf = fp.mf6.ModflowGwfnpf(self.model, icelltype=ict, k=k, save_flows=True)
			#
			if restart:
				head_file = self.name+"_"+(self.times.stmp["local_start"]-timedelta(seconds=self.times.dt["exchange"]))\
				.strftime("%Y%m%d_%H%M%S")+".hds"
				if res:
					head_path = os.path.join(res.mf_path, head_file)
				else:
					head_path = os.path.join(self.path, head_file)
				hds = fp.utils.binaryfile.HeadFile(head_path)
				start = np.zeros((self.disv.nlay.data, self.Params.ncells))
				start = hds.get_data(kstpkper=(self.Params.internalTS-1, len(hds.get_kstpkper())-1)).mean(axis=1)
				self.ic = fp.mf6.ModflowGwfic(self.model, strt=start)
				self.Fort.Ini.gws = self.ic.strt.data[0,:]*1000
				logging.info("MF: Copied initial heads from previous run and passed to FORT as well")
			else:
				start = np.ones((self.Params.nlay,elements))
				start[0,:] = bottom_elev[0,:]
				if self.Params.ini_h[0][0] == "ratio":
					start[0,:] = (self.Params.top[:] - bottom_elev[0,:])*self.Params.ini_h[0][1] + bottom_elev[0,:]
				elif self.Params.ini_h[0][0] == "abs":
					start[0,:] = np.array(self.Params.ini_h[0][1])
				else:
					raise ValueError("Type"+self.Params.ini_h[0][0]+"invalid")
				self.Fort.Ini.gws = start[0,:]*1000
				start[1,:] = bottom_elev[1,:]
				if self.Params.ini_h[1][0] == "rel":
					start[1,:] = self.Params.top + self.Params.ini_h[1][1]
				elif self.Params.ini_h[1][0] == "abs":
					start[1,:] = np.array(self.Params.ini_h[1][1])
				else:
					raise ValueError("Type"+self.Params.ini_h[1][0]+"invalid")
				start[2,:] = bottom_elev[2,:]
				if self.Params.ini_h[2][0] == "rel":
					start[2,:] = self.Params.top + self.Params.ini_h[2][1]
				elif self.Params.ini_h[2][0] == "abs":
					start[2,:] = np.array(self.Params.ini_h[2][1])
				else:
					raise ValueError("Type"+self.Params.ini_h[2][0]+"invalid")
				self.ic = fp.mf6.ModflowGwfic(self.model, strt=start)
			#
			ss = np.ones((self.Params.nlay,elements))
			ss[0,:] = self.Params.ss[0]
			ss[1,:] = self.Params.ss[1]
			ss[2,:] = self.Params.ss[2]
			sy = np.ones((self.Params.nlay,elements))
			sy[0,:] = self.Params.sy[0]
			sy[1,:] = self.Params.sy[1]
			sy[2,:] = self.Params.sy[2]
			icv = np.ones((self.Params.nlay,elements))
			icv[1:3] = 0
			self.sto = fp.mf6.ModflowGwfsto(self.model, save_flows=True, iconvert=icv, ss=ss, sy=sy)
			#
			chd_sp_data = []
			self.chd_cells = []
			for i in self.Params.ridge_elems:
				chd_sp_data.append([2, int(i), self.Delft.NC.elems.z[int(i)]+1])
			if not self.Delft.prepBC:
				graben_head = self.Delft.NC.elems.z.min()+0.2
			else:
				graben_head = self.Delft.bc_val.mean()
			for elem in self.Params.graben_elems:
				if graben_head >= self.Params.top[elem]:
					chd_sp_data.append([0, int(elem), graben_head])
					self.chd_cells.append(elem)
			self.chd = fp.mf6.ModflowGwfchd(self.model, maxbound=len(chd_sp_data), stress_period_data=chd_sp_data, save_flows=True)
			self.Fort.Ini.chd_cells = self.chd_cells
			#
			self.rch = fp.mf6.ModflowGwfrch(self.model, maxbound=self.Params.ncells, stress_period_data=0, save_flows=True)
			#
			headrec = [self.head_file]
			budgetrec = ["{}.cbb".format(self.name)]
			saverec = [("HEAD", "LAST"), ("BUDGET", "LAST")]
			printrec = [("HEAD", "LAST")]
			self.oc = fp.mf6.ModflowGwfoc(self.model, saverecord=saverec, head_filerecord=headrec,\
			budget_filerecord=budgetrec, printrecord=printrec)
			#
			self.Params.obs_pts = [("P7","head",1,303) ,("P6","head",1,36), ("PGraben","head",1,273), ("P12","head",1,619),\
			("P2","head",1,2183), ("P1","head",1,1993)]
			self.obs = fp.mf6.ModflowUtlobs(self.model, continuous=self.Params.obs_pts, digits=12)
			#
			self.sim.write_simulation(silent=True)
			logging.info("MF: Built actively. k and gws_ini values passed to FORT")
		self.times.run["MF.build"].append(time.time() - start_time)
		if run:
			self.Run(load=True)

	def update(self, run=False):
		start_time = time.time()
		def prepRch(gw_dis, elem_index):
			sp_data = {}
			for t in range(gw_dis.shape[1]):
				data = []
				for e in range(gw_dis.shape[0]):
					data.append([0, elem_index[e], gw_dis[e,t]])
				sp_data[t] = {"filename": "gw_dis_sp{0}.dat".format(t), "data": data}
			return sp_data
		#
		if self.Delft.prepBC:
			chd_sp_data = []
			self.chd_cells = []
			for elem in self.Params.ridge_elems:
				chd_sp_data.append([2, int(elem), self.Delft.NC.elems.z[int(elem)]+1])
			graben_head = self.Delft.Res.sws_abs[-1,:]
			for elem in self.Params.graben_elems:
				if graben_head[elem] >= self.Params.top[elem]:
					chd_sp_data.append([0, int(elem), graben_head[elem]])
					self.chd_cells.append(elem)
			self.Fort.Ini.chd_cells = self.chd_cells
			self.chd = fp.mf6.ModflowGwfchd(self.model, maxbound=len(chd_sp_data), stress_period_data=chd_sp_data, save_flows=True)
			self.chd.write()
		#
		period_data = []
		for i in range(self.times.nTS["exchange"]):
			period_data.append((self.times.dt["FortRun"], self.Params.internalTS, 1))
		self.tdis.nper = len(period_data)
		self.tdis.perioddata = period_data
		self.tdis.write()
		#
		if os.path.exists(os.path.join(self.path, self.name+".rch")): os.remove(os.path.join(self.path, self.name+".rch"))
		files = os.listdir(self.path)
		for item in files:
			if item.endswith(".dat"):
				os.remove(os.path.join(self.path, item))
		sp_data = prepRch(self.Fort.Res.gw_dis[:,1:], np.arange(self.Params.ncells))
		self.rch = fp.mf6.ModflowGwfrch(self.model, maxbound=self.Params.ncells, stress_period_data=sp_data, save_flows=True)
		self.rch.write()
		self.ic.strt = self.Res.h[:,-1,:]
		self.ic.write()
		logging.info("MF: Initial heads copied from previous sim and gw_dis written as RCH")
		self.times.run["MF.update"].append(time.time() - start_time)
		logging.info("MF: gw_dis written as RCH")
		if run:
			self.Run(load=True)

	def Run(self, silent=True, rewrite=False, load=False, update=False, err_rerun=False):
		start_time = time.time()
		if rewrite:
			self.sim.write_simulation(silent=silent)
		if update:
			self.update()
		logging.info("MF: Run initialized")
		self.success, buff = self.sim.run_simulation(silent=silent)
		if not self.success:
			logging.error("MF: Run failed. Check LST file")
			if err_rerun:
				self.Run(silent=False, err_rerun=False)
			err_rerun = False
			raise RuntimeError("ModFlow Run failed. Check LST file for more info")
		self.times.run["MF.Run"].append(time.time() - start_time)
		logging.info("MF: Run complete")
		if load:
			self.load()

	def load(self):
		start_time = time.time()
		self.hds = fp.utils.binaryfile.HeadFile(os.path.join(self.path,self.head_file))
		self.Res.h = np.zeros((self.disv.nlay.data, self.times.nTS["exchange"], self.Params.ncells))
		for i in range(len(self.hds.get_kstpkper())):
			self.Res.h[:,i,:] = self.hds.get_data(kstpkper=(self.Params.internalTS-1, i)).mean(axis=1)
		self.times.run["MF.load"].append(time.time() - start_time)
		self.cbb = fp.utils.binaryfile.CellBudgetFile(os.path.join(self.path,self.name+".cbb"))
		self.lst = fp.utils.mflistfile.Mf6ListBudget(os.path.join(self.path,self.name+".lst"))
		logging.info("MF: Loaded results")

	def plot(self, layer=1, ts=-1):
		fig = plt.figure()
		fig.tight_layout()
		fig.suptitle("heads for layer "+str(layer))
		pmv = fp.plot.PlotMapView(model=self.model, layer=layer-1)
		pmv.plot_grid(linewidth=0.1)
		hplot = pmv.plot_array(self.Res.h[:,ts,:], cmap=mpl.cm.turbo_r)
		plt.colorbar(hplot)


#%%
class resNC:
	def __init__(self, times, Fort, Delft, MF6, op_path="../output/", op_name="GWSWEX"):
		self.path = op_path
		if not os.path.exists(self.path):
			os.mkdir(op_path)
		self.delft_path = os.path.join(self.path, "delft/")
		self.mf_path = os.path.join(self.path, "modflow/")
		self.fort_path = os.path.join(self.path, "fortran/")
		self.name = op_name
		self.times = times
		self.Fort = Fort
		self.Delft = Delft
		self.MF6 = MF6
		self.save = None
		self.plt_path = "../../../Desktop"

	def build(self):
		self.OP = nc.Dataset(os.path.join(self.path, self.name+".nc"), "w")
		self.OP.source = "GWSWEX: Groundwater-Surfacewater Exchange Model"
		self.OP.createDimension("element", self.Delft.NC.elems.n)
		self.OP.createDimension("time", None)

		timeV = self.OP.createVariable("time", "i", ("time", ))
		timeV.units = "seconds since " + self.times.str["global_start"]
		time_unixV = self.OP.createVariable("time_unix", "i", ("time", ))
		time_unixV.units = "unix time, seconds"
		runtimeV = self.OP.createVariable("runtime", "i", ("time", ))
		runtimeV.units = "seconds"
		nTSV = self.OP.createVariable("nTS", "i", ("time", ))
		runtimeV[0] = 0
		nTSV[0] = 0

		areaV =  self.OP.createVariable("area", "i", ("element", ))
		areaV.units = "m2"
		areaV[:] = self.Delft.NC.elems.area

		pV = self.OP.createVariable("p", "d", ("time", ), fill_value=np.nan)
		pV.long_name = "Precipitation rate"
		pV.units = "mm/s"
		pV[0] = 0
		eV = self.OP.createVariable("et", "d", ("time", ), fill_value=np.nan)
		eV.long_name = "Evapotranspiration rate"
		eV.units = "mm/s"
		eV[0] = 0
		pinV = self.OP.createVariable("influx_p", "d", ("time", ), fill_value=np.nan)
		pinV.long_name = "Precipitation Influx"
		pinV.units = "m3/dt"
		pinV[0] = 0
		einV = self.OP.createVariable("outflux_et", "d", ("time", ), fill_value=np.nan)
		einV.long_name = "Evapotranspiration Outflux"
		einV.units = "m3/dt"
		einV[0] = 0

		gwsV = self.OP.createVariable("gws", "d", ("time", "element"), fill_value=np.nan)
		gwsV.long_name = "Groundwater Levels (Layer 1: Unsaturated) from reference sea level"
		gwsV.units = "m"
		gws_l2V = self.OP.createVariable("gws_l2", "d", ("time", "element"), fill_value=np.nan)
		gws_l2V.long_name = "Groundwater Levels (Layer 2: Aquitard) from reference sea level"
		gws_l2V.units = "m"
		gws_l3V = self.OP.createVariable("gws_l3", "d", ("time", "element"), fill_value=np.nan)
		gws_l3V.long_name = "Groundwater Levels (Layer 3: Quaternary Aquifer) from reference sea level"
		gws_l3V.units = "m"
		gwsFV = self.OP.createVariable("gws_fort", "d", ("time", "element"), fill_value=np.nan)
		gwsFV.long_name = "Groundwater Levels (Layer 1: Unsaturated) from reference sea level: as reported by Fort"
		gwsFV.units = "m"
		swsV = self.OP.createVariable("sws", "d", ("time", "element"), fill_value=np.nan)
		swsV.long_name = "Surfacewater Levels above ground level"
		swsV.units = "m"
		swsFV = self.OP.createVariable("sws_fort", "d", ("time", "element"), fill_value=np.nan)
		swsFV.long_name = "Surfacewater Levels above ground level: as reported by Fort"
		swsFV.units = "m"
		smV = self.OP.createVariable("sm", "d", ("time", "element"), fill_value=np.nan)
		smV.long_name = "Soil Moisture"
		smV.units = "mm"
		epvV = self.OP.createVariable("epv", "d", ("time", "element"), fill_value=np.nan)
		epvV.long_name = "Field Capacity"
		epvV.units = "mm"

		fort_residualV = self.OP.createVariable("fort_residual", "d", ("time", ), fill_value=np.nan)
		fort_residualV.long_name = "Residual from FortRun"
		fort_residualV.units = "mm"
		delft_influxV = self.OP.createVariable("influx_delft_lat", "d", ("time", ), fill_value=np.nan)
		delft_influxV.long_name = "Influx of water into Delft FM: Laterals"
		delft_influxV.units = "m3/dt"
		delft_outfluxV = self.OP.createVariable("outflux_delft_lat", "d", ("time", ), fill_value=np.nan)
		delft_outfluxV.long_name = "Outflux of water into Delft FM: Laterals"
		delft_outfluxV.units = "m3/dt"
		delft_influxV = self.OP.createVariable("influx_delft_bnd", "d", ("time", ), fill_value=np.nan)
		delft_influxV.long_name = "Influx of water into Delft FM: Boundaries"
		delft_influxV.units = "m3/dt"
		delft_outfluxV = self.OP.createVariable("outflux_delft_bnd", "d", ("time", ), fill_value=np.nan)
		delft_outfluxV.long_name = "Outflux of water into Delft FM: Boundaries"
		delft_outfluxV.units = "m3/dt"
		mf6_influxV = self.OP.createVariable("influx_mf6_rch", "d", ("time", ), fill_value=np.nan)
		mf6_influxV.long_name = "Influx of water into Modflow 6: RCH"
		mf6_influxV.units = "m3/dt"
		mf6_outfluxV = self.OP.createVariable("outflux_mf6_rch", "d", ("time", ), fill_value=np.nan)
		mf6_outfluxV.long_name = "Outflux of water from Modflow 6: RCH"
		mf6_outfluxV.units = "m3/dt"
		mf6_influxV = self.OP.createVariable("influx_mf6_ss", "d", ("time", ), fill_value=np.nan)
		mf6_influxV.long_name = "Influx of water into Modflow 6: SS"
		mf6_influxV.units = "m3/dt"
		mf6_outfluxV = self.OP.createVariable("outflux_mf6_ss", "d", ("time", ), fill_value=np.nan)
		mf6_outfluxV.long_name = "Outflux of water from Modflow 6: SS"
		mf6_outfluxV.units = "m3/dt"
		mf6_influxV = self.OP.createVariable("influx_mf6_sy", "d", ("time", ), fill_value=np.nan)
		mf6_influxV.long_name = "Influx of water into Modflow 6: SY"
		mf6_influxV.units = "m3/dt"
		mf6_outfluxV = self.OP.createVariable("outflux_mf6_sy", "d", ("time", ), fill_value=np.nan)
		mf6_outfluxV.long_name = "Outflux of water from Modflow 6: SY"
		mf6_outfluxV.units = "m3/dt"
		mf6_influxV = self.OP.createVariable("influx_mf6_chd", "d", ("time", ), fill_value=np.nan)
		mf6_influxV.long_name = "Influx of water into Modflow 6: CHD"
		mf6_influxV.units = "m3/dt"
		mf6_outfluxV = self.OP.createVariable("outflux_mf6_chd", "d", ("time", ), fill_value=np.nan)
		mf6_outfluxV.long_name = "Outflux of water from Modflow 6: CHD"
		mf6_outfluxV.units = "m3/dt"
		sm_fluxV = self.OP.createVariable("flux_sm", "d", ("time", ), fill_value=np.nan)
		sm_fluxV.long_name = "Flux of Soil Moisture"
		sm_fluxV.units = "m3/dt"


		sw_flux_fortV = self.OP.createVariable("sw_dis_fort", "d", ("time", "element"), fill_value=np.nan)
		sw_flux_fortV.long_name = "SW Storage Discharge as reported by Fort"
		sw_flux_fortV.standard_name = "SW Discharge: Fort"
		sw_flux_fortV.units = "m3/dt"
		sw_flux_fortV[0] = 0
		sw_flux_delftV = self.OP.createVariable("sw_dis_delft", "d", ("time", "element"), fill_value=np.nan)
		sw_flux_delftV.long_name = "SW Storage Discharge as reported by Delft FM"
		sw_flux_delftV.standard_name = "SW Discharge: Delft FM"
		sw_flux_delftV.units = "m3/dt"
		sw_flux_delftV[0] = 0
		gw_flux_fortV = self.OP.createVariable("gw_dis_fort", "d", ("time", "element"), fill_value=np.nan)
		gw_flux_fortV.long_name = "GW Storage Discharge as reported by Fort"
		gw_flux_fortV.standard_name = "GW Discharge: Fort"
		gw_flux_fortV.units = "m3/dt"
		gw_flux_fortV[0] = 0
		gw_flux_mfV = self.OP.createVariable("gw_dis_mf6", "d", ("time", "element"), fill_value=np.nan)
		gw_flux_mfV.long_name = "GW Storage Discharge as reported by Modflow 6"
		gw_flux_mfV.standard_name = "GW Discharge: MF6"
		gw_flux_mfV.units = "m3/dt"
		gw_flux_mfV[0] = 0
		sm_fluxV = self.OP.createVariable("sm_dis", "d", ("time", "element"), fill_value=np.nan)
		sm_fluxV.long_name = "Soil Moisture Discharge"
		sm_fluxV.units = "m3/dt"
		sm_fluxV[0] = 0

		storage_smV = self.OP.createVariable("storage_sm", "d", ("time", ), fill_value=np.nan)
		storage_smV.long_name = "Soil Moisture Storage for the Model Area"
		storage_smV.units = "m3/dt"
		storage_swV = self.OP.createVariable("storage_sw", "d", ("time", ), fill_value=np.nan)
		storage_swV.long_name = "Surface-water Storage for the Model Area"
		storage_swV.units = "m3/dt"
		storage_gwV = self.OP.createVariable("storage_gw", "d", ("time", ), fill_value=np.nan)
		storage_gwV.long_name = "Groundwater Storage (Layer1: Unconfined) for the Model Area"
		storage_gwV.units = "m3/dt"
		gwl2_storageV = self.OP.createVariable("storage_gw_l2", "d", ("time", ), fill_value=np.nan)
		gwl2_storageV.long_name = "Groundwater Storage (Layer2: Aquitard) for the Model Area"
		gwl2_storageV.units = "m3/dt"
		gwl3_storageV = self.OP.createVariable("storage_gw_l3", "d", ("time", ), fill_value=np.nan)
		gwl3_storageV.long_name = "Groundwater Storage (Layer3: Quaternary) for the Model Area"
		gwl3_storageV.units = "m3/dt"

		obs_g2V = self.OP.createVariable("obs_g2", "d", ("time"), fill_value=np.nan)
		obs_g2V.long_name = "Water Level at Observation Station: graben_002_corr, i.e. Upstream"
		obs_g2V.standard_name = "graben_002_corr"
		obs_g2V.units = "m"
		obs_g1V = self.OP.createVariable("obs_g1", "d", ("time"), fill_value=np.nan)
		obs_g1V.long_name = "Water Level at Observation Station: graben_001_corr, i.e. Downstream"
		obs_g1V.standard_name = "graben_001_corr"
		obs_g1V.units = "m"
		obs_p7 = self.OP.createVariable("obs_p7", "d", ("time"), fill_value=np.nan)
		obs_p7.long_name = "Groundwater Level at Observation Station P7"
		obs_p7.standard_name = "Pegel 7"
		obs_p7.units = "m"
		obs_p6 = self.OP.createVariable("obs_p6", "d", ("time"), fill_value=np.nan)
		obs_p6.long_name = "Groundwater Level at Observation Station P6"
		obs_p6.standard_name = "Pegel 6"
		obs_p6.units = "m"
		obs_pG = self.OP.createVariable("obs_pG", "d", ("time"), fill_value=np.nan)
		obs_pG.long_name = "Groundwater Level at Observation Station PGraben"
		obs_pG.standard_name = "Pegel am Graben"
		obs_pG.units = "m"
		obs_p12 = self.OP.createVariable("obs_p12", "d", ("time"), fill_value=np.nan)
		obs_p12.long_name = "Groundwater Level at Observation Station P12"
		obs_p12.standard_name = "Pegel 12"
		obs_p12.units = "m"
		obs_p2 = self.OP.createVariable("obs_p2", "d", ("time"), fill_value=np.nan)
		obs_p2.long_name = "Groundwater Level at Observation Station P2"
		obs_p2.standard_name = "Pegel 2"
		obs_p2.units = "m"
		obs_p1 = self.OP.createVariable("obs_p1", "d", ("time"), fill_value=np.nan)
		obs_p1.long_name = "Groundwater Level at Observation Station P1"
		obs_p1.standard_name = "Pegel 1"
		obs_p1.units = "m"
		obs_p2_sm = self.OP.createVariable("obs_p2_sm", "d", ("time"), fill_value=np.nan)
		obs_p2_sm.long_name = "Soil Moisture Storage at Observation Station P2"
		obs_p2_sm.standard_name = "Pegel 2: SM"
		obs_p2_sm.units = "m"

		timeV[0] = 0
		time_unixV[0] = self.times.unix["global_start"]
		gwsV[0] = self.MF6.Res.h[0,-1,:]
		gws_l2V[0] = self.MF6.Res.h[1,-1,:]
		gws_l3V[0] = self.MF6.Res.h[2,-1,:]
		gwsFV[0] = gwsV[0]
		swsV[0] = self.Delft.Res.sws[0,:]
		swsFV[0] = self.Delft.Res.sws[0,:]
		smV[0] = self.Fort.Ini.sm
		epvV[0] = self.Fort.Ini.epv

		storage_smV[0] = np.multiply(self.Fort.Ini.sm[:], self.Delft.NC.elems.area.T).sum()/1000
		storage_swV[0] = np.multiply(self.Fort.Ini.sws[:], self.Delft.NC.elems.area.T).sum()/1000
		storage_gwV[0] = np.multiply(self.MF6.Res.h[0,0,:]-self.MF6.disv.botm.data[0], self.Delft.NC.elems.area.T).sum()\
		*self.MF6.sto.sy.data[0].mean()
		gwl2_storageV[0] = np.multiply(self.MF6.Res.h[1,0,:]-self.MF6.disv.botm.data[1], self.MF6.Params.aq_volume[1].T).sum()*\
		self.MF6.sto.ss.data[1].mean()
		gwl3_storageV[0] = np.multiply(self.MF6.Res.h[2,0,:]-self.MF6.disv.botm.data[2], self.MF6.Params.aq_volume[2].T).sum()*\
		self.MF6.sto.ss.data[2].mean()
		logging.info("RES: NC file initialized")

	def update(self, save=False):
		self.save = save
		idx = self.times.nTS["run_num"] + 1
		self.OP["nTS"][idx] = self.times.nTS["exchange"]
		self.OP["time"][idx] = self.times.rel["local_end"]
		self.OP["time_unix"][idx] = self.times.unix["local_end"]
		self.OP["gws"][idx] = self.MF6.Res.h[0,-1,:]
		self.OP["gws_fort"][idx] = self.Fort.Res.gws[:,-1].T/1000
		self.OP["gws_l2"][idx] = self.MF6.Res.h[1,-1,:]
		self.OP["gws_l3"][idx] = self.MF6.Res.h[2,-1,:]
		self.OP["sws"][idx] = self.Delft.Res.sws[1:,:]
		self.OP["sws_fort"][idx] = self.Fort.Res.sws[:,-1].T/1000
		self.OP["sm"][idx] = self.Fort.Res.sm[:,-1].T
		self.OP["epv"][idx] = self.Fort.Res.epv[:,-1].T
		self.OP["fort_residual"][idx] = np.abs(self.Fort.Res.Qdiff[1:]).sum()
		self.OP["p"][idx] = self.Fort.Ini.p[1:].sum()
		self.OP["et"][idx] = self.Fort.Ini.et[1:].sum()

		self.OP["sw_dis_fort"][idx] = self.Delft.sw_flux[:,1:].sum(axis=1)*self.times.dt["FortRun"]
		mapNC = nc.Dataset(self.Delft.mapNC_path)
		self.OP["sw_dis_delft"][idx] = mapNC["mesh2d_qin"][:].sum(axis=0)*self.times.dt["FortRun"]
		mapNC.close()
		self.gw_flux = np.zeros(self.Fort.Res.gw_dis.shape)
		for i in range(self.gw_flux.shape[1]):
			self.gw_flux[:,i] = np.multiply(self.Fort.Res.gw_dis[:,i], self.Delft.NC.elems.area.T).T
		self.OP["gw_dis_fort"][idx] = (self.gw_flux[:,1:].sum(axis=1)*self.times.dt["FortRun"])
		rch_cbb = []
		mf6_dis = np.zeros(self.Fort.Res.gw_dis.T[1:].shape)
		for t in range(self.times.nTS["exchange"]):
			rch_cbb.append(self.MF6.cbb.get_data(text="RCH")[t])
			for el in range(rch_cbb[t].size):
				mf6_dis[t,el] = rch_cbb[t][el][2]
		self.OP["gw_dis_mf6"][idx] = mf6_dis.sum(axis=0)*self.times.dt["FortRun"]
		self.OP["influx_mf6_rch"][idx] = self.MF6.lst.get_cumulative(names="RCH_IN")[t][3]
		self.OP["outflux_mf6_rch"][idx] = self.MF6.lst.get_cumulative(names="RCH_OUT")[-1][3]

		self.OP["influx_mf6_ss"][idx] = self.MF6.lst.get_cumulative(names="STO-SS_IN")[-1][3]
		self.OP["outflux_mf6_ss"][idx] = self.MF6.lst.get_cumulative(names="STO-SS_OUT")[-1][3]
		self.OP["influx_mf6_sy"][idx] = self.MF6.lst.get_cumulative(names="STO-SY_IN")[-1][3]
		self.OP["outflux_mf6_sy"][idx] = self.MF6.lst.get_cumulative(names="STO-SY_OUT")[-1][3]
		self.OP["influx_mf6_chd"][idx] =self.MF6.lst.get_cumulative(names="CHD_IN")[-1][3]
		self.OP["outflux_mf6_chd"][idx] = self.MF6.lst.get_cumulative(names="CHD_OUT")[-1][3]

		self.sm_flux = np.zeros(self.Fort.Res.sm_dis.shape)
		for i in range(self.sm_flux.shape[1]):
			self.sm_flux[:,i] = np.multiply(self.Fort.Res.sm_dis[:,i], self.Delft.NC.elems.area.T).T
		self.sm_flux = self.sm_flux
		self.OP["sm_dis"][idx] = (self.sm_flux[:,1:].sum(axis=1)*self.times.dt["FortRun"])
		self.OP["flux_sm"][idx] = (self.sm_flux[:,1:].sum()*self.times.dt["FortRun"])

		self.OP["storage_sm"][idx] = np.multiply(self.Fort.Res.sm[:,-1], self.Delft.NC.elems.area.T).sum()/1000
		self.OP["storage_sw"][idx] = np.multiply(self.Delft.Res.sws[-1,:], self.Delft.NC.elems.area.T).sum()
		self.OP["storage_gw"][idx] = np.multiply(self.MF6.Res.h[0,-1,:]-self.MF6.disv.botm.data[0], self.Delft.NC.elems.area.T).sum()\
		*self.MF6.sto.sy.data[0].mean()
		self.OP["storage_gw_l2"][idx] = np.multiply(self.MF6.Res.h[1,-1,:]-self.MF6.disv.botm.data[1], self.MF6.Params.aq_volume[1].T).sum()\
		*self.MF6.sto.ss.data[1].mean()
		self.OP["storage_gw_l3"][idx] = np.multiply(self.MF6.Res.h[2,-1,:]-self.MF6.disv.botm.data[2], self.MF6.Params.aq_volume[2].T).sum()\
		*self.MF6.sto.ss.data[2].mean()

		self.OP["influx_p"][idx] = self.OP["p"][idx]/1000*self.Delft.NC.elems.area.sum()*self.times.dt["FortRun"]
		self.OP["outflux_et"][idx] = self.OP["et"][idx]/1000*self.Delft.NC.elems.area.sum()*self.times.dt["FortRun"]

		hisNC = nc.Dataset(self.Delft.hisNC_path)
		self.OP["obs_g1"][idx] = hisNC["waterlevel"][-1,1]
		self.OP["obs_g2"][idx] = hisNC["waterlevel"][-1,0]
		self.OP["influx_delft_lat"][idx] = hisNC["water_balance_laterals_in"][-1]
		self.OP["outflux_delft_lat"][idx] = hisNC["water_balance_laterals_out"][-1]
		self.OP["influx_delft_bnd"][idx] = hisNC["water_balance_boundaries_in"][-1]
		self.OP["outflux_delft_bnd"][idx] = hisNC["water_balance_boundaries_out"][-1]

		hisNC.close()
		mf_obs = pd.read_csv(self.MF6.path+"0")
		self.OP["obs_p7"][idx] = mf_obs["P7"].to_numpy()[-1]
		self.OP["obs_p6"][idx] = mf_obs["P6"].to_numpy()[-1]
		self.OP["obs_pG"][idx] = mf_obs["PGRABEN"].to_numpy()[-1]
		self.OP["obs_p12"][idx] = mf_obs["P12"].to_numpy()[-1]
		self.OP["obs_p2"][idx] = mf_obs["P2"].to_numpy()[-1]
		self.OP["obs_p1"][idx] = mf_obs["P1"].to_numpy()[-1]
		self.OP["obs_p2_sm"][idx] = self.Fort.Res.sm[2182,-1].T

		logging.info("*RES: Wrote OP to NC file\n")
		if save:
			self.saveOPs()
		else:
			self.times.runtimes.append(time.time() - self.times.run["start"])
			self.OP["runtime"][idx] = self.times.runtimes[-1]
			logging.info("*RUNTIME{}: {}\n\n".format(self.times.nTS["run_num"], self.times.runtimes[int(self.times.nTS["run_num"]-1)]))
		self.times.nTS["ran"] = self.times.nTS["ran"] + self.times.nTS["exchange"]
		self.times.nTS["run_num"] = self.times.nTS["run_num"] + 1

	def saveOPs(self):
		if not os.path.exists(self.delft_path):
			os.mkdir(self.delft_path)
		his_file = os.path.join(self.delft_path, self.Delft.model_name+"_"+self.times.stmp["local_start"]\
		.strftime("%Y%m%d_%H%M%S")+"_his.nc")
		map_file = os.path.join(self.delft_path, self.Delft.model_name+"_"+self.times.stmp["local_start"]\
		.strftime("%Y%m%d_%H%M%S")+"_map.nc")
		shutil.copy(self.Delft.mapNC_path, map_file)
		shutil.copy(self.Delft.hisNC_path, his_file)
		logging.info("*RES: Copied DELFT outputs to output folder")
		if not os.path.exists(self.mf_path):
			os.mkdir(self.mf_path)
		hds_file = os.path.join(self.mf_path, self.MF6.name+"_"+self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S")+".hds")
		cbb_file = os.path.join(self.mf_path, self.MF6.name+"_"+self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S")+".cbb")
		lst_file = os.path.join(self.mf_path, self.MF6.name+"_"+self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S")+".lst")
		shutil.copy(os.path.join(self.MF6.path, self.MF6.name+".lst"), lst_file)
		shutil.copy(os.path.join(self.MF6.path, self.MF6.name+".cbb"), cbb_file)
		shutil.copy(os.path.join(self.MF6.path, self.MF6.name+".hds"), hds_file)
		logging.info("*RES: Copied MF outputs to output folder")
		if not os.path.exists(self.fort_path):
			os.mkdir(self.fort_path)
		fort_file = os.path.join(self.fort_path, "wasenmoos_fort_"+self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S"))
		np.savez(fort_file, gws=self.Fort.Res.gws, sws=self.Fort.Res.sws, sm=self.Fort.Res.sm, epv=self.Fort.Res.epv,\
		Qin=self.Fort.Res.Qin, Qout=self.Fort.Res.Qout, Qdiff=self.Fort.Res.Qdiff,\
		gw_dis=self.Fort.Res.gw_dis, sw_dis=self.Fort.Res.sw_dis, sm_dis=self.Fort.Res.sm_dis)
		logging.info("*RES: Dumped FORT outputs to output folder")
		self.times.runtimes.append(time.time() - self.times.run["start"])
		logging.info("*RUNTIME{}: {}\n\n".format(self.times.nTS["run_num"], self.times.runtimes[int(self.times.nTS["run_num"]-1)]))
		idx = self.times.nTS["run_num"] + 1
		self.OP["runtime"][idx] = self.times.runtimes[-1]

		nxt_rst_file = self.Delft.model_name+"_"+self.times.stmp["local_end"].strftime("%Y%m%d_%H%M%S")+"_rst.nc"
		if not os.path.exists(os.path.join(self.Delft.OP_path_abs, nxt_rst_file)):
			self.OP.close()
			logging.error("Restart file for next timestep not generated. Check DIA/MDU files.\n\n\n*EXITING*")
			raise RuntimeError("Restart file for next timestep not generated. Check DIA/MDU files.")

	def close(self):
		self.OP.close()
		if self.save:
			rst_file = self.Delft.model_name+"_"+self.times.stmp["local_start"].strftime("%Y%m%d_%H%M%S")+"_rst.nc"
			rst_Cpath = os.path.join(self.delft_path, rst_file)
			rst_Opath = os.path.join(self.Delft.OP_path_abs, rst_file)
			shutil.copy(rst_Opath, rst_Cpath)
			logging.info("*RES: Copied the final DELFT rst file")
		logging.info("*RES: Closed NC file\nMaking a copy of the log")
		shutil.copy(os.path.join(self.path,"GWSWEX.log"), os.path.join(self.path,"GWSWEX_"+self.times.stmp["global_start"].\
		strftime("%Y%m%d_%H%M%S")+"-"+self.times.stmp["local_end"].strftime("%Y%m%d_%H%M%S")+".log"))
		logging.info("*LOG END*\n\n\n")
		logging.shutdown()

	def loadNC(self, nc_path=None):
		if not nc_path:
			nc_path = os.path.join(self.path, self.name+".nc")
		self.OP = nc.Dataset(nc_path)
		return self.OP

	def report(self, drop=None, plot=False, closeNC=True):
		self.loadNC()
		if drop:
			n = self.OP["time_unix"][:].size - drop
		else:
			n = self.OP["time_unix"][:].size
		self.rep = bundleIT()
		self.rep.p = self.OP["influx_p"][:n].sum()
		self.rep.et = self.OP["outflux_et"][:n].sum()
		self.rep.sm = self.OP["flux_sm"][:n].sum()
		self.rep.mf6_in, self.rep.mf6_out = {}, {}
		self.rep.mf6_in["rch"]= self.OP["influx_mf6_rch"][:n].sum()
		self.rep.mf6_in["ss"] = self.OP["influx_mf6_ss"][:n].sum()
		self.rep.mf6_in["sy"] = self.OP["influx_mf6_sy"][:n].sum()
		self.rep.mf6_in["chd"] = self.OP["influx_mf6_chd"][:n].sum()
		self.rep.mf6_in["total"] = sum(self.rep.mf6_in.values())
		self.rep.mf6_out["rch" ] = self.OP["outflux_mf6_rch"][:n].sum()
		self.rep.mf6_out["ss" ]= self.OP["outflux_mf6_ss"][:n].sum()
		self.rep.mf6_out["sy"] = self.OP["outflux_mf6_sy"][:n].sum()
		self.rep.mf6_out["chd"] = self.OP["outflux_mf6_chd"][:n].sum()
		self.rep.mf6_out["total"] = sum(self.rep.mf6_out.values())
		self.rep.delft_in, self.rep.delft_out= {}, {}
		self.rep.delft_in["lat"] = self.OP["influx_delft_lat"][:n].sum()
		self.rep.delft_in["bnd"] = self.OP["influx_delft_bnd"][:n].sum()
		self.rep.delft_in["total"] = sum(self.rep.delft_in.values())
		self.rep.delft_out["lat"] = self.OP["outflux_delft_lat"][:n].sum()
		self.rep.delft_out["bnd"] = self.OP["outflux_delft_bnd"][:n].sum()
		self.rep.delft_out["total"] = sum(self.rep.delft_out.values())
		self.rep.err = {}
		self.rep.err["fort"] = np.abs(self.OP["fort_residual"][:n]).sum()*self.Delft.NC.elems.area.sum()/10000
		self.rep.fluxes = {}
		self.rep.fluxes["sw_fort"] = self.OP["sw_dis_fort"][:n].sum()
		self.rep.fluxes["sw_delft"] = self.OP["sw_dis_delft"][:n].sum()
		self.rep.fluxes["gw_fort"] = self.OP["gw_dis_fort"][:n].sum()
		self.rep.fluxes["gw_mf6"] = self.OP["gw_dis_mf6"][:n].sum()
		self.rep.err["fort-delft"] = (abs(self.OP["sw_dis_fort"][:n]).sum() - (abs(self.OP["sw_dis_delft"][:n])).sum())
		self.rep.err["fort-mf6"] = (self.rep.fluxes["gw_fort"] - self.rep.fluxes["gw_mf6"])
		self.rep.storage_delta = {}
		self.rep.storage_delta["sm"] = (self.OP["storage_sm"][-1] - self.OP["storage_sm"][0])
		self.rep.storage_delta["sw"] = (self.OP["storage_sw"][-1] - self.OP["storage_sw"][0])
		self.rep.storage_delta["gw"] = (self.OP["storage_gw"][-1] - self.OP["storage_gw"][0])
		self.rep.storage_delta["gw_l2"] = (self.OP["storage_gw_l2"][-1] - self.OP["storage_gw_l2"][0])
		self.rep.storage_delta["gw_l3"] = (self.OP["storage_gw_l2"][-1] - self.OP["storage_gw_l3"][0])
		self.rep.storage_delta["total_ini"] = (self.OP["storage_sm"][0] + self.OP["storage_sw"][0] + self.OP["storage_gw"][0] + \
		self.OP["storage_gw_l2"][0] + self.OP["storage_gw_l3"][0]).data
		self.rep.storage_delta["total_fin"] = (self.OP["storage_sm"][-1] + self.OP["storage_sw"][-1] + self.OP["storage_gw"][-1] +\
		self.OP["storage_gw_l2"][-1] + self.OP["storage_gw_l3"][-1]).data
		self.rep.storage_delta["diff"] = self.rep.storage_delta["total_fin"] - self.rep.storage_delta["total_ini"]
		self.rep.Qin = self.rep.p - self.rep.et
		self.rep.Qout = self.rep.sm + self.rep.mf6_in["rch"]-self.rep.mf6_out["rch"] + self.rep.delft_in["lat"]-self.rep.delft_out["lat"]
		Qdiff_abs = self.rep.Qin - self.rep.Qout
		Qdiff_rel = Qdiff_abs/self.rep.Qin
		Qdiff_perc = Qdiff_rel*100
		self.rep.err["abs"] = Qdiff_abs
		self.rep.err["rel"] = Qdiff_rel
		self.rep.err["perc"] = Qdiff_perc
		print("Absolute Error:", Qdiff_abs, "m3\nRelative Error:", Qdiff_rel,"\nPercentage Error: %.2f" % Qdiff_perc, "\b%")
		if plot:
			if plot == "delft":
				fig = plt.figure()
				fig.suptitle("Sum of residuals over time steps of FortRun")
				plt.plot(self.OP["sw_dis_fort"][:].sum(axis=1), label="fort")
				plt.plot(self.OP["sw_dis_delft"][:].sum(axis=1), label="delft")
				plt.legend()
			if plot == "mf6":
				fig = plt.figure()
				fig.suptitle("Sum of residuals over time steps of FortRun")
				plt.plot(self.OP["gw_dis_fort"][:].sum(axis=1), label="fort")
				plt.plot(self.OP["gw_dis_mf6"][:].sum(axis=1), label="mf6")
				plt.legend()
			if plot == "fort":
				fig = plt.figure()
				fig.suptitle("Sum of residuals over time steps of FortRun")
				plt.plot(self.OP["fort_residual"][:], label="fort")
				plt.legend()
		return self.rep

	def loadObs(self, path=None):
		self.Obs = {}
		if not path:
			self.obs_path = "../data/obs/"
		#
		def readOBS(file):
			DF = pd.read_table(os.path.join(self.obs_path, file))
			DF[DF.columns[0]] = pd.to_datetime(DF[DF.columns[0]])
			DF.index = DF[DF.columns[0]]
			DF = DF.drop(DF.columns[0], 1)
			start = datetime(self.times.stmp["global_start"].year, self.times.stmp["global_start"].month, self.times.stmp["global_start"].day,\
			self.times.stmp["global_start"].hour)
			end = datetime(self.times.stmp["global_end"].year, self.times.stmp["global_end"].month, self.times.stmp["global_end"].day,\
			self.times.stmp["global_end"].hour)
			p7_mask = (DF.index >= start) & (DF.index <= end)
			DF = DF.loc[p7_mask]
			return DF, DF[DF.columns[0]].to_numpy()
		#
		self.Obs["p7DF"], self.Obs["p7"] = readOBS("P7_20200501-0601.dat")
		self.Obs["p6DF"], self.Obs["p6"] = readOBS("P6_20200501-0601.dat")
		self.Obs["pGDF"], self.Obs["pG"] = readOBS("PG_20200501-0601.dat")
		self.Obs["p12DF"], self.Obs["p12"] = readOBS("P12_20200501-0601.dat")
		self.Obs["p2DF"], self.Obs["p2"] = readOBS("P2_20200501-0601.dat")
		self.Obs["p1DF"], self.Obs["p1"] = readOBS("P1_20200501-0601.dat")
		self.Obs["p2_smDF"], self.Obs["p2_sm"] = readOBS("P2_SM_20200501-0601.dat")

	def loadFort(self, load_time=None):
		if not load_time:
			load_time = (self.times.stmp["local_start"] - timedelta(seconds=self.times.dt["exchange"])).strftime("%Y%m%d_%H%M%S")
		fort_path = os.path.join(self.fort_path, "wasenmoos_fort_"+load_time+".npz")
		fort_file = np.load(fort_path)
		self.Fort.Res.gws = fort_file["gws"]
		self.Fort.Res.gw_dis = fort_file["gw_dis"]
		self.Fort.Res.sws = fort_file["sws"]
		self.Fort.Res.sw_dis = fort_file["sw_dis"]
		self.Fort.Res.sm = fort_file["sm"]
		self.Fort.Res.sm_dis = fort_file["sm_dis"]
		self.Fort.Res.epv = fort_file["epv"]
		self.Fort.Res.Qin = fort_file["Qin"]
		self.Fort.Res.Qout = fort_file["Qout"]
		self.Fort.Res.Qdiff = fort_file["Qdiff"]

	def loadMF6(self, load_time=None):
		if not load_time:
			load_time = (self.times.stmp["local_start"] - timedelta(seconds=self.times.dt["exchange"])).strftime("%Y%m%d_%H%M%S")
		head_file = self.MF6.name+"_"+load_time+".hds"
		self.MF6.hds = fp.utils.binaryfile.HeadFile(os.path.join(self.mf_path, head_file))
		self.MF6.Res.h = self.MF6.hds.get_data(kstpkper=(self.MF6.Params.internalTS-1, len(self.MF6.hds.get_kstpkper())-1))
		cbb_file = self.MF6.name+"_"+load_time+".cbb"
		self.MF6.cbb = fp.utils.binaryfile.CellBudgetFile(os.path.join(self.mf_path, cbb_file))
		lst_file = self.MF6.name+"_"+load_time+".lst"
		self.MF6.lst = fp.utils.mflistfile.Mf6ListBudget(os.path.join(self.mf_path, lst_file))

	def loadDelft(self, load_time=None):
		if not load_time:
			load_time = map_file = (self.times.stmp["local_start"]-timedelta(seconds=self.times.dt["exchange"])\
			).strftime("%Y%m%d_%H%M%S")
		map_file = self.delft_path + self.Delft.model_name + "_" + load_time + "_map.nc"
		mapNC = nc.Dataset(map_file)
		self.Delft.Res.sws_abs = np.array(mapNC["mesh2d_s1"][:])
		self.Delft.Res.sws = np.array(mapNC["mesh2d_waterdepth"][:])
		self.Delft.Res.sws_prev = np.array(mapNC["mesh2d_s0"][:])
		mapNC.close()

	def plotMesh(self, element_wise_data, fromNC=None, ts=-1, title=None, item=None, save=False, pt=None):
		fig = plt.figure()
		fig.tight_layout()
		pal = {"gw":"#E74C3C", "sm":"#2ECC71", "fc":"#5EFCA1", "sv":"#E7AC3C", "sw":"#2980B9", "p":"#1A3BFF", "et":"#FF6600"\
		, "red":"red", "green":"green", "grey":"grey", "cyan":"c", "pink":"m", "jet": mpl.cm.jet}
		if item in pal and item != "jet":
			col = mpl.colors.LinearSegmentedColormap.from_list("", ["#e1e1e1", pal[item]])
		elif item == "jet":
			col = mpl.cm.jet
		else:
			col = mpl.cm.turbo
		if title:
			fig.suptitle(title)
		if fromNC == "key":
			element_wise_data = self.OP[element_wise_data][:]
		elif fromNC == "sum":
			element_wise_data = self.OP[element_wise_data][:].sum(axis=0)
		elif fromNC == "mean":
			element_wise_data = self.OP[element_wise_data][:].mean(axis=0)
		elif fromNC == "ts":
			element_wise_data = self.OP[element_wise_data][ts]
		pmv = fp.plot.PlotMapView(model=self.MF6.model)
		pmv.plot_grid(linewidth=0.05)
		hplot = pmv.plot_array(element_wise_data, cmap=col)
		plt.colorbar(hplot)
		plt.tight_layout()
		if save:
			if pt:
				title = pt
			plt.savefig(os.path.join(self.plt_path, title+".svg"), format="svg", dpi=600, bbox_inches='tight')

	def plot1D(self, data, fromNC=None, items=None, labels=None, title=None, to_idx=None, save=False, dual=True, pt=None):
		pal = {"gw":"#E74C3C", "sm":"#2ECC71", "fc":"#5EFCA1", "sv":"#E7AC3C", "sw":"#2980B9", "p":"#1A3BFF", "et":"#FF6600"\
		, "red":"red", "green":"green", "grey":"grey", "cyan":"c", "pink":"m", "default":"black", "d":"black"}
		if not to_idx and type(data) != str:
			to_idx = min(data[0].size, self.OP["time_unix"][:].size)
		else:
			to_idx = self.OP["time_unix"][:].size
		fig, ax = plt.subplots()
		manager = plt.get_current_fig_manager()
		manager.window.showMaximized()
		fig.autofmt_xdate()
		ax.ticklabel_format(useOffset=False)
		ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%d-%m-%y'))
		plt.tight_layout()
		if title:
			fig.suptitle(title)
		c = 0
		posix_times = self.OP["time_unix"][:]
		x = [datetime.fromtimestamp(x) for x in posix_times]
		if items and not labels:
			labels = items
		if not items:
			items = ["default", "default", "default", "default", "default", "default"]
			labels = [data, ]
		if fromNC == "key":
			data = [self.OP[data][:], ]
		elif fromNC == "sum":
			data = [self.OP[data][:].sum(axis=1), ]
		elif fromNC == "mean":
			data = [self.OP[data][:].mean(axis=1), ]
		if len(data) == 2 and dual:
			for lin in data:
				if c != 0:
					ax=ax.twinx()
				ax.plot(x[:to_idx], lin[:to_idx], color=pal[items[c]], linewidth=0.75)
				ax.set_ylabel(labels[c], color=pal[items[c]])
				ax.tick_params(axis='y', labelcolor=pal[items[c]])
				plt.legend()
				c += 1
		else:
			for lin in data:
				plt.plot(x[:to_idx], lin[:to_idx], color=pal[items[c]], label=labels[c], linewidth=0.75)
				plt.legend()
				plt.tight_layout()
				c += 1
		if save:
			if not title and not pt:
				title = ""
				if items:
					for x in items:
						title = title + x
				elif labels:
					for x in items:
						title = title + x
				else:
					title = "plot"
			if pt:
				title = pt
			fig.savefig(os.path.join(self.plt_path, title+".svg"), format="svg", dpi=600, bbox_inches='tight')