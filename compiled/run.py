#%%
import GWSWEX
import numpy as np

Fort = GWSWEX.Fort("../fortran/", None)

Fort.Ini.elems = 1
Fort.Ini.ts = 300
Fort.Ini.ts_size = 60
Fort.Ini.k = np.full(Fort.Ini.elems, 0)
Fort.Ini.n = 0.3
Fort.Ini.n_gw = 0.3
Fort.Ini.m = 0.3
Fort.Ini.beta = 0.7
Fort.Ini.alpha = 0.15
Fort.Ini.sw_th = 5e-2
Fort.Ini.chd_cells = np.full(Fort.Ini.elems, 0)
Fort.Ini.gok = np.random.default_rng().uniform(-3, 3, Fort.Ini.elems) + 100
Fort.Ini.gws = np.random.default_rng().uniform(-1, 1, Fort.Ini.elems) + 80
Fort.Ini.epv = (Fort.Ini.gok - Fort.Ini.gws)*Fort.Ini.n
Fort.Ini.sm = Fort.Ini.epv*0.8

Fort.build(restart=False, res=None, run=True)

# %%
