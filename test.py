import os, sys, psutil
sys.path.append('libs/')
from gwswex_wrapper import gwswex

os.environ['OMP_NUM_THREADS'] = str(psutil.cpu_count(logical = False))

gwswex.initialize()

gwswex.finalize()