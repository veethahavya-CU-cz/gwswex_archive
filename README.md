# GWSWEX f2py variant

# Dependencies
* python 3: numpy, scipy
* fortran: gfortran, OMP

# Installation
Run build.sh to recompile the shared object library file (to import into py) if using a different OS or py version.  
The included binary is compiled for: Debian GNU/Linux 11 (bullseye) with GNU Fortran (Debian 10.2.1-6) 10.2.1 20210110, and Python 3.9.2.  
Build Flags:  
-q: quiet  
-v: verbose  
-p: with OMP support + -q  
-o: with OMP support + -v  
