# GWSWEX standalone variant

# Dependencies
* python 3: numpy, scipy
* fortran: gfortran, OMP, fgsl

# Installation
Run build.sh (tweak the fgsl lib and include locations as applicable) to compile the standalone program.   
Build Flags:  
-p: with OMP support (parallelization)
-c: quiet, no OMP
-n: regular verbosity, with OMP
