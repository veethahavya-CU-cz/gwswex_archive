#!/bin/bash
cd src/fort
clear
gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX && echo "successfully compiled GWSWEX"
mv GWSWEX ../../exe/fort