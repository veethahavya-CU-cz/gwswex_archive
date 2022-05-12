#!/bin/bash
cd src/fort
clear
gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -o GWSWEX && echo "successfully compiled GWSWEX"
mv GWSWEX ../../exe/fort || cd ../.. && mkdir -p exe/fort/ && clear && mv src/fort/GWSWEX exe/fort/ && echo "successfully compiled GWSWEX"