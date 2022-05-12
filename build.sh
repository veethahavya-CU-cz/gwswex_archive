#!/bin/bash
cd src/fort
clear
gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX && echo "successfully compiled GWSWEX"
mv GWSWEX ../../exe/fort



while getopts ":cpn" opt; do
    case $opt in
        c)
            cd src/fort/
            gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -o GWSWEX && echo "successfully compiled GWSWEX"
            mv GWSWEX ../../exe/fort
            cd ~/GWSWEX
            ;;
        p)
            cd src/fort/
            gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX && clear && echo "successfully compiled GWSWEX"
            mv GWSWEX ../../exe/fort
            cd ~/GWSWEX
            ;;
        n)
            cd src/fort/
            gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX && echo "successfully compiled GWSWEX"
            mv GWSWEX ../../exe/fort
            cd ~/GWSWEX
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        *)
            echo "Error in command line parsing" >&2
            exit 1
            ;;
    esac
done

if (( $OPTIND == 1 )); then
    cd src/fort/ && clear
    gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX && echo "successfully compiled GWSWEX"
    mv GWSWEX ../../exe/fort
    cd ~/GWSWEX
fi