#!/bin/bash

while getopts ":qvpol" opt; do
    case $opt in
        q)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature
            f2py -c GWSWEX.pyf GWSWEX.f90 --quiet
            mv *.so ~/GWSWEX/
            cd ~/GWSWEX
            ;;
        v)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature
            f2py -c GWSWEX.pyf GWSWEX.f90
            mv *.so ~/GWSWEX/
            cd ~/GWSWEX
            ;;
        p)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature
            f2py -c GWSWEX.pyf GWSWEX.f90 --f90flags="-fopenmp" -lgomp --quiet
            mv *.so ~/GWSWEX/
            cd ~/GWSWEX
            ;;
        o)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature
            f2py -c GWSWEX.pyf GWSWEX.f90 --f90flags="-fopenmp" -lgomp
            mv *.so ~/GWSWEX/
            cd ~/GWSWEX
            ;;
        l)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature
            f2py -c -L/usr/local/lib/ GWSWEX.pyf GWSWEX.f90 --f90flags="-fopenmp" -lgomp -lfgsl #--quiet -I/usr/local/include/fgsl/
            mv *.so ~/GWSWEX/
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
    cd src/fort/
    f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature
    f2py -c GWSWEX.pyf GWSWEX.f90 --quiet
    cd ~/GWSWEX
fi