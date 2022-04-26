#!/bin/bash

while getopts ":qv" opt; do
    case $opt in
        q)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature && f2py -c GWSWEX.pyf GWSWEX.f90 --quiet
            mv *.so ~/GWSWEX/
            cd ~/GWSWEX
            ;;
        v)
            cd src/fort/
            f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature && f2py -c GWSWEX.pyf GWSWEX.f90
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
    f2py GWSWEX.f90 -m fortrapper -h GWSWEX.pyf --overwrite-signature && f2py -c GWSWEX.pyf GWSWEX.f90 --quiet
    mv *.so ~/GWSWEX/
    cd ~/GWSWEX
fi