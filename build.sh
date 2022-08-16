#!/bin/bash

function help() {
    echo "Usage: ./build.sh [-p] [-n]"
    echo " -p: compile with openMP"
    echo " -n: compile without openMP"
    echo " -h: print this help"
    exit 0
}

while getopts ":pnh" opt; do
    case $opt in
        p)
            echo "compiling with openMP"
            cd src/fort/
            gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX
            if [ $? -eq 0 ]; then
                echo "successfully compiled GWSWEX"
                (rm -rf ../../exe/fort/* && mv GWSWEX ../../exe/fort) || (mkdir -p ../../exe/fort && mv GWSWEX ../../exe/fort && echo "dir created")
                echo "placed the program in exe/fort"
                cd ../../
            else
                echo "GWSWEX compilation failed"
            fi
            ;;
        n)
            echo "compiling without openMP"
            cd src/fort/
            gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -o GWSWEX
            if [ $? -eq 0 ]; then
                echo "successfully compiled GWSWEX"
                (rm -rf ../../exe/fort/* && mv GWSWEX ../../exe/fort) || (mkdir -p ../../exe/fort && mv GWSWEX ../../exe/fort && echo "dir created")
                echo "placed the program in exe/fort"
                cd ../../
            else
                echo "GWSWEX compilation failed"
            fi
            ;;
        h)
            help
            ;;
        \?)
            echo "invalid option: -$OPTARG" >&2
            exit 1
            ;;
        *)
            echo "error in command line parsing" >&2
            exit 1
            ;;
    esac
done

# Defaults to compiling with openMP
if (( $OPTIND == 1 )); then
	cd src/fort/
	gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -o GWSWEX
	if [ $? -eq 0 ]; then
		echo "successfully compiled GWSWEX"
		mv GWSWEX ../../exe/fort || mkdir -p ../../exe/fort && mv GWSWEX ../../exe/fort && echo "dir created"
		echo "placed the program in exe/fort"
		cd ~/GWSWEX
	else
		echo "GWSWEX compilation failed"
	fi
fi