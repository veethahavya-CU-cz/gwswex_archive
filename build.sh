#!/bin/bash

function help() {
    echo 'Usage: ./build.sh [-p] [-n]'
    echo ' -p: compile with OpenMP'
    echo ' -n: compile without OpenMP'
    echo ' -h: print this help'
    exit 0
}

while getopts ':pnh' opt; do
    case $opt in
        p)
            cd src/fort/
            echo '========================================  Compiling GWSWEX fortran module  ========================================' &> ../../build.log
            gfortran -Wall -Wno-tabs -pedantic -c -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -O3 -fPIC -march=native -mtune=znver2 preGWSWEX.f90 >> ../../build.log 2>&1
            if [ $? -eq 0 ]; then
                echo $'========================================  Successfully compiled GWSWEX fortran module  ======================================== \n\n\n' >> ../../build.log
                echo '========================================  Compiling GWSWEX python wrapper  ========================================' >> ../../build.log
                export LDFLAGS=-Wl,-rpath=../../libs/
                export NPY_DISTUTILS_APPEND_FLAGS=1
                (f2py --verbose -c -m gwswex_wrapper --no-lower --build-dir f2py_scratch --fcompiler=gnu95 --f90flags='-fopenmp -march=native -mtune=znver2' \
                    --opt='-O3' -I. -I/usr/local/include/fgsl/ -L/usr/local/lib/ -lfgsl -lgomp preGWSWEX.o GWSWEX.f90) >> ../../build.log 2>&1
                    if [ $? -eq 0 ]; then
                        echo $'========================================  Successfully compiled GWSWEX python wrapper  ======================================== \n\n\n' >> ../../build.log
                        (rm -f ../../libs/*.so 2>/dev/null && mv gwswex_wrapper*.so ../../libs/) \
                            || (mkdir -p ../../libs && mv gwswex_wrapper*.so ../../libs/ && echo 'library directory created' >> ../../build.log)
                        echo 'placed the program in libs/' >> ../../build.log
                        echo 'cleaning up' >> ../../build.log
                        rm -rf *.mod *.o *.pyf f2py_scratch 2>/dev/null
                        cd ../../
                        echo 'Build Successful'
                    else
                        echo $'========================================  Failed to compile GWSWEX python wrapper  ======================================== \n\n\n' >> build.log
                        echo 'Build Failed!'
                        exit 1
                    fi
            else
                echo $'========================================  Failed to compile GWSWEX fortran module  ======================================== \n\n\n' >> build.log
                echo 'Build Failed!'
                exit 1
            fi
            ;;
        n)
            echo 'compiling without OpenMP'
            cd src/fort/
            gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -o GWSWEX -O3 -xtarget=amd64 -xarch=amd64
            if [ $? -eq 0 ]; then
                echo 'successfully compiled GWSWEX'
                (rm -f ../../exe/fort/* 2>/dev/null && mv GWSWEX ../../exe/fort) || (mkdir -p ../../exe/fort && mv GWSWEX ../../exe/fort && echo 'dir created')
                echo 'placed the program in exe/fort'
                echo 'cleaning up'
                rm -f *.mod
                cd ../../
            else
                echo 'GWSWEX compilation failed'
            fi
            ;;
        h)
            help
            ;;
        \?)
            echo 'invalid option: -$OPTARG' >&2
            exit 1
            ;;
        *)
            echo 'error in command line parsing' >&2
            exit 1
            ;;
    esac
done

# Defaults to compiling with OpenMP
if (( $OPTIND == 1 )); then
	cd src/fort/
	gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX -O3 -xtarget=amd64 -xarch=amd64
	if [ $? -eq 0 ]; then
		echo 'successfully compiled GWSWEX (with OpenMP by default)'
		(rm -f ../../exe/fort/* 2>/dev/null && mv GWSWEX ../../exe/fort) || (mkdir -p ../../exe/fort && mv GWSWEX ../../exe/fort && echo 'dir created')
		echo 'placed the program in exe/fort'
        echo 'cleaning up'
        rm -f *.mod
		cd ../../
	else
		echo 'GWSWEX compilation failed'
	fi
fi