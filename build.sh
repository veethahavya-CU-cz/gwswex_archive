#!/bin/bash
while getopts ":cpn" opt; do
	case $opt in
		n)
			cd src/fort/
			gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -o GWSWEX
			if [ $? -eq 0 ]; then
				echo "successfully compiled GWSWEX"
				mv GWSWEX ../../exe/fort
				cd ../../
			else
				echo "GWSWEX compilation failed"
			fi
			;;
		p)
			cd src/fort/
			gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -lgomp -fopenmp -o GWSWEX
			if [ $? -eq 0 ]; then
				echo "successfully compiled GWSWEX"
				mv GWSWEX ../../exe/fort
				cd ../../
			else
				echo "GWSWEX compilation failed"
			fi
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
	gfortran GWSWEX.f90 -L/usr/local/lib/ -lfgsl -I/usr/local/include/fgsl/ -o GWSWEX
	if [ $? -eq 0 ]; then
		echo "successfully compiled GWSWEX"
		mv GWSWEX ../../exe/fort
		cd ~/GWSWEX
	else
		echo "GWSWEX compilation failed"
	fi
fi