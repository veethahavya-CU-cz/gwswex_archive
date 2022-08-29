#!/bin/bash

function help() {
    echo "usage: ./build.sh [-p] [-n]"
    echo " -d: delete all output files/directories"
    echo " -r: rename output files"
    echo " -h: print help"
    exit 0
}

while getopts ":drh" opt; do
    case $opt in
        d)
            echo "removing old figs and I/O files"
            rm -rf output/fort/figs/*.*
            rm -rf exe/fort/input
            rm -rf exe/fort/output
            rm output/*.npz
            ;;
        r)
            cd output/figs/
            echo "renaming and saving outputs"
            echo "enter the name of the simulation: "
            read simname
            for f in *;do
                mv -v "$f" "${f%.*}_$simname.${f##*.}"
            done
            cd ..
            for f in *.npz;do
                mv -v "$f" "${f%.*}_$simname.${f##*.}"
            done
            ;;
        g)
            git fetch --prune
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

if (( $OPTIND == 1 )); then
    echo "choose a valid flag to perform option"
    help
fi