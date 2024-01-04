#!/bin/bash

# exit immediately if errors
set -e

version=v`grep "^VERSION" Makefile|cut -c 9-`

dir=releases/$version
mkdir -p $dir

platform=$(uname -s)

if [ $platform == "Darwin" ];then
    echo "Publishing releases on MacOS avx2";
    echo "Publishing PCAone without libiomp5 threading";
    export CC="gcc"
    export CXX="g++"
    make clean && make -j6 MKLROOT=/opt/intel/oneapi/mkl/latest STATIC=1 && zip PCAone-avx2-Mac.zip PCAone manual.pdf && mv PCAone-avx2-Mac.zip $dir && echo "Publishing PCAone-avx2-Mac.zip done";

    # echo "Publishing PCAone with libiomp5 threading";
    # export CC=$(find $(brew --prefix)/bin/ -name "gcc-[0-9]*" | tail -1)
    # export CXX=$(find $(brew --prefix)/bin/ -name "g++-[0-9]*" | tail -1)
    # export LDFLAGS="-L/usr/local/lib/"

    # awk ' NR==7 {$0="#"$0}; NR==6 {sub("# ", "")}; NR==16 || NR==18 {sub("0", "1")}; 1' Makefile >mac2.makefile

elif [ $platform == "Linux" ];then
    echo "Publishing releases on Linux x86_x64";
    make clean && make -j6 MKLROOT=$HOME/intel/oneapi/mkl/latest STATIC=1 AVX=0 && zip PCAone-x64-Linux.zip PCAone manual.pdf && mv PCAone-x64-Linux.zip $dir && echo "Publishing PCAone_Linux.zip PCAone-x64-Linux.zip done";
    make clean && make -j6 MKLROOT=$HOME/intel/oneapi/mkl/latest STATIC=1 AVX=1  && zip PCAone-avx2-Linux.zip PCAone manual.pdf && mv PCAone-avx2-Linux.zip $dir && echo "Publishing PCAone_Linux.zip PCAone-avx2-Linux.zip done";
    # awk ' NR==7 {$0="#"$0}; NR==5 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >linux.makefile
    # awk 'NR==7 {$0="#"$0}; NR==5 {sub("# ", "")}; NR==16 || NR==18 {sub("0", "1")}; 1' Makefile >linux2.makefile

fi

