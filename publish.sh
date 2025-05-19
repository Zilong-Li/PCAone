#!/bin/bash

# exit immediately if errors
set -e

version=v`grep "^VERSION" Makefile|cut -c 9-`

dir=releases/$version
mkdir -p $dir

platform=$(uname -s)

if [ $platform == "Darwin" ];then
    
    echo "Publishing releases on MacOS Silicon";
    export LDFLAGS="-L"$(brew --prefix libomp)/lib
    export CPPFLAGS="-I"$(brew --prefix libomp)/include
    make clean && make -j6 STATIC=1 && zip PCAone-Mac.zip PCAone PCAone.pdf && mv PCAone-Mac.zip $dir && echo "Publishing PCAone-Mac.zip done";

elif [ $platform == "Linux" ];then

    echo "Publishing releases on Linux ";
    make clean && make -j6 MKLROOT=/opt/intel/oneapi/mkl/latest ONEAPI_COMPILER=/opt/intel/oneapi/compiler/latest STATIC=1 AVX=1 && zip PCAone-Linux.zip PCAone PCAone.pdf && mv PCAone-Linux.zip $dir && echo "Publishing PCAone-Linux.zip done";

fi

