#!/bin/bash

set -e

version=v`grep "^VERSION" Makefile|cut -c 9-`

dir=releases/$version
mkdir -p $dir

platform=$(uname -s)

if [ $platform == "Darwin" ];then
    echo "Publishing releases on MacOS";

    echo "Publishing PCAone without libiomp5 threading";
    export CC="gcc"
    export CXX="g++"
    awk 'NR==7 {$0="#"$0}; NR==6 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >mac.makefile
    make -f mac.makefile clean && make -f mac.makefile AVX=1 2>/dev/null && zip PCAone-avx2-Mac.zip PCAone && mv PCAone-avx2-Mac.zip $dir && echo "Publishing PCAone-x64-Mac.zip done";

    echo "Publishing PCAone with libiomp5 threading";
    export CC=$(find $(brew --prefix)/bin/ -name "gcc-[0-9]*" | tail -1)
    export CXX=$(find $(brew --prefix)/bin/ -name "g++-[0-9]*" | tail -1)

    awk ' NR==7 {$0="#"$0}; NR==6 {sub("# ", "")}; NR==16 || NR==18 {sub("0", "1")}; 1' Makefile >mac2.makefile
    make -f mac2.makefile clean && make -f mac2.makefile AVX=1 2>/dev/null && zip PCAone-avx2-Mac-iomp5.zip PCAone && mv PCAone-avx2-Mac-iomp5.zip $dir && echo "Publishing PCAone-x64-Mac-iomp5.zip done";

elif [ $platform == "Linux" ];then
    echo "Publishing releases on Linux";

    echo "Publishing PCAone without libiomp5 threading";
    awk ' NR==7 {$0="#"$0}; NR==5 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >linux.makefile
    make -f linux.makefile clean && make -f linux.makefile 2>/dev/null && zip PCAone-x64-Linux.zip PCAone && mv PCAone-x64-Linux.zip $dir && echo "Publishing PCAone_Linux.zip PCAone-x64-Linux.zip done";
    make -f linux.makefile clean && make -f linux.makefile AVX=1 2>/dev/null && zip PCAone-avx2-Linux.zip PCAone && mv PCAone-avx2-Linux.zip $dir && echo "Publishing PCAone_Linux.zip PCAone-avx2-Linux.zip done";

    echo "Publishing PCAone with libiomp5 threading";
    awk 'NR==7 {$0="#"$0}; NR==5 {sub("# ", "")}; NR==16 || NR==18 {sub("0", "1")}; 1' Makefile >linux2.makefile
    make -f linux2.makefile clean && make -f linux2.makefile 2>/dev/null && zip PCAone-x64-Linux-iomp5.zip PCAone && mv PCAone-x64-Linux-iomp5.zip $dir && echo "Publishing PCAone-x64-Linux-iomp5.zip done";
    make -f linux2.makefile clean && make -f linux2.makefile AVX=1 2>/dev/null && zip PCAone-avx2-Linux-iomp5.zip PCAone && mv PCAone-avx2-Linux-iomp5.zip $dir && echo "Publishing PCAone-avx2-Linux-iomp5.zip done";

fi

