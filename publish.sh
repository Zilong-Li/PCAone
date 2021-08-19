#!/bin/bash

platform=$(uname -s)

if [ $platform == "Darwin" ];then
    echo "Publishing releases on MacOS";
    awk 'NR==9 || NR==75 {$0="#"$0}; NR==8 || NR==76 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >mac.makefile
    make -f mac.makefile clean && make -f mac.makefile 2>/dev/null && zip PCAone_MacOS.zip PCAone
    echo "Publishing PCAone_MacOS.zip done";
    awk 'NR==9 {$0="#"$0}; NR==8 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >mac.makefile
    make -f mac.makefile clean && make -f mac.makefile 2>/dev/null && zip PCAone_MacOS_MKL.zip PCAone
    echo "Publishing PCAone_MacOS_MKL.zip done";

elif [ $platform == "Linux" ];then
    echo "Publishing releases on Linux";
    awk 'NR==9 || NR==46 {$0="#"$0}; NR==7 || NR==47 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >linux.makefile
    make -f linux.makefile clean && make -f linux.makefile 2>/dev/null && zip PCAone_Linux.zip PCAone
    echo "Publishing PCAone_Linux.zip done";
    awk 'NR==9 {$0="#"$0}; NR==7 {sub("# ", "")}; NR==16 {sub("0", "1")}; 1' Makefile >linux.makefile
    make -f linux.makefile clean && make -f linux.makefile 2>/dev/null && zip PCAone_Linux_MKL.zip PCAone
    echo "Publishing PCAone_Linux_MKL.zip done";

fi

