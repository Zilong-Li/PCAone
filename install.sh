#!/usr/bin/env bash

##########
# first, this script will install mkl dependency via conda
# second, it will try to add or link libiomp5 into runtime path

set -e

version=0.1.6

abort() {
  printf "%s\n" "$@"
  exit 1
}

# get conda prefix
if [ ${CONDA_PREFIX} ];then
    PREFIX=${CONDA_PREFIX}
else
    PREFIX=$(conda info --base)
fi
echo "try to install dependency in the environment $PREFIX "
# install mkl via conda
#
conda install -c conda-forge -y mkl

# check OS.
OS="$(uname)"

if [[ "${OS}" == "Linux" ]]
then
    system="Linux"
    echo "download PCAone for ${system}"
    url="https://github.com/Zilong-Li/PCAone/releases/download/v${version}/PCAone-v${version}_x64-${system}-iomp5.zip"
    curl -OL $url && unzip "PCAone-v${version}_x64-${system}-iomp5.zip"
    # download patchelf
    conda install -c conda-forge -y patchelf
    # add rpath
    patchelf --set-rpath ${PREFIX}/lib PCAone
    # try to link
    [[ $? != 0 ]] && link -sf ${PREFIX}/lib/libiomp5.so /usr/local/lib

elif [[ "${OS}" == "Darwin" ]]
then
    system="Mac"
    echo "download PCAone for ${system}"
    url="https://github.com/Zilong-Li/PCAone/releases/download/v${version}/PCAone-v${version}_x64-${system}-iomp5.zip"
    curl -OL $url && unzip "PCAone-v${version}_x64-${system}-iomp5.zip"
    # add rpath
    install_name_tool -add_rpath ${PREFIX}/lib PCAone
    # try to link
    [[ $? != 0 ]] && link -sf ${PREFIX}/lib/libiomp5.dylib /usr/local/lib
fi
