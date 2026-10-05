# Installation

There are 3 ways to install PCAone.

## Download compiled binary

There are compiled binaries provided for both Linux and Mac platform. Check
[the releases page](https://github.com/Zilong-Li/PCAone/releases) to download one.

```shell
pkg=https://github.com/Zilong-Li/PCAone/releases/latest/download/PCAone-Linux.zip
wget $pkg || curl -LO $pkg
unzip -o PCAone-Linux.zip
```

## Via Conda

PCAone is also available from [bioconda](https://anaconda.org/bioconda/pcaone).

```shell
conda config --add channels bioconda
conda install pcaone
PCAone --help
```

## Build from source

PCAone has been tested on both `Linux` and `MacOS` system. To build PCAone from the source code, the following dependencies are required:

- GCC/Clang compiler with C++17 support
- GNU make
- zlib

On Linux, we **recommend** building the software from source with MKL as backend to maximize the performance.

### With MKL or OpenBLAS as backend

Build PCAone dynamically with MKL can maximize the performance for large
dataset particularly, because the faster threading layer `libiomp5` will be
linked at runtime. There are two options to obtain MKL library:

- download `MKL` from [the website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)

After having `MKL` installed, find the `MKL` root path and replace the path below with your own.

```shell
make -j4 MKLROOT=/opt/intel/oneapi/mkl/latest \
         ONEAPI_COMPILER=/opt/intel/oneapi/compiler/latest
```

Alternatively, for advanced user, modify variables directly in `Makefile` and run `make` to use MKL or OpenBlas as backend.

- install `MKL` by conda

```shell
conda install -c conda-forge -c anaconda -y mkl mkl-include intel-openmp
git clone https://github.com/Zilong-Li/PCAone.git
cd PCAone
# if mkl is installed by conda then use ${CONDA_PREFIX} as mklroot
make -j4 MKLROOT=${CONDA_PREFIX}
./PCAone -h
```

### Without MKL or OpenBLAS dependency

If you don't want any optimized math library as backend, just run:

```shell
git clone https://github.com/Zilong-Li/PCAone.git
cd PCAone
make -j4
./PCAone -h
```

### macOS

For macOS users, install OpenMP first and build as above. The
[mac workflow](https://github.com/Zilong-Li/PCAone/blob/main/.github/workflows/mac.yml) shows a full build.

```shell
brew install libomp
make -j4
```

Commands in this guide using `./PCAone` assume a local binary; use `PCAone` if it is on your PATH.
