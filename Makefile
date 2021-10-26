###########################

VERSION=0.1.3

# for mkl
# make sure libiomp5 can be found
# MKLROOT       = /home/zilong/intel/oneapi/mkl/latest
# MKLROOT       = /opt/intel/oneapi/mkl/latest
MKLROOT       =

# install openblas lapack on mac with brew install openblas lapack
# OPENBLAS_ROOT = /usr/local/opt/openblas
# LAPACK_ROOT   = /usr/local/opt/lapack
OPENBLAS_ROOT =
LAPACK_ROOT   = /usr/local
STATIC       := 0

####### INC, LPATHS, LIBS, MYFLAGS
program       = PCAone
# for mac user, please change this to gnu gcc instead of the default clang version
# brew install gcc && ln -s $(which g++-11) /usr/local/bin/g++
CXX          ?= g++       # use default g++ only if not set in env
CXXFLAGS	  = -O3 -Wall -std=c++11 -mavx2 -mfma -ffast-math -m64
MYFLAGS       = -DVERSION=\"$(VERSION)\" -DNDEBUG
INC           = -I./external -I/usr/local/include
LPATHS        = -L/usr/local/lib
# CURRENT_DIR   = $(shell pwd)
SLIBS         = ./external/zstd/lib/libzstd.a  
# detect OS architecture and add flags
Platform     := $(shell uname -s)
ifeq ($(Platform),Linux)
###### for linux
	ifeq ($(strip $(STATIC)),1)
		CXXFLAGS += -static-libstdc++ 
		# CXXFLAGS += -static
        SLIBS    += /usr/lib/x86_64-linux-gnu/libz.a
	else
		CXXFLAGS += -march=native
        DLIBS    += -lz 
	endif

	MYFLAGS  += -fopenmp
	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		LPATHS  += -L${MKLROOT}/lib
		ifeq ($(strip $(STATIC)),1)
			SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
			# SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread
		else
			DLIBS += -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
		endif

	else ifneq ($(strip $(OPENBLAS_ROOT)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC     += -I${OPENBLAS_ROOT}/include -I${LAPACK_ROOT}/include
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib
		DLIBS   += -llapack -llapacke -lopenblas -lgfortran -lgomp -lpthread

	else
		DLIBS   += -lgomp -lpthread
	endif

else ifeq ($(Platform),Darwin)
###### for mac
	ifeq ($(strip $(STATIC)),1)
		SLIBS += /usr/local/opt/zlib/lib/libz.a /usr/local/lib/libomp.a  # clang needs libomp.a
		CXXFLAGS += -stdlib=libc++ -Xpreprocessor -fopenmp 
	else
        DLIBS += -lz
		CXXFLAGS += -march=native -fopenmp
	endif

	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		LPATHS  += -L${MKLROOT}/lib
		ifeq ($(strip $(STATIC)),1)
			SLIBS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread
			# SLIBS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread
		else
			DLIBS += -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
		endif

	else ifneq ($(strip $(OPENBLAS_ROOT)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC     += -I${OPENBLAS_ROOT}/include -I${LAPACK_ROOT}/include
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib
		DLIBS   += -llapack -llapacke -lopenblas -lgfortran -lgomp -lpthread

	else
		DLIBS   += -lgomp -lpthread
	endif

endif

OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))

LIBS += ${SLIBS} ${DLIBS} -lm -ldl

.PHONY: clean

all: ${program}

${program}: zstdlib bgenlib ${OBJ}
	$(CXX) $(CXXFLAGS) -o $(program) ${OBJ} ./external/bgen/bgenlib.a ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

zstdlib:
	(cd ./external/zstd/lib/; $(MAKE))

bgenlib:
	(cd ./external/bgen/; $(MAKE))

clean:
	(rm -f $(OBJ) $(program))
	(cd ./external/bgen/; $(MAKE) clean)
	(cd ./external/zstd/lib/; $(MAKE) clean)
