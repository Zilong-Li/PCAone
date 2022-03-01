######################### configure ################

# for mkl
# make sure libiomp5 can be found
# MKLROOT       = /home/zilong/intel/oneapi/mkl/latest
# MKLROOT       = /opt/intel/oneapi/mkl/latest
MKLROOT       =

# install openblas lapack on mac with brew install openblas lapack
# OPENBLAS_ROOT = /usr/local/opt/openblas
# LAPACK_ROOT   = /usr/local/opt/lapack
OPENBLAS_ROOT =
LAPACK_ROOT   =

# by default dynamical linking
STATIC       := 0
# only if static = 1, IOMP5 works
IOMP5        := 0

########################### end ###########################

VERSION=0.1.6
# detect OS architecture and add flags
Platform     := $(shell uname -s)

$(info "building PCAone on ${Platform} -- version ${VERSION}")


####### INC, LPATHS, LIBS, MYFLAGS
program       = PCAone
# for mac user, please change this to gnu gcc instead of the default clang version
# brew install gcc && ln -s $(which g++-11) /usr/local/bin/g++
# use default g++ only if not set in env
CXX           ?= g++
CXXFLAGS	  = -O3 -Wall -std=c++11 -mavx2 -mfma -ffast-math -m64 -fPIC -pipe
MYFLAGS       = -DVERSION=\"$(VERSION)\" -DNDEBUG
LINKFLAGS     = -s
CFLAGS        =
# CURRENT_DIR   = $(shell pwd)
INC           = -I./external -I./external/zstd/lib
LPATHS        = -L/usr/local/lib

ifeq ($(strip $(STATIC)),1)
	INC  += -I./external/zstd/lib
	SLIBS += ./external/zstd/lib/libzstd.a

	ifeq ($(Platform), Darwin)
		SLIBS += /usr/local/opt/zlib/lib/libz.a
		ifeq ($(strip $(IOMP5)), 1)
			SLIBS += /usr/local/opt/gcc/lib/gcc/11/libgomp.a  # gcc need libgomp.a
			CXXFLAGS += -static-libgcc
			MYFLAGS  += -fopenmp
		else
			SLIBS += /usr/local/lib/libomp.a  # clang needs libomp.a
			CXXFLAGS += -stdlib=libc++
			CFLAGS  += -Xpreprocessor -fopenmp
		endif
	else
		SLIBS    += /usr/lib/x86_64-linux-gnu/libz.a
		ifeq ($(strip $(IOMP5)), 1)
			CXXFLAGS += -static-libgcc -static-libstdc++
		else
			CXXFLAGS += -static
		endif
	endif
else
	CXXFLAGS += -march=native
	DLIBS    += -lz -lzstd
	INC      += -I/usr/local/include
endif


ifeq ($(Platform),Linux)
###### for linux
	MYFLAGS  += -fopenmp
	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		LPATHS  += -L${MKLROOT}/lib
		ifeq ($(strip $(STATIC)),1)
			ifeq ($(strip $(IOMP5)), 1)
				SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
			else
				SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread
			endif
		else
			DLIBS += -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
		endif

	else ifneq ($(strip $(OPENBLAS_ROOT)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC     += -I${OPENBLAS_ROOT}/include -I${LAPACK_ROOT}/include
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib
		DLIBS   += -llapack -llapacke -lopenblas -lgfortran -lgomp -lpthread

	else
		CXXFLAGS += -static-libgcc -static-libstdc++  # helpful to fix some glibcxx issue
		DLIBS   += -lgomp -lpthread
	endif

else ifeq ($(Platform),Darwin)
###### for mac
	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		LPATHS  += -L${MKLROOT}/lib
		ifeq ($(strip $(STATIC)),1)
			ifeq ($(strip $(IOMP5)), 1)
				SLIBS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread
			else
				SLIBS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread
			endif
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
		MYFLAGS  += -fopenmp
	endif

endif

OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))

SLIBS += ./external/bgen/bgenlib.a

LIBS += ${SLIBS} ${DLIBS} -lm -ldl

.PHONY: clean

all: ${program}

${program}: zstdlib bgenlib ${OBJ}
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LINKFLAGS) -o $(program) ${OBJ}  ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

zstdlib:

ifeq ($(STATIC),1)
	(cd ./external/zstd/lib/; $(MAKE))
else
	@echo "no building zstd manually"
endif

bgenlib:
	(cd ./external/bgen/; $(MAKE))

rm:
	(rm -f $(OBJ) $(program))
	(cd ./external/bgen/; $(MAKE) clean)

clean:
	(rm -f $(OBJ) $(program))
	(cd ./external/bgen/; $(MAKE) clean)
	(cd ./external/zstd/lib/; $(MAKE) clean)
