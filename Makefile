###########################

VERSION=0.1.0

BGEN_PATH     =
# for mkl on mac
# make sure libiomp5 can be found on your computer
MKLROOT       = /opt/intel/oneapi/mkl/latest
LIBIOMP5      = /usr/local/lib

OPENBLAS_ROOT =
LAPACK_ROOT   =
STATIC       := 1

#######
program       = PCAone
# for mac user, please change this to gnu gcc instead of the default clang version
CXX           = g++
CXXFLAGS	  = -O3 -Wall -std=c++11 -ftree-vectorize -ffast-math -fPIC -fopenmp
MYFLAGS       = -DVERSION=\"$(VERSION)\" -DNDEBUG
INC           = -I./external
#######
# INC, LPATHS, LIBS, MYFLAGS

# detect OS architecture and add flags
Platform      := $(shell uname -s)
ifeq ($(Platform),Linux)
###### for linux
	ifeq ($(strip $(STATIC)),1)
		CXXFLAGS += -Wl,-Bdynamic -static-libgcc -static-libstdc++
	else
		CXXFLAGS += -march=native
	endif

	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		CXXFLAGS += -m64
		ifeq ($(strip $(STATIC)),1)
			LIBS += -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_gnu_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lpthread -ldl -lm -lgomp
		else
			LIBS    += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread
			LPATHS  += -L${MKLROOT}/lib
		endif

	else ifneq ($(strip $(OPENBLAS_ROOT)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC      += -I${OPENBLAS_ROOT}/include/
		LIBS    += -Wl,-rpath=${OPENBLAS_ROOT}/lib/ -llapack -llapacke -lopenblas -lgfortran
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib

	endif

else ifeq ($(Platform),Darwin)
###### for mac
	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		LPATHS  += -L${MKLROOT}/lib -L${LIBIOMP5}
		ifeq ($(strip $(STATIC)),1)
			LIBS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm
			# LIBS += ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm
		else
			LIBS    +=  -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${LIBIOMP5} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
		endif
		CXXFLAGS += -m64

	else ifneq ($(strip $(OPENBLAS_ROOT)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		INC     += -I${OPENBLAS_ROOT}/include/
		LIBS    += -llapack -llapacke -lopenblas
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib

	endif
    CXXFLAGS    += -static-libstdc++

endif

# OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))
OBJ = \
	src/Main.o \
	src/Data.o \
	src/Utils.o \
	src/Arnoldi.o \
	src/Halko.o \
	src/FilePlink.o \
	src/FileBeagle.o \

ifneq ($(strip $(BGEN_PATH)), )
    MYFLAGS += -DWITH_BGEN
    INC     += -I${BGEN_PATH}/genfile/include -I${BGEN_PATH}/3rd_party/zstd-1.1.0 -I${BGEN_PATH}/3rd_party/zstd-1.1.0/lib
    LPATHS  += -L${BGEN_PATH}/build/ -L${BGEN_PATH}/build/3rd_party/zstd-1.1.0/
    LIBS    += -lbgen -lzstd
	OBJ     += src/FileBgen.o
endif

LIBS        += -lz

PGEN_PATH     = ./external/pgenlib/
PGEN_OBJECTS  = $(patsubst %.cc,%.o,$(wildcard ${PGEN_PATH}include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard ${PGEN_PATH}*.cpp))
OBJECTS       = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) ${PGEN_OBJECTS}

.PHONY: clean

all: ${program}

${program}: ${OBJ}
	$(CXX) $(CXXFLAGS) ${MYFLAGS} -o $(program) ${OBJ} ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f $(OBJ) $(program)