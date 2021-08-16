###########################

VERSION=0.1.0

# for mkl on mac
# make sure libiomp5 can be found on your computer
# MKLROOT       = /opt/intel/oneapi/mkl/latest
# LIBIOMP5      = /opt/intel/oneapi/compiler/latest/mac/compiler/lib
MKLROOT       =
LIBIOMP5      =

# for openblas on mac
# brew install openblas lapack
# OPENBLAS_ROOT = /usr/local/opt/openblas
# LAPACK_ROOT   = /usr/local/opt/lapack
OPENBLAS_ROOT =
LAPACK_ROOT   =
STATIC       := 0

####### INC, LPATHS, LIBS, MYFLAGS
program       = PCAone
# for mac user, please change this to gnu gcc instead of the default clang version
CXX           = g++-11
CXXFLAGS	  = -O3 -Wall -std=c++11 -mavx -mavx2 -ffast-math -fopenmp
MYFLAGS       = -DVERSION=\"$(VERSION)\" -DNDEBUG -DWITH_BGEN
INC           = -I./external -I/usr/include -I/usr/local/include
LPATHS        = -L/usr/lib -L/usr/local/lib
#######

# detect OS architecture and add flags
Platform      := $(shell uname -s)
ifeq ($(Platform),Linux)
###### for linux
	ifeq ($(strip $(STATIC)),1)
		CXXFLAGS += -Wl,-Bdynamic -static-libgcc -static-libstdc++
        LIBS     += /usr/lib/x86_64-linux-gnu/libz.a /usr/lib/x86_64-linux-gnu/libzstd.a
	else
		CXXFLAGS += -march=native
        LIBS     += -lz -lzstd
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
		INC     += -I${OPENBLAS_ROOT}/include -I${LAPACK_ROOT}/include
		LIBS    += -Wl,-rpath=${OPENBLAS_ROOT}/lib/ -llapack -llapacke -lopenblas -lgfortran
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib

	endif

else ifeq ($(Platform),Darwin)
###### for mac
	ifeq ($(strip $(STATIC)),1)
		LIBS += /usr/local/opt/zlib/lib/libz.a /usr/local/lib/libzstd.a  # path to static lib
	else
        LIBS += -lz -lzstd
	endif

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
		INC     += -I${OPENBLAS_ROOT}/include -I${LAPACK_ROOT}/include
		LIBS    += -llapack -llapacke -lopenblas
		LPATHS  += -L${OPENBLAS_ROOT}/lib -L${LAPACK_ROOT}/lib

	endif
    CXXFLAGS    += -static-libstdc++

endif

OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))


# PGEN_PATH     = ./external/pgenlib/
# PGEN_OBJECTS  = $(patsubst %.cc,%.o,$(wildcard ${PGEN_PATH}include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard ${PGEN_PATH}*.cpp))
# OBJ       = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) ${PGEN_OBJECTS}

.PHONY: clean

all: ${program}

${program}: bgenlib ${OBJ}
	$(CXX) $(CXXFLAGS) ${MYFLAGS} -o $(program) ${OBJ} ./external/bgen/bgenlib.a ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

bgenlib:
	(cd ./external/bgen/; $(MAKE) clean && $(MAKE))

clean:
	rm -f $(OBJ) $(program)