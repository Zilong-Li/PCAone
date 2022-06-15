######################### configure ################

# for mkl
# make sure libiomp5 can be found
# MKLROOT=/home/zilong/intel/oneapi/mkl/latest
# MKLROOT=/opt/intel/oneapi/mkl/latest
MKLROOT       =

# install openblas lapack on mac with brew install openblas lapack
# OPENBLAS_ROOT = /usr/local/opt/openblas
# LAPACK_ROOT   = /usr/local/opt/lapack
OPENBLAS_ROOT =
LAPACK_ROOT   =

# by default dynamical linking
STATIC        = 0
# only if static = 1, IOMP5 works
IOMP5         = 0

########################### end ###########################

VERSION=0.1.9
# detect OS architecture and add flags
Platform     := $(shell uname -s)

$(info "building PCAone on ${Platform} -- version ${VERSION}")


####### INC, LPATHS, LIBS, MYFLAGS
program       = PCAone
# use default g++ only if not set in env
CXX           ?= g++
CXXFLAGS	  += -O3 -Wall -std=c++11 -ffast-math -m64 -fPIC
MYFLAGS        = -DVERSION=\"$(VERSION)\" -DNDEBUG
LDFLAGS       += -s  # this is obsolete and igonored on mac
# CURRENT_DIR   = $(shell pwd)
INC           = -I./external -I./external/zstd/lib
# LPATHS        = -L/usr/local/lib
PCALIB = libpcaone.a
AVX = 1

ifeq ($(strip $(AVX)),1)
  $(info "use -mavx2 for PCAone")
  CXXFLAGS += -mavx2 -mfma
endif

ifeq ($(strip $(STATIC)),1)
	INC  += -I./external/zstd/lib
	ifeq ($(Platform), Darwin)
		SLIBS += /usr/local/opt/zlib/lib/libz.a
		ifeq ($(strip $(IOMP5)), 1)
			SLIBS += /usr/local/opt/gcc/lib/gcc/11/libgomp.a  # gcc need libgomp.a
			CXXFLAGS += -static-libgcc
			MYFLAGS  += -fopenmp
		else
			SLIBS += /usr/local/lib/libomp.a  # clang needs libomp.a
			CXXFLAGS += -stdlib=libc++
			# CFLAGS += -Xpreprocessor -fopenmp
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
	DLIBS    += -lz
endif


ifeq ($(Platform),Linux)
###### for linux
	MYFLAGS  += -fopenmp
	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		LPATHS  += -L${MKLROOT}/lib -L${MKLROOT}/lib/intel64
		ifeq ($(strip $(STATIC)),1)
			ifeq ($(strip $(IOMP5)), 1)
				SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
			else
				SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread
			endif
		else
			DLIBS += -Wl,--no-as-needed -Wl,-rpath,${MKLROOT}/lib,-rpath,${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
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
	MYFLAGS  += -Xpreprocessor -fopenmp
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
		DLIBS += -lomp -lpthread
	endif

endif

# OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))
# tar OBJ as libpcaone.a
OBJ = src/Arnoldi.o src/Halko.o src/Data.o src/Utils.o \
		src/FileBeagle.o src/FileCsv.o src/FileBgen.o src/FilePlink.o


SLIBS += ./external/zstd/lib/libzstd.a ./external/bgen/bgenlib.a

LIBS += ${SLIBS} ${DLIBS} -lm -ldl

.PHONY: all clean

all: ${program}

${program}: zstdlib bgenlib pcaonelib src/Main.o
	$(CXX) $(CXXFLAGS) -o $(program) src/Main.o ${PCALIB} ${LPATHS} ${LIBS} ${LDFLAGS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

zstdlib:
	(cd ./external/zstd/lib/; $(MAKE) ZSTD_LIB_COMPRESSION=0 ZSTD_LIB_DICTBUILDER=0)

bgenlib:
	(cd ./external/bgen/; $(MAKE))

pcaonelib:$(OBJ)
	ar -rcs $(PCALIB) $?

rm:
	(rm -f src/*.o $(program))
	(cd ./external/bgen/; $(MAKE) clean)

clean:
	(rm -f src/*.o $(program))
	(cd ./external/bgen/; $(MAKE) clean)
	(cd ./external/zstd/lib/; $(MAKE) clean)

