######################### configure ################

# for mkl
# MKLROOT=$HOME/intel/oneapi/mkl/latest # for linux
# MKLROOT=/opt/intel/oneapi/mkl/latest  # for mac
# make sure libiomp5 can be found on mac
# ln -s /opt/intel/oneapi/compiler/latest/mac/compiler/lib/libiomp5.dylib /opt/intel/oneapi/mkl/latest/lib/
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

VERSION=0.4.5
# detect OS architecture and add flags
Platform     := $(shell uname -s)

$(info "building PCAone on ${Platform} -- version ${VERSION}")


####### INC, LPATHS, LIBS, MYFLAGS
program       = PCAone
CXX           ?= g++    # use default g++ only if not set in env
CXXSTD         = c++11  # default c++11 if not set by the user
CXXFLAGS	  += -O3 -Wall -std=$(CXXSTD) -m64 -fPIC
MYFLAGS        = -DVERSION=\"$(VERSION)\"
LDFLAGS       += -s  # this is obsolete and will be igonored on mac
# CURRENT_DIR   = $(shell pwd)
INC           = -I./external -I./external/zstd/lib
PCALIB = libpcaone.a
AVX = 1
DEBUG = 0

ifeq ($(strip $(AVX)),1)
  $(info "use -mavx2 for PCAone")
  CXXFLAGS += -mavx2 -mfma
endif

ifeq ($(strip $(DEBUG)),1)
  $(info "use -mavx2 for PCAone")
  MYFLAGS += -DDEBUG
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
		SLIBS    += /usr/lib64/libz.a
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
	INC      += -I/usr/local/include
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
		DLIBS += -lomp -lpthread -L/usr/local/lib
	endif

endif

# OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))
# tar OBJ as libpcaone.a
OBJ = src/Arnoldi.o src/Halko.o src/Data.o src/Utils.o src/Cmd.o \
		src/FileBeagle.o src/FileCsv.o src/FileBgen.o src/FilePlink.o \
		src/FileBinary.o src/LD.o

SLIBS += ./external/bgen/bgenlib.a ./external/zstd/lib/libzstd.a

LIBS += ${SLIBS} ${DLIBS} -lm -ldl

.PHONY: all clean ld_matrix ld_prune ld_clump

all: ${program}

${program}: zstdlib bgenlib pcaonelib src/Main.o
	$(CXX) $(CXXFLAGS) -o $(program) src/Main.o ${PCALIB} ${LPATHS} ${LIBS} ${LDFLAGS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

zstdlib:
	(cd ./external/zstd/lib/; $(MAKE) lib-nomt)

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

data:
	wget http://popgen.dk/zilong/datahub/pca/example.tar.gz
	tar -xzf example.tar.gz && rm -f example.tar.gz

ld_matrix:
	./PCAone -b example/plink -k 3 --ld -o adj -m 1
	./PCAone -b example/plink -k 3 --ld -o pcaone
	diff adj.kept.bim pcaone.kept.bim
	cut -f1 adj.kept.bim | sort -cn  ## check if sorted
	awk '$$1==3' adj.kept.bim | cut -f4 | sort -cn

ld_prune:
	./PCAone -B adj.residuals --ld-bim adj.kept.bim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --ld-bim adj.kept.bim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 1
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null

ld_clump:
	./PCAone -B adj.residuals \
	--ld-bim adj.kept.bim \
	--clump example/plink.pheno0.assoc,example/plink.pheno1.assoc  \
	--clump-p1 0.01 \
	--clump-p2 0.05 \
	--clump-r2 0.1 \
	--clump-bp 10000000 \
	-o adj_clump_m0

	./PCAone -B adj.residuals \
	--ld-bim adj.kept.bim \
	--clump example/plink.pheno0.assoc,example/plink.pheno1.assoc  \
	--clump-p1 0.01 \
	--clump-p2 0.05 \
	--clump-r2 0.1 \
	--clump-bp 10000000 \
	-o adj_clump_m1 -m 1
