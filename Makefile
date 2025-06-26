######################### configure ################

VERSION=0.5.2
# Detect OS and CPU architecture to add flags conditionally
Platform        := $(shell uname -s)
ARCH            := $(shell uname -m)

$(info "building PCAone on ${Platform} (${ARCH}) -- version ${VERSION}")

# MKLROOT=/opt/intel/oneapi/mkl/latest
MKLROOT            =
# /opt/intel/oneapi/compiler/latest/lib/
ONEAPI_COMPILER    := ${MKLROOT}
ONEAPI_OMP5        := $(ONEAPI_COMPILER)/lib/libiomp5.a
# make sure libiomp5 can be found on mac
# ln -s /opt/intel/oneapi/compiler/latest/mac/compiler/lib/libiomp5.dylib /opt/intel/oneapi/mkl/latest/lib/

# for MacOS use accelerate framework in default
# install openblas lapack on mac with brew install openblas lapack
OPENBLAS    =


STATIC          = 0  # by default dynamical linking
IOMP5           = 1  # use libiomp5 for mkl
DEBUG           = 0  # debug
AVX             =    # 1: enable avx2 fma, 0: do nothing

########################### end ###########################


####### INC, LPATHS, LIBS, MYFLAGS
program         = PCAone
CXX             ?= g++    # use default g++ only if not set in env
CXXSTD          = c++17
# Base compiler flags, applied to all architectures
CXXFLAGS        += -O3 -Wall -std=$(CXXSTD)

# Apply the -m64 flag only on the x86_64 architecture.
# This flag is not valid on aarch64.
ifeq ($(ARCH), x86_64)
	CXXFLAGS += -m64
endif

MYFLAGS         = -DVERSION=\"$(VERSION)\" -DNDEBUG
INC             = -I./external -I./external/zstd/lib
PCALIB = libpcaone.a

ifeq ($(Platform), Darwin)
	CXXFLAGS += -march=native
else
	ifeq ($(strip $(AVX)),)
		# Only apply -march=native on x86_64
		ifeq ($(ARCH), x86_64)
			CXXFLAGS += -march=native
		endif
	else ifeq ($(strip $(AVX)),1)
		# Only add AVX flags on x86_64, as they are not valid on aarch64
		ifeq ($(ARCH), x86_64)
			CXXFLAGS += -mavx2 -mfma
		endif
	endif
endif

ifeq ($(strip $(DEBUG)),1)
	MYFLAGS += -DDEBUG
endif


ifeq ($(Platform),Linux)
###### for linux
	MYFLAGS  += -fopenmp
	ifneq ($(strip $(MKLROOT)),)
		MYFLAGS += -DWITH_MKL -DEIGEN_USE_MKL_ALL
		INC     += -I${MKLROOT}/include/
		ifeq ($(strip $(STATIC)),1)
			ifeq ($(strip $(IOMP5)), 1)
				SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${ONEAPI_OMP5} -Wl,--end-group -L${ONEAPI_COMPILER}/lib
			else
				SLIBS += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp
			endif
		else
			DLIBS += -Wl,--no-as-needed -Wl,-rpath,${ONEAPI_COMPILER}/lib,-rpath,${MKLROOT}/lib/intel64 -L${ONEAPI_COMPILER}/lib -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
		endif

	else ifneq ($(strip $(OPENBLAS)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		DLIBS   += -llapack -llapacke -lopenblas -lgfortran -lgomp

	else
		DLIBS   += -lgomp
	endif

else ifeq ($(Platform),Darwin)
###### for mac
	MYFLAGS  += -Xpreprocessor -fopenmp

	ifneq ($(strip $(OPENBLAS)),)
		MYFLAGS += -DWITH_OPENBLAS -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
		ifeq ($(strip $(STATIC)),1)
			# SLIBS += /usr/local/opt/gcc/lib/gcc/11/libgomp.a  # gcc need libgomp.a
			SLIBS += /opt/homebrew/opt/openblas/lib/libopenblas.a /opt/homebrew/opt/libomp/lib/libomp.a -L/opt/homebrew/opt/gcc/lib/gcc/14/ -lgfortran
		else
			DLIBS   += -lopenblas -lomp
		endif

	else
		MYFLAGS += -DEIGEN_USE_BLAS
		ifeq ($(strip $(STATIC)),1)
			SLIBS += /opt/homebrew/opt/libomp/lib/libomp.a -framework Accelerate
		else
			DLIBS   += -lomp -framework Accelerate
		endif

	endif

endif

## libz
ifeq ($(strip $(STATIC)),1)
	ifeq ($(Platform), Darwin)
		SLIBS += /opt/homebrew/opt/zlib/lib/libz.a
		# CXXFLAGS += -stdlib=libc++
	else
		# Previous path: /usr/lib/x86_64-linux-gnu/libz.a
		# Default dynamic link below (-lz)
		SLIBS    += -lz
		CXXFLAGS += -static
		# CXXFLAGS += -static-libgcc -static-libstdc++
	endif
else
	DLIBS    += -lz
endif

OBJ = src/Arnoldi.o src/Halko.o src/Data.o src/Utils.o src/Cmd.o \
		src/FileBeagle.o src/FileCsv.o src/FileBgen.o src/FilePlink.o \
		src/FileBinary.o src/FileUSV.o src/LD.o src/Projection.o src/Inbreeding.o \
		src/kfunc.o

SLIBS += ./external/bgen/bgenlib.a ./external/zstd/lib/libzstd.a

LIBS += $(SLIBS) $(DLIBS) -lpthread -ldl -lm

.PHONY: all clean ld_matrix ld_r2 ld_prune ld_clump ld_tests test_full test_aarch64

all: ${program}

${program}: zstdlib bgenlib pcaonelib src/Main.o
	$(CXX) $(CXXFLAGS) -o $(program) src/Main.o $(PCALIB) $(LPATHS) $(LIBS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(MYFLAGS) -o $@ -c $< $(INC) $(CPPFLAGS)

zstdlib:
	$(MAKE) -C external/zstd/lib lib-nomt

bgenlib:
	$(MAKE) -C external/bgen CFLAGS='$(CFLAGS)' CXXFLAGS='$(CXXFLAGS)'

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

example_tests:
	./PCAone -g example/test.bgen -n4 -o m0
	./PCAone -g example/test.bgen -n4 -m0.1 -o m1
	diff m0.eigvals m1.eigvals
	./PCAone --csv example/BrainSpinalCord.csv.zst -k 10 -m 2 --scale 2 -S

hwe:
	./PCAone -b example/plink -k 3 -V
	./PCAone -b example/plink --USV pcaone --inbreed 1 -o inbreed_m0
	./PCAone -b example/plink --USV pcaone --inbreed 1 -o inbreed_m1 -m 1
	diff inbreed_m0.hwe inbreed_m1.hwe
	rm -f inbreed* pcaone.*

ld_matrix:
	./PCAone -b example/plink -k 3 --ld -o adj -d 2
	./PCAone -b example/plink -k 3 --ld -o pcaone -d 2 -m 2
	diff adj.mbim pcaone.mbim
	cut -f1 adj.mbim | sort -cn  ## check if sorted
	awk '$$1==3' adj.mbim | cut -f4 | sort -cn
	rm -f pcaone.*

ld_r2:
	./PCAone -B adj.residuals --match-bim adj.mbim --ld-bp 1000 --print-r2 -o adj_r2


ld_prune:
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 1
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null

ld_clump:
	./PCAone -B adj.residuals --match-bim adj.mbim --clump example/plink.pheno0.assoc --clump-p1 0.01 --clump-p2 0.05 --clump-r2 0.1 --clump-bp 10000000 -m 0 -o adj_clump_m0
	./PCAone -B adj.residuals --match-bim adj.mbim --clump example/plink.pheno0.assoc --clump-p1 0.01 --clump-p2 0.05 --clump-r2 0.1 --clump-bp 10000000 -m 1 -o adj_clump_m1
	diff adj_clump_m0.p0.clump adj_clump_m1.p0.clump > /dev/null

ld_tests:
	./PCAone -b example/plink -k 3 --ld -o adj -d 0 --maf 0.1
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 2
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null
	./PCAone -b example/plink -k 3 --ld -o adj -d 0 -m 4
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 2
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null
	./PCAone -b example/plink -k 3 --ld -o adj -d 1 --maf 0.1
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 2
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null
	./PCAone -b example/plink -k 3 --ld -o adj -d 1 -m 4
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 2
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null
	./PCAone -b example/plink -k 3 --ld -o adj -d 2 --maf 0.1
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 2
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null
	./PCAone -b example/plink -k 3 --ld -o adj -d 2 -m 4
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m0 -m 0
	./PCAone -B adj.residuals --match-bim adj.mbim  --ld-r2 0.8  --ld-bp 1000000 -o adj_prune_m1 -m 2
	diff adj_prune_m0.ld.prune.out adj_prune_m1.ld.prune.out > /dev/null

#################################################################
# CI/CD Test Suites
# Use 'make test_full' for complete testing (e.g., on x86_64)
# Use 'make test_aarch64' for a faster suite to avoid timeouts
#################################################################

# Define a faster version of example_tests that excludes the long-running CSV test
example_tests_fast:
	./PCAone -g example/test.bgen -n4 -o m0
	./PCAone -g example/test.bgen -n4 -m0.1 -o m1
	diff m0.eigvals m1.eigvals

# The fast test suite for aarch64
test_aarch64: data example_tests_fast hwe ld_matrix ld_r2 ld_clump ld_tests
	@echo "SUCCESS: aarch64 fast test suite completed."

# The complete test suite for other architectures like x86_64
test_full: data example_tests hwe ld_matrix ld_r2 ld_clump ld_tests
	@echo "SUCCESS: Full test suite completed."
