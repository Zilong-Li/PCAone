###########################

VERSION=0.1.3

MKLROOT       =

####### INC, LPATHS, LIBS, MYFLAGS
program       = PCAone
# for mac user, please change CXX to g++ gcc instead of the default clang version
# brew install gcc && ln -s $(which g++-11) /usr/local/bin/g++
# export CC="gcc-11"
# export CXX="g++-11"
# export LDFLAGS="-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/11/"
# conda install -c anaconda mkl mkl-include intel-openmp
# CXX           = g++
CXXFLAGS	  = -O3 -Wall -std=c++11 -march=native -ffast-math -m64
MYFLAGS       = -fopenmp -DVERSION=\"$(VERSION)\" -DNDEBUG -DWITH_BGEN -DWITH_MKL -DEIGEN_USE_MKL_ALL
INC           = -I./external -I/usr/local/include
LPATHS        = -L/usr/local/lib

INC          += -I${MKLROOT}/include/
LPATHS       += -L${MKLROOT}/lib
DLIBS        += -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lz -lzstd


#######

OBJ = $(patsubst %.cpp, %.o, $(wildcard ./src/*.cpp))

LIBS += ${SLIBS} ${DLIBS} -lm -ldl


# PGEN_PATH     = ./external/pgenlib/
# PGEN_OBJECTS  = $(patsubst %.cc,%.o,$(wildcard ${PGEN_PATH}include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard ${PGEN_PATH}*.cpp))
# OBJ       = $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) ${PGEN_OBJECTS}

.PHONY: clean

all: ${program}

${program}: bgenlib ${OBJ}
	$(CXX) $(CXXFLAGS) -o $(program) ${OBJ} ./external/bgen/bgenlib.a ${LPATHS} ${LIBS}

%.o: %.cpp
	${CXX} ${CXXFLAGS} ${MYFLAGS} -o $@ -c $< ${INC}

bgenlib:
	(cd ./external/bgen/; $(MAKE))

clean:
	(rm -f $(OBJ) $(program);cd ./external/bgen/; $(MAKE) clean)
