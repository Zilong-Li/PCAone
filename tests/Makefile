# detect OS architecture and add flags

VERSION=0.3.2
Platform     := $(shell uname -s)

CXX      = g++
CXXFLAGS = -std=c++11 -Wall -O3 -mavx2 -mfma
MYFLAGS  = -DVERSION=\"$(VERSION)\"
INC      = -I. -I../src -I../external
DLIBS   += -lz

ifeq ($(Platform),Linux)
	MYFLAGS += -fopenmp
	DLIBS   += -lgomp -lpthread
else ifeq ($(Platform),Darwin)
	MYFLAGS += -Xpreprocessor -fopenmp
	INC     += -I/usr/local/include
	DLIBS   += -lomp -lpthread -L/usr/local/lib
endif

SLIBS += ../external/zstd/lib/libzstd.a ../external/bgen/bgenlib.a ../src/libpcaone.a
LIBS  += ${SLIBS} ${DLIBS} -lm -ldl

OBJS     = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
# OBJS = test-main.o test-io.o test-phaseless.o test-fastphase.o test-forward-backward.o

BINS = test-shuffle-csvzstd.bin

.PHONY: all test clean

all: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) ${MYFLAGS} -o $@ -c $< ${INC}

%.bin: %.o zstdlib bgenlib test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

zstdlib:
	(cd ../external/zstd/lib/; $(MAKE) lib-nomt)

bgenlib:
	(cd ../external/bgen/; $(MAKE))

test: $(BINS)
	./test-shuffle-csvzstd.bin --success

clean:
	rm -f *.o *.bin
