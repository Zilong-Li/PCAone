#
# Makefile to compile bgen library
#
Platform    := $(shell uname -s)
CXX         ?= g++
CXXFLAGS     = -O2 -Wall -std=c++11
AVX = 1

ifeq ($(Platform), Darwin)
	CXXFLAGS += -march=native
else
	CXXFLAGS += -mavx -mavx2
endif

OBJECTS      = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
INC          = -I../zstd/lib

bgenlib.a: ${OBJECTS}
	ar -rcs $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
