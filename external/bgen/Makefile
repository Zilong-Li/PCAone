#
# Makefile to compile bgen library
#
# CXX         ?= g++

# Platform    := $(shell uname -s)
# AVX = 1
# ifeq ($(Platform), Darwin)
# 	CXXFLAGS += -march=native
# else
# 	CXXFLAGS += -mavx -mavx2
# endif

OBJECTS      = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
INC          = -I../zstd/lib

bgenlib.a: ${OBJECTS}
	ar -rcs $@ $^

%.o: %.cpp
	${CXX} ${CFLAGS} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
