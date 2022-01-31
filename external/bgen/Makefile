#
# Makefile to compile bgen library
#
CXX         ?= g++
CXXFLAGS     = -O3 -Wall -std=c++11 -mavx -mavx2

OBJECTS      = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
INC          = -I../zstd/lib

bgenlib.a: ${OBJECTS}
	ar -rcs $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
