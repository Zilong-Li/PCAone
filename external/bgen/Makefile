#
# Makefile to compile bgen library
#
CXX         ?= g++
CXXFLAGS    += -O3 -Wall -std=c++11
AVX = 1

ifeq ($(strip $(AVX)), 1)
  $(info "use -mavx2 for bgen")
  CXXFLAGS += -mavx -mavx2 -DAVX2
endif

OBJECTS      = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
INC          = -I../zstd/lib

bgenlib.a: ${OBJECTS}
	ar -rcs $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
