#
# Makefile to compile PGEN library
#
AR         ?= ar
CXX        ?= g++
CXXFLAGS    = -O3 -Wall -std=c++11
CFLAGS      = 

UNAME_S := $(shell uname -s)
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_S),Darwin)
	ifeq ($(UNAME_P),x86_64)
		CXXFLAGS += -arch x86_64
	endif
endif

OBJECTS       = $(patsubst %.cc,%.o,$(wildcard ./include/*.cc)) $(patsubst %.cpp,%.o,$(wildcard *.cpp))  $(patsubst %.cc,%.o,$(wildcard *.cc))
INC           =  -I./simde/ -I./include/

pgenlib.a: ${OBJECTS}
	${AR} rcs $@ $^ 

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a ./include/*.o
