#
# Makefile to compile bgen library
#
CXX          = g++-11
CXXFLAGS     = -O3 -Wall -std=c++11 -mavx -mavx2

OBJECTS      = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
INC          = -I/usr/include -I/usr/local/include

bgenlib.a: ${OBJECTS}
	ar -rcs $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} -o $@ -c $< ${INC}

clean:
	rm -f *.o *.a
