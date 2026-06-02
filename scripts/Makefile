program = PCAone-merge-plink
CXX ?= g++
CXXSTD ?= c++17
CXXFLAGS ?= -O3 -Wall -std=$(CXXSTD)

.PHONY: all clean

all: $(program)

$(program): MergePlink.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(program)
