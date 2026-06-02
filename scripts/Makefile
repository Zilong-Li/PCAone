programs = PCAone-merge-plink PCAone-merge-pgen
CXX ?= g++
CXXSTD ?= c++17
CXXFLAGS ?= -O3 -Wall -std=$(CXXSTD)

.PHONY: all clean

all: $(programs)

PCAone-merge-plink: MergePlink.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

PCAone-merge-pgen: MergePgen.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(programs)
