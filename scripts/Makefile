programs = PCAone-merge-plink PCAone-merge-pgen summarise_ld_r2bin
CXX ?= g++
CXXSTD ?= c++17
CXXFLAGS ?= -O3 -Wall -std=$(CXXSTD)

.PHONY: all clean

all: $(programs)

PCAone-merge-plink: MergePlink.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

PCAone-merge-pgen: MergePgen.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

summarise_ld_r2bin: summarise_ld_r2bin.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lz

clean:
	rm -f $(programs)
