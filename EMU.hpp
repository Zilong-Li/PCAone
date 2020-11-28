#ifndef __EMU__
#define __EMU__

#include <chrono>
#include <iostream>

using namespace std;

typedef unsigned int uint;
typedef unsigned long long uint64;

class MeasureTime {

  public:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    time_t start_time_info, end_time_info;

    void init() {
      auto start = std::chrono::system_clock::now(); // wall clock
      start_time_info = std::chrono::system_clock::to_time_t( start );
      begin = std::chrono::steady_clock::now(); // to measure elapsed time
    }

    void stop(){
      auto endtime = std::chrono::system_clock::now();
      end_time_info = std::chrono::system_clock::to_time_t( endtime );
      end = std::chrono::steady_clock::now();
    }

    MeasureTime(void);
    ~MeasureTime(void);
};

struct Param {
    string intype = ""; // bfile, pfile or bgen
    string bed_prefix;
    string pgen_prefix;
    string bgen;
    uint k;
    bool batch = false; // if load all genotypes into RAM.
    uint threads = 1;
    uint iter = 100;
    // for arnoldi iteration;
    uint maxiter = 500;
    double tol = 1e-6;
};

void parse_params(int argc, char* argv[], struct Param* params);

#endif
