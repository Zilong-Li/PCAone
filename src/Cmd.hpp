#ifndef CMD_H_
#define CMD_H_

#include <iostream>
#include <iterator>
#include "popl/popl.hpp"

using namespace popl;

typedef unsigned int uint;
typedef unsigned long long uint64;

enum class FileType
{
    None,
    PLINK,
    CSV,
    BEAGLE,
    BGEN
};


class Param
{
public:
    Param(int argc, char** argv);
    ~Param();

    FileType intype = FileType::None; // PLINK, CSV, BEAGLE, BGEN
    std::string bed_prefix = "";
    std::string pgen_prefix = "";
    std::string bgen = "";
    std::string beagle = "";
    std::string csvfile = "";
    std::string outfile = "pcaone";
    std::string tmpfile = "";
    //
    uint64 nsamples = 0;
    uint64 nsnps = 0;
    uint k = 10;
    uint maxp = 20; // maximum number of power iterations
    uint threads = 10;
    uint blocksize = 0;
    uint bands = 64;
    // for emu iteration
    uint maxiter = 100;
    double alpha = 0.001;
    // can be tol_emu or tol_pcangsd
    double tolem = 1e-4;
    double tolmaf = 1e-4;
    double maf = 0.0;
    // for arnoldi
    uint ncv = 20; // max(20, 2*k + 1)
    uint imaxiter = 1000;
    double itol = 1e-6;
    // for halko
    uint oversamples = 10;
    double tol = 1e-4;
    uint buffer = 2;

    double memory = 0; // 0 for disable
    bool cpmed = false;
    bool printv = false;
    bool runem = false;
    bool batch = true; // if load all matrix into RAM.
    bool noshuffle = false;
    bool fast = true;
    bool emu = false;
    bool pcangsd = false; // read GP field for PCAngsd instead of GT.
    bool halko = false;
    bool arnoldi = false;
    bool verbose = false;
    bool printu = false;
    bool mev = false;
    uint rand = 0;

    std::ostringstream ss;

};


#endif // CMD_H_
