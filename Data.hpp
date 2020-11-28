#ifndef __DATA_H__
#define __DATA_H__

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <vector>
#include "Utils.hpp" // count_lines
#include "EMU.hpp"

using namespace Eigen;
using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

/**
* Recode genotype codes to allelic dosages of first allele in .bim file (A1),
* similarly to .raw files generated with '--recode A' in PLINK. A coding for
* the missing value needs to be provided in 'na_value'.
* 00 ->  2 (copies of A1)
* 10 ->  1 (copy of A1)
* 11 ->  0 (copy of A1)
* 01 -> -9 (missing)
*/
#define MISS_INTEGER -9
const int MAP2GENO[4] = {2, MISS_INTEGER, 1, 0};

class Data
{
public:
    Data(const Param& params_): params(params_) {}

    ~Data() {}

    const Param& params;

    void run();
    void read_bed_whole_genomat();
    void init_whole_E();
    void run_whole_EM();
    void update_whole_E();

    size_t nsamples, nsnps;
    MatrixXd G;
    VectorXd F;

private:
    vector<uchar> inbed;
    vector<bool> C;
};


#endif