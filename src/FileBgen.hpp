#ifndef __FileBgen__
#define __FileBgen__

#include "Data.hpp"
#include "bgen/bgen.h"

// const double GENOTYPE_THRESHOLD = 0.9;
// const float BGEN_MISSING_VALUE = -9;
// const float BGEN2GENO[4] = {0, 0.5, 1, BGEN_MISSING_VALUE};

using namespace bgen;

class FileBgen : public Data
{
public:
    // use Data default constructor;
    using Data::Data;

    ~FileBgen() {}

    virtual void get_matrix_dimensions();
    virtual void read_all_and_centering();
    // for blockwise
    // virtual void estimate_F() {}
    virtual void read_snp_block_initial(uint start_idx, uint stop_idx, bool standardize = false) {}
    virtual void read_snp_block_update(uint start_idx, uint stop_idx, const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT, bool standardize = false) {}


    virtual void open_check_file() {}
    virtual void close_check_file() {}

private:
    Bgen* bg;

};



#endif