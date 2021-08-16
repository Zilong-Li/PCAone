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
    // using Data::Data;
    FileBgen(const Param& params_) : Data(params_)
        {
            bg = new Bgen(params.bgen, "", true);
            nsamples = bg->header.nsamples;
            nsnps = bg->header.nvariants;
            cout << timestamp() << "the layout of bgen file is " << bg->header.layout << ". N samples is " << nsamples << ". M snps is " << nsnps << endl;
        }

    ~FileBgen() { delete bg; }

    virtual void read_all_and_centering();
    // for blockwise
    // virtual void estimate_F() {}
    virtual void read_snp_block_initial(uint start_idx, uint stop_idx, bool standardize = false);
    virtual void read_snp_block_update(uint start_idx, uint stop_idx, const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT, bool standardize = false) {}


    virtual void check_file_offset_first_var() { bg->set_offset_first_var(); }

private:
    Bgen* bg;
    Variant var;
    float* dosages;
    float* probs1d;
    bool frequency_was_estimated = false;

};



#endif