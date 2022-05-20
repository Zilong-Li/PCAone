#ifndef __FileBgen__
#define __FileBgen__

#include "Data.hpp"
#include "bgen/bgen.h"

// const double GENOTYPE_THRESHOLD = 0.9;
// const double BGEN_MISSING_VALUE = -9;
// const double BGEN2GENO[4] = {0, 0.5, 1, BGEN_MISSING_VALUE};

class FileBgen : public Data
{
public:
    // using Data::Data;
    FileBgen(const Param& params_) : Data(params_)
        {
            llog << timestamp() << "start parsing BGEN format" << std::endl;
            bg = new bgen::Bgen(params.bgen, "", true);
            nsamples = bg->header.nsamples;
            nsnps = bg->header.nvariants;
            llog << timestamp() << "the layout of bgen file is " << bg->header.layout << ". N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
        }

    ~FileBgen() { delete bg; }

    virtual void read_all_and_centering();
    // for blockwise
    virtual void check_file_offset_first_var() { bg->set_offset_first_var(); }

    virtual void read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);

    virtual void read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize = false) {}

private:
    bgen::Bgen* bg;
    bgen::Variant var;
    float* dosages;
    float* probs1d;
    bool frequency_was_estimated = false;

};



#endif
