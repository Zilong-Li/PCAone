#ifndef __FilePlink__
#define __FilePlink__

#include "Data.hpp"

/**
* Recode genotype codes to allelic dosages of first allele in .bim file (A1),
* similarly to .raw files generated with '--recode A' in PLINK. A coding for
* the missing value needs to be provided in 'na_value'.
* 00 ->  2 (copies of A1)
* 10 ->  1 (copy of A1)
* 11 ->  0 (copy of A1)
* 01 ->  3 (missing)
*/
const float BED_MISSING_VALUE = -9;
const float BED2GENO[4] = {1, BED_MISSING_VALUE, 0.5, 0};

class FileBed : public Data
{
public:
    //
    using Data::Data;

    ~FileBed() {}

    virtual void get_matrix_dimensions();
    virtual void read_all_and_centering();
    // for blockwise
    virtual void read_snp_block_initial(uint start_idx, uint stop_idx, bool standardize = false);
    virtual void read_snp_block_update(uint start_idx, uint stop_idx, const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT, bool standardize = false);


    virtual void open_check_file();
    virtual void close_check_file();

private:
    std::ifstream bed_ifstream;
    uint64 bed_bytes_per_snp;
    bool frequency_was_estimated = false;
    vector<uchar> inbed;
};

// bgen repo: https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk

#endif