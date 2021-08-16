#ifndef __DATA_H__
#define __DATA_H__

#include "Utils.hpp"

class Data
{
public:
    Data(const Param& params_): params(params_) {}

    virtual ~Data() {}

    virtual void read_all_and_centering() = 0;
    // for blockwise
    virtual void check_file_offset_first_var() = 0;
    virtual void read_snp_block_initial(uint start_idx, uint stop_idx, bool standardize = false) = 0;
    virtual void read_snp_block_update(uint start_idx, uint stop_idx, const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT, bool standardize = false) = 0;

    bool snpmajor = true;
    bool nsamples_ge_nsnps = false;  // if nsamples greater than or equal to nsnps
    bool initialFonly = false;
    uint nsamples, nsnps;
    uint nblocks = 1;
    uint bandFactor = 1;
    vector<uint> start, stop;
    MatrixXf G;  // genotype matrix, can be initial E or centered E, which is nsamples x nsnps;
    MatrixXf P;  // genotype probability, nsnps x nsamples x 3.
    VectorXf F;  // observed or estimated population allele frequency
    VectorXf Dc; // diagnal vector of covariance matrix
    ArrayXXf centered_geno_lookup;
    vector<bool> C; // 1 or true indicates a ind's snp is missing and need to be predicted.
    const Param& params;


    void prepare(uint& blocksize);
    void read_bed_batch();
    void standardize_E();
    void pcangsd_standardize_E(const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT);
    void update_batch_E(const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT);
    void write_eigs_files(const VectorXf& vals, const MatrixXf& vecs);

    // for blockwise
    void calcu_vt_initial(const MatrixXf& T, MatrixXf& VT);
    void calcu_vt_update(const MatrixXf& T, const MatrixXf& U, const VectorXf& svals, MatrixXf& VT, bool standardize);
    // update Eb, using V as predictor and Db as input
    // void update_block_E(uint start_idx, uint stop_idx, const MatrixXf& U, bool standardize = false);
    // MatrixXf calcu_vt_from_Eb(const MatrixXf& T, const MatrixXf& U, bool standardize);

    // calculate G * X or X * G by block
    // MatrixXf calcu_block_matmul(const MatrixXf& X, bool rightside);

    // MatrixXf calcu_block_matmul(const MatrixXf& X, bool rightside, const MatrixXf& U, const VectorXf& S, const MatrixXf& V, bool standardize = false);
    // calculate G' * X or X * G' by block
    // MatrixXf calcu_block_matmul_trans(const MatrixXf& X, bool rightside);

    // MatrixXf calcu_block_matmul_trans(const MatrixXf& X, bool rightside, const MatrixXf& U, const VectorXf& S, const MatrixXf& V, bool standardize = false);

};


#endif