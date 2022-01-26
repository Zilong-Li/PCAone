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
    virtual void read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) = 0;
    virtual void read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MatrixXd& U, const VectorXd& svals, const MatrixXd& VT, bool standardize = false) = 0;

    bool snpmajor = true;
    bool nsamples_ge_nsnps = false;  // if nsamples greater than or equal to nsnps
    bool initialFonly = false;
    uint64 nsamples, nsnps;  // prevent from potential integer overflow
    uint nblocks = 1;
    uint bandFactor = 1;
    std::vector<uint> start, stop;
    MatrixXd G;  // genotype matrix, can be initial E or centered E, which is nsamples x nsnps;
    MatrixXd P;  // genotype probability, nsnps x nsamples x 3.
    VectorXd F;  // observed or estimated population allele frequency
    VectorXd Dc; // diagnal vector of covariance matrix
    ArrayXXd centered_geno_lookup;
    std::vector<bool> C; // 1 or true indicates a ind's snp is missing and need to be predicted.
    const Param& params;


    void prepare(uint& blocksize);
    void read_bed_batch();
    void standardize_E();
    void pcangsd_standardize_E(const MatrixXd& U, const VectorXd& svals, const MatrixXd& VT);
    void update_batch_E(const MatrixXd& U, const VectorXd& svals, const MatrixXd& VT);
    void write_eigs_files(const VectorXd& S, const MatrixXd& U, const MatrixXd& V);

    // for blockwise
    void calcu_vt_initial(const MatrixXd& T, MatrixXd& VT);
    void calcu_vt_update(const MatrixXd& T, const MatrixXd& U, const VectorXd& svals, MatrixXd& VT, bool standardize);
    // update Eb, using V as predictor and Db as input
    // void update_block_E(uint start_idx, uint stop_idx, const MatrixXd& U, bool standardize = false);
    // MatrixXd calcu_vt_from_Eb(const MatrixXd& T, const MatrixXd& U, bool standardize);

    // calculate G * X or X * G by block
    // MatrixXd calcu_block_matmul(const MatrixXd& X, bool rightside);

    // MatrixXd calcu_block_matmul(const MatrixXd& X, bool rightside, const MatrixXd& U, const VectorXd& S, const MatrixXd& V, bool standardize = false);
    // calculate G' * X or X * G' by block
    // MatrixXd calcu_block_matmul_trans(const MatrixXd& X, bool rightside);

    // MatrixXd calcu_block_matmul_trans(const MatrixXd& X, bool rightside, const MatrixXd& U, const VectorXd& S, const MatrixXd& V, bool standardize = false);

};


#endif
