#ifndef PCAONE_DATA_
#define PCAONE_DATA_

#include "Cmd.hpp"
#include "Utils.hpp"

const double VAR_TOL = 1e-9;

class Data {
 public:
  Data(const Param& params_) : params(params_) {}

  virtual ~Data() = default;

  virtual void read_all() = 0;
  // for blockwise
  virtual void check_file_offset_first_var() = 0;
  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) = 0;
  virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals,
                                 const Mat2D& VT, bool standardize) = 0;

  void prepare();
  void standardize_E();
  void filter_snps_resize_F();
  void save_snps_in_bim();
  void pcangsd_standardize_E(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  void update_batch_E(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  void write_eigs_files(const Mat1D& E, const Mat1D& S, const Mat2D& U, const Mat2D& V);
  void write_residuals(const Mat1D& S, const Mat2D& U, const Mat2D& VT);
  // for blockwise
  void calcu_vt_initial(const Mat2D& T, Mat2D& VT, bool standardize);
  void calcu_vt_update(const Mat2D& T, const Mat2D& U, const Mat1D& svals, Mat2D& VT, bool standardize);
  // given PCs, predict the missing values, then update in place by block
  void predict_missing_E(const Mat2D& U,uint64 start_idx, uint64 stop_idx);

  const Param& params;

  double readtime = 0;
  bool snpmajor = true;
  bool nsamples_ge_nsnps = false;  // if nsamples greater than or equal to nsnps
  bool impute = false;  // true if emu is enabled and there are missing values
  uint blocksize = 0, nsamples = 0, nsnps = 0;
  uint nblocks = 1;
  uint bandFactor = 1;
  uint nops = 0;
  std::vector<uint> start, stop;
  PermMat perm;  // permuation order of SNPs
  Mat2D G;       // genotype matrix, can be initial E or centered E, which is nsamples x nsnps;
  Mat2D P;       // normalized genotype likelihoods, (nsamples x 2) x nsnps.
  Mat1D F;       // observed or estimated population allele frequency
  Mat1D Dc;      // diagnal vector of covariance matrix
  ArrBool C;     // nsnps x nsample, if there is missing value
  Arr2D centered_geno_lookup;
  std::vector<int> keepSNPs;  // store index of SNPs to keep
};

#endif  // PCAONE_DATA_
