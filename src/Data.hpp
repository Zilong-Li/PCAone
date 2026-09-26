#ifndef PCAONE_DATA_
#define PCAONE_DATA_

#include "Cmd.hpp"
#include "Common.hpp"
#include "Utils.hpp"  // UsvTransform

const double VAR_TOL = 1e-9;

class Data {
 public:
  Data(const Param& params_)
      : params(params_) {}

  virtual ~Data() = default;

  virtual void read_all() = 0;
  // for blockwise
  virtual void check_file_offset_first_var() = 0;
  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) = 0;
  virtual void read_block_update(uint64 start, uint64 stop, const Mat2D& U, const Mat1D& svals, const Mat2D& VT, bool standardize) = 0;

  void prepare();
  // one warning for all sites with MAF=0. The readers count them inside their
  // OpenMP loops and call this afterwards, because cao is not thread-safe.
  void warn_monomorphic(uint64 n) const {
    if (params.maf == 0 && n > 0) cao.warn(std::to_string(n) + " sites with MAF=0 found! remove them first!");
  }
  void standardize_E();
  void filter_snps_resize_F();  // filter first, then update nsnps
  void save_snps_in_mbim();
  void pcangsd_standardize_E(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  void fit_with_pi(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  // In-core winSVD keeps G and V in shuffled order, but F and C retain
  // filtered input order. P additionally needs the keepSNPs mapping.
  uint unpermuted_snp_index(uint j) const {
    return in_core_permuted ? static_cast<uint>(perm.indices()(j)) : j;
  }
  void write_eigs_files(const Mat1D& E, const Mat1D& S, const Mat2D& U, const Mat2D& V);
  // Record how the matrix that was just decomposed relates to the 0..1
  // allele-frequency scale, so write_eigs_files() can put it in .sigvals and
  // -P/--USV can invert U*S*V' later. Call right before write_eigs_files(),
  // passing whether the FINAL decomposition standardized.
  void set_svd_transform(bool standardized);
  // Projection and selection read U/S/V from a reference PCA and compare them
  // against a genotype matrix they scale themselves, so that matrix has to
  // carry the transform the REFERENCE applied -- not the one this run's
  // -C/--scale asks for. resolve_ref_scaling() says whether to standardize and
  // rejects transforms that cannot be replayed; the other two apply it, in
  // place of standardize_E(), which is gated on this run's --scale.
  // allow_dosage accepts a reference that decomposed 0..2 dosages (gscale=2);
  // only --project 3 can replay that, since it rescales the target per site.
  bool resolve_ref_scaling(const UsvTransform& t, const std::string& src, bool allow_dosage = false) const;
  void standardize_E_ref(const UsvTransform& t);
  void standardize_block_ref(const UsvTransform& t, uint64 start_idx, uint block_cols);
  // for blockwise
  // void fit_with_pi_block(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  void calcu_vt_initial(const Mat2D& T, Mat2D& VT, bool standardize);
  void calcu_vt_update(const Mat2D& T, const Mat2D& U, const Mat1D& svals, Mat2D& VT, bool standardize);
  // given PCs, predict the missing values, then update in place by block
  void predict_missing_E(const Mat2D& U, uint64 start_idx, uint64 stop_idx);

 public:
  const Param& params;
  double readtime = 0;
  bool snpmajor = true;
  bool nsamples_ge_nsnps = false;  // if nsamples greater than or equal to nsnps
  uint blocksize = 0, nsamples = 0, nsnps = 0;
  uint nsnps_all = 0;  // sites in the input before --maf; keepSNPs indexes them
  uint nblocks = 1;
  uint bandFactor = 1;
  uint nops = 0;
  std::vector<uint> start, stop;
  double p_miss = 0.0;         // proportion of genotype missingness
  PermMat perm;                // permuation order of SNPs
  bool in_core_permuted = false;  // true only after G has actually been shuffled
  Mat2D G;                     // genotype matrix, can be initial E or centered E, which is nsamples x nsnps;
  Mat2D P;                     // normalized genotype likelihoods, (nsamples x 2) x nsnps.
  Mat1D F;                     // observed or estimated population allele frequency
  Mat1D Dc;                    // diagnal vector of covariance matrix
  ArrBool C;                   // nsnps x nsample, if there is missing value
  Arr2D centered_geno_lookup;  // lookup table for centering genotypes
  int svd_scale = SCALE_STANDARDIZE_GENETIC;  // scaling applied to the decomposed matrix; 0 = none
  int svd_gscale = 1;                         // 1: genotypes coded 0..1; 2: dosages coded 0..2
  Int1D keepSNPs;              // store index of SNPs to keep
  Int1D keepRefSNPs;           // store matching SNP indices in the reference .mbim
  Int1D flipSNPs;              // local SNP indices (into keepSNPs order) with flipped alleles
};

#endif  // PCAONE_DATA_
