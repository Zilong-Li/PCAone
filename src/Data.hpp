#ifndef PCAONE_DATA_
#define PCAONE_DATA_

#include "Cmd.hpp"
#include "Common.hpp"
#include "PackedPca.hpp"
#include "TextIO.hpp"  // UsvTransform
#include "Utils.hpp"

const double VAR_TOL = 1e-9;

// the 2-bit codes of read_block_packed(): the allele frequency (0..1) each
// stands for, and the code of a missing call
struct PackedCoding {
  double value[4];
  int missing;
};

class Data {
 public:
  Data(const Param& params_)
      : params(params_) {}

  virtual ~Data() = default;

  virtual void read_all() = 0;
  // for blockwise
  virtual void check_file_offset_first_var() = 0;
  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) = 0;
  virtual void read_block_update(
      uint64 start, uint64 stop, const Mat2D& U, const Mat1D& svals, const Mat2D& VT, bool standardize) = 0;

  // out-of-core --project: the 2-bit calls of the sites [start_idx, stop_idx],
  // a record of (nsamples + 3) / 4 bytes each, 4 samples a byte, the first in
  // the low bits, valid until the next read. The caller centres and scales
  // them, by packed_coding() and F. Only if packed_calls(), else the blocks are
  // read by read_block_initial().
  virtual bool packed_calls() { return false; }
  virtual const uint8_t* read_block_packed(uint64 start_idx, uint64 stop_idx) { return nullptr; }
  virtual PackedCoding packed_coding() const { return {}; }

  void prepare();
  // blocks of b sites over the nsnps sites, for a caller that sizes its own
  // blocks: out-of-core --project (run_projection_blockwise)
  void set_blocks(uint b);
  // one warning for all sites with MAF=0. The readers count them inside their
  // OpenMP loops and call this afterwards, because cao is not thread-safe.
  void warn_monomorphic(uint64 n) const {
    if (params.maf == 0 && n > 0) cao.warn(std::to_string(n) + " sites with MAF=0 found! remove them first!");
  }
  void standardize_E();
  // whether the RSVD of this run works on the 2-bit calls (packed): a PCA by
  // --svd 1 or 2 of a BED, centred, and standardized or not, with no EM and no
  // --no-packed. In-core, read_all() keeps them; out-of-core, the blocks are
  // read by read_block_packed().
  bool packed_pca() const;
  // the centring and scaling of the sites [j0, j0 + n) of packed (all from j0
  // by default), as read_all() and standardize_E() give G, with the missing
  // calls at the mean. Out-of-core, packed is made here, without records.
  void set_packed_sites(bool standardize, uint64 j0 = 0, uint64 n = (uint64)-1);
  void filter_snps_resize_F();  // filter first, then update nsnps
  void save_snps_in_mbim();
  void pcangsd_standardize_E(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  void fit_with_pi(const Mat2D& U, const Mat1D& svals, const Mat2D& VT);
  // In-core winSVD keeps G and V in shuffled order, but F and C retain
  // filtered input order. P additionally needs the keepSNPs mapping.
  uint unpermuted_snp_index(uint j) const { return in_core_permuted ? static_cast<uint>(perm.indices()(j)) : j; }
  // .sigvals, .eigvals, .eigvecs, .eigvecs2 (PLINK, PGEN) and, unless
  // --no-loadings, .loadings and .mbim
  void write_eigs_files(const Mat1D& E, const Mat1D& S, const Mat2D& U, const Mat2D& V);
  // .cov and .eigvecs2 of the PCAngsd path: the covariance of the standardized
  // expected genotypes G with its corrected diagonal Dc, and its top -k eigenvectors
  void write_pcangsd_cov() const;
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
  // VTout = T * G, where G is imputed from U, svals and VT as it is read. VTout
  // may be VT itself: each block of VT is read before it is overwritten.
  void calcu_vt_update(
      const Mat2D& T, const Mat2D& U, const Mat1D& svals, const Mat2D& VT, Mat2D& VTout, bool standardize);
  // winEM's E-step (--winem): impute the missing calls of the columns of Z, the
  // sites [j0, j0 + Z.cols()) of G, by least squares of each site's observed
  // calls on the rows of U (nsamples x k, orthonormal) that observe them, then
  // clipped to [-f, 1 - f] and refined by EMU's step until its fixed point.
  // Stateless: it reads only U and the observed calls, not the old imputations.
  // In-core the missing calls are those of C; out-of-core they are read as NaN
  // (mark_missing). With an empty U the missing calls are set to the mean, 0.
  void estep_reg(const Mat2D& U, Eigen::Ref<Mat2D> Z, uint64 j0);
  // winEM's E-step for --pcangsd (in-core): the expected genotypes of the
  // columns [j0, j0 + Z.cols()) of G, given their genotype likelihoods P and
  // the individual allele frequencies of the fit U U'e of each site's own
  // expected genotypes e. That fit depends on e, so each site iterates
  // e = E[g | P, pi(U U'e)] - 2f from its expectations in G towards its fixed
  // point under U, a few steps per pass, where the EM of --pcangsd takes one.
  void estep_pcangsd(const Mat2D& U, Eigen::Ref<Mat2D> Z, uint64 j0);
  // estep_pcangsd() for --pcangsd, else estep_reg()
  void estep(const Mat2D& U, Eigen::Ref<Mat2D> Z, uint64 j0) {
    if (params.pcangsd)
      estep_pcangsd(U, Z, j0);
    else
      estep_reg(U, Z, j0);
  }
  // given PCs, predict the missing values, then update in place by block
  void predict_missing_E(const Mat2D& U, uint64 start_idx, uint64 stop_idx);

 public:
  const Param& params;
  double readtime = 0;
  double esteptime = 0;  // seconds in estep_reg(), winEM's E-step
  uint64 estep_unconverged = 0;  // sites whose fixed point estep_reg() did not reach
  // winEM's over-relaxation: the missing information of each sample, as a
  // fraction of its sites. EMU: the missing calls, counted by estep_reg() in
  // the first pass (in-core from C, out-of-core from the NaN). PCAngsd: the
  // sum over sites of Var(g | GL, pi) / (2 pi (1 - pi)) of the last pass,
  // accumulated by estep_pcangsd() and read by rbar_pass().
  Mat1D miss_count, rbar_sum;
  bool miss_counted = false;
  // the change of the expectations over the pass, sum of (e_new - e_old)^2,
  // and their size, sum of e_new^2: PCAngsd's E-step lags behind the basis,
  // so winEM also tests their relative change per pass before it stops
  double de_sum = 0, e_sum = 0;
  void reset_pass_sums() {
    rbar_sum.setZero();
    de_sum = e_sum = 0;
  }
  Mat1D rbar_pass() const { return rbar_sum / (double)nsnps; }
  bool snpmajor = true;
  bool nsamples_ge_nsnps = false;  // if nsamples greater than or equal to nsnps
  uint blocksize = 0, nsamples = 0, nsnps = 0;
  uint nsnps_all = 0;  // sites in the input before --maf; keepSNPs indexes them
  uint nblocks = 1;
  uint bandFactor = 1;
  uint nops = 0;
  std::vector<uint> start, stop;
  double p_miss = 0.0;             // proportion of genotype missingness
  PermMat perm;                    // permuation order of SNPs
  bool in_core_permuted = false;   // true only after G has actually been shuffled
  bool metadata_permuted = false;  // the .bim read is the permuted copy's (BED rewritten out-of-core)
  // --svd 1 and 2 of a BED: the 2-bit calls in place of G (PackedPca.hpp), its
  // sites in the order of G's columns. In-core read_all() makes it when
  // packed_pca(); out-of-core set_packed_sites(), without the calls.
  std::unique_ptr<PCAone::PackedGenotypes> packed;
  Mat2D G;                         // genotype matrix, can be initial E or centered E, which is nsamples x nsnps;
  Mat2D P;                         // normalized genotype likelihoods, (nsamples x 2) x nsnps.
  Mat1D F;                         // observed or estimated population allele frequency
  Mat1D Dc;                        // diagnal vector of covariance matrix
  ArrBool C;                       // nsnps x nsample, if there is missing value
  // --winem in-core: the missing samples of each column of G, in G's (shuffled)
  // order: those of column j are miss_idx[miss_ptr[j] .. miss_ptr[j + 1])
  std::vector<uint64> miss_ptr;
  std::vector<uint> miss_idx;
  Arr2D centered_geno_lookup;      // lookup table for centering genotypes
  // out-of-core --project: read_block_initial() returns a missing call as NaN
  // instead of 0 (the reference mean), so the caller can tell it from a called
  // genotype. Only with frequencies taken from the reference (!params.dopca).
  bool mark_missing = false;
  int svd_scale = SCALE_STANDARDIZE_GENETIC;  // scaling applied to the decomposed matrix; 0 = none
  int svd_gscale = 1;                         // 1: genotypes coded 0..1; 2: dosages coded 0..2
  Int1D keepSNPs;                             // store index of SNPs to keep
  Int1D keepRefSNPs;                          // store matching SNP indices in the reference .mbim
  Int1D flipSNPs;                             // local SNP indices (into keepSNPs order) with flipped alleles
};

#endif  // PCAONE_DATA_
