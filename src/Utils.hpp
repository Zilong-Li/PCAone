#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include "Common.hpp"
#include "Logger.hpp"
#include "Timer.hpp"
#include "zlib.h"

// MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
Logger cao;  // logger
Timer tick;  // Timer
#else
extern Timer tick;
extern Logger cao;
#endif

/// get the machine information
std::string get_machine();

void standardize(Mat2D& X, double tol = 1e-10);

/* @brief gets a line of gzfile
 * @param gz   file hander returned by gzopen
 * @param buf  buffer used for storing data
 * @param size buffer size for realloc buffer
 * @return extended buf length
 */
int tgets(gzFile gz, char** buf, uint64* size);

void fcloseOrDie(FILE* file);

FILE* fopenOrDie(const char* filename, const char* instruction);

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file);

size_t fwriteOrDie(const void* buffer, size_t sizeToWrite, FILE* file);

size_t count_lines(const std::string& fpath);

std::string timestamp();

void flip_UV(Mat2D& U, Mat2D& V, bool ubase = true);

void flip_Y(const Mat2D& X, Mat2D& Y);

double rmse(const Mat2D& X, const Mat2D& Y);

double rmse1d(const Mat1D& x, const Mat1D& y);

Mat1D minSSE(const Mat2D& X, const Mat2D& Y);

double mev(const Mat2D& X, const Mat2D& Y);

// out = X' * Y for a wide X (samples x sites) and a thin Y, so out is tall.
// Eigen's multithreaded GEMM packs a kc x X.cols() panel of X' for a tall
// result, up to min(n, 256) * m doubles: nearly a second copy of X when n is
// small. Column panels of X, one thread each, need a few MB and run 2-4x
// faster. Every row of out comes from one panel, over the whole of n, so the
// result does not depend on the number of threads.
void mul_Xt_Y(const Eigen::Ref<const Mat2D>& X, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out);

// K += X * X' on the lower triangle of the square K (the strict upper triangle
// is left untouched), for X of K.rows() x b. Eigen's rankUpdate runs on one
// thread and its GEMM does twice the flops, so the lower triangle is cut into
// tiles, each its own GEMM on one thread. A BLAS build hands it to ?syrk.
void syrk_lower_add(Eigen::Ref<Mat2D> K, const Eigen::Ref<const Mat2D>& X);

// copy the strict lower triangle of the square K onto its upper triangle
void mirror_lower(Eigen::Ref<Mat2D> K);

// Call f(j, x) for j in [0, count), where x is column `first + j` of A * B
// (A is n x k, B is k x m). The product is formed a panel of ~1 MB at a time in
// a thread-local buffer: one GEMM per panel in place of the scalar triple loops
// the reconstructions U * S * V' were written as. Runs in parallel over panels.
template <typename DerivedB, typename Func>
void for_each_product_column(
    const Mat2D& A, const Eigen::MatrixBase<DerivedB>& B, Eigen::Index first, Eigen::Index count, Func f) {
  const Eigen::Index bs = std::max<Eigen::Index>(1, (Eigen::Index(1) << 17) / std::max<Eigen::Index>(1, A.rows()));
  const Eigen::Index nb = (count + bs - 1) / bs;
#pragma omp parallel
  {
    Mat2D panel;
#pragma omp for schedule(static)
    for (Eigen::Index b = 0; b < nb; ++b) {
      const Eigen::Index c = b * bs, w = std::min<Eigen::Index>(bs, count - c);
      panel.noalias() = A * B.middleCols(first + c, w);
      for (Eigen::Index t = 0; t < w; ++t) f(c + t, panel.col(t));
    }
  }
}

void mev_rmse_byk(const Mat2D& X, const Mat2D& Y, Mat1D& Vm, Mat1D& Vr);

String1D split_string(const std::string& s, const std::string& separators);

// the EM-PCA loops ran out of --maxiter without a word
template <typename P>
void warn_em_not_converged(const P& params, double diff) {
  cao.warn("EM-PCA did not converge in --maxiter " + std::to_string(params.maxiter) +
           " iterations (diff = " + std::to_string(diff) + ", --tol-em = " + std::to_string(params.tolem) + ")");
}

// Median of a buffer the caller is done with: reorders it in place.
//
// std::nth_element is O(n) where std::sort is O(n log n), and the value is
// identical -- after nth_element every element below `mid` is <= v[mid], so the
// largest of them is the other central order statistic. robust_cov_gk() calls
// this millions of times, so the difference is the runtime of --selection 2.
template <typename T>
auto median_inplace(std::vector<T>& v) {
  static_assert(!std::is_same_v<T, bool>, "Boolean type is not supported");

  if (v.empty()) {
    throw std::invalid_argument("Cannot calculate median of an empty vector");
  }

  const size_t n = v.size();
  const size_t mid = n / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  if (n % 2 != 0) return v[mid];
  const T lo = *std::max_element(v.begin(), v.begin() + mid);
  return static_cast<T>((lo + v[mid]) / static_cast<T>(2));
}

template <typename T>
auto get_median(std::vector<T> v) {
  return median_inplace(v);
}

// Write a matrix as tab-separated rows, NA for a non-finite value. Finite
// values come out byte for byte as Eigen::IOFormat(6, DontAlignCols, "\t", "\n")
// prints them; that prints a NaN as "nan", which R does not read as missing.
template <typename Derived>
void write_rows(std::ostream& os, const Eigen::DenseBase<Derived>& M) {
  const std::streamsize old = os.precision(6);
  for (Eigen::Index i = 0; i < M.rows(); ++i) {
    for (Eigen::Index j = 0; j < M.cols(); ++j) {
      if (j) os << '\t';
      const double v = M(i, j);
      if (std::isfinite(v))
        os << v;
      else
        os << "NA";
    }
    os << '\n';
  }
  os.precision(old);
}

void make_plink2_eigenvec_file(int K, std::string fout, const std::string& fin, const std::string& fam);
void make_plink2_eigenvec_from_psam(int K, const std::string& fout, const std::string& fin, const std::string& fpsam);

bool isZstdCompressed(const char* filename);

Mat2D read_usv(const std::string& path);

// How the matrix that a PCA decomposed relates to the 0..1 allele-frequency
// scale, recorded in .sigvals so that -P/--USV can invert U*S*V' later.
//
// Nothing in .eigvecs/.sigvals/.loadings used to say this, so --inbreed had to
// assume the reference run used the defaults. It does not always: --missme
// without --emu leaves the final decomposition unstandardised, and the pcangsd
// path decomposes centred *dosages* rather than PCAone's 0..1 coding.
struct UsvTransform {
  bool known = false;                     // false for a .sigvals written before this was recorded
  int scale = SCALE_STANDARDIZE_GENETIC;  // scaling actually applied: -9 standardised, 0 none
  int ploidy = 2;                         // ploidy of the reference run
  int gscale = 1;                         // 1: genotypes coded 0..1; 2: dosages coded 0..2
};

void read_sigvals(const std::string& path, uint& N, uint& M, Mat1D& S, UsvTransform* transform = nullptr);

Mat1D read_eigvals(const std::string& path);

Mat2D read_eigvecs(const std::string& path, int n, int k);

Mat1D read_frq(const std::string& path);

struct BimMatch {
  bool identical = false;
  Int1D bim_indices;
  Int1D mbim_indices;
  std::vector<bool> flip;  // parallel to bim_indices: true if ref/alt are swapped vs mbim
};

BimMatch match_bim_to_mbim(const std::string& bim_file, const std::string& mbim_file);

BimMatch match_pvar_to_mbim(const std::string& pvar_file, const std::string& mbim_file);

BimMatch match_beagle_to_mbim(const std::string& beagle_file, const std::string& mbim_file);

std::string pvar_line_to_bim_line(const std::string& line, const std::string& path);

std::string decode_beagle_allele(const std::string& allele);

void parse_beagle_file(Mat2D& P, gzFile fp, const int nsamples, const int nsnps);

String1D parse_beagle_samples(const std::string& fin);

void write_eigvecs2_beagle(const Mat2D& U, const std::string& fin, const std::string& fout);

// .cov and .eigvecs2 of the PCAngsd path: the covariance of the standardized
// expected genotypes E with its corrected diagonal Dc, and its top -k eigenvectors
class Param;
void write_pcangsd_cov(const Mat2D& E, const Mat1D& Dc, uint nsnps, const Param& params);

/// return the p-value of 1-degreed chi-squared
double chisq1d(const double x);

/// stream compress a file by zstd
void zstd_compress_file(const std::string& fname, std::string outname, int level);

/// estimate population allele frequency with genotype likelihoods
void emMAF_with_GL(Mat1D& F, const Mat2D& P, int maxiter, double tolmaf);

// zstd deccompression buffer
struct ZstdDS {
  ZstdDS() {
    buffInTmp.reserve(buffInSize);
    buffOutTmp.reserve(buffOutSize);
  }
  ~ZstdDS() {
    ZSTD_freeDCtx(dctx);
    fcloseOrDie(fin);
  }
  FILE* fin = nullptr;
  size_t const buffInSize = ZSTD_DStreamInSize();
  size_t const buffOutSize = ZSTD_DStreamOutSize();
  ZSTD_DCtx* const dctx = ZSTD_createDCtx();
  size_t lastRet = 1;
  std::string buffInTmp, buffOutTmp;
};

// zstd compression buffer
struct ZstdCS {
  ZstdCS() {
    buffInTmp.reserve(buffInSize);
    buffOutTmp.reserve(buffOutSize);
  }
  ~ZstdCS() {
    ZSTD_freeCCtx(cctx);
    fcloseOrDie(fout);
  }
  FILE* fout = nullptr;
  size_t const buffInSize = ZSTD_CStreamInSize();
  size_t const buffOutSize = ZSTD_CStreamOutSize();
  ZSTD_CCtx* const cctx = ZSTD_createCCtx();
  size_t lastRet = 1;
  std::string buffInTmp, buffOutTmp;
};

// Chi-square CDF: P(X <= x) where X ~ chi-sq(df)
// chi-sq CDF = regularized lower incomplete gamma P(df/2, x/2)
// pchisq(x, df, lower_tail=true)  = kf_gammap(df/2, x/2)
// pchisq(x, df, lower_tail=false) = kf_gammaq(df/2, x/2)
double pchisq(double x, int df, bool lower_tail);

// Chi-square quantile function via Newton's method
// Find x such that P(X <= x) = p, where X ~ chi-sq(df)
// Uses the chi-square PDF for the Newton update:
//   f(x) = x^(df/2-1) * exp(-x/2) / (2^(df/2) * Gamma(df/2))
double qchisq(double p, int df);

void galinsky_selection_stat(Mat2D& V);
// `keep[i] == 0` marks a site with no residual variance: it is left out of the
// robust fit and the inflation factor, and its outputs come back as NaN.
void pcadapt_selection_stats(const Mat2D& Z, const std::vector<char>& keep, Mat1D& stat, Mat1D& chi2_stat,
                             Mat1D& pval, double& gif);

#endif  // PCAONE_UTILES_
