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

void mev_rmse_byk(const Mat2D& X, const Mat2D& Y, Mat1D& Vm, Mat1D& Vr);

String1D split_string(const std::string& s, const std::string& separators);

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
