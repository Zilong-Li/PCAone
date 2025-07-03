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

template <typename T>
auto get_median(std::vector<T> v) {
    static_assert(!std::is_same_v<T, bool>, "Boolean type is not supported");

    if (v.empty()) {
        throw std::invalid_argument("Cannot calculate median of an empty vector");
    }

    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2 == 0) {
        return (v[n / 2 - 1] + v[n / 2]) / static_cast<T>(2);
    } else {
        return v[n / 2];
    }
}

void make_plink2_eigenvec_file(int K, std::string fout, const std::string& fin, const std::string& fam);

bool isZstdCompressed(const char* filename);

Mat2D read_usv(const std::string& path);

void read_sigvals(const std::string& path, uint& N, uint& M, Mat1D& S);

Mat1D read_eigvals(const std::string& path);

Mat2D read_eigvecs(const std::string& path, int n, int k);

Mat1D read_frq(const std::string& path);

void check_bim_vs_mbim(const std::string& bim_file, const std::string& mbim_file);

void parse_beagle_file(Mat2D& P, gzFile fp, const int nsamples, const int nsnps);

String1D parse_beagle_samples(const std::string& fin);

void write_eigvecs2_beagle(const Mat2D& U, const std::string& fin, const std::string& fout);

/// return the p-value of 1-degreed chi-squared
double chisq1d(const double x);

/// stream compress a file by zstd
void zstd_compress_file(const std::string& fname, std::string outname, int level);

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

#endif  // PCAONE_UTILES_
