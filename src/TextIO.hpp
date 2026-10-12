/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/TextIO.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 *
 * Text input and output: the files of a PCA (.eigvecs, .eigvals, .sigvals,
 * .loadings, .mbim), the sample IDs of the input, and the fast, parallel
 * writer that every large text output goes through.
 ******************************************************************************/
#ifndef PCAONE_TEXTIO_
#define PCAONE_TEXTIO_

#include <omp.h>

#include <cstring>
#include <memory>
#include <ostream>
#include <thread>

#include "Common.hpp"

class Param;

// ---------------------------------------------------------------------------
// numbers to text
// ---------------------------------------------------------------------------
//
// An ostream formats a double through the locale, ~230 ns a value, and Eigen's
// operator<< adds a stringstream per coefficient when it aligns the columns.
// These write the same bytes, rounded exactly as printf rounds (to nearest,
// ties to even, from the exact binary value), in ~35 ns, and need no locale,
// so they run on every thread. tests/test_textio.cpp checks them against
// printf and Eigen.
//
// Each returns the end of what it wrote at p, and may write scratch beyond it:
// TEXT_NUM_MAX bytes at p must be writable. A NaN is written "nan" whatever
// its sign bit (glibc's printf writes some as "-nan").

constexpr size_t TEXT_NUM_MAX = 24;

// printf("%.6g"): what an ostream writes at its default precision, and so what
// Eigen::IOFormat(6) and `os << x` wrote
char* put_g6(double x, char* p);

// printf("%.6f"): what std::to_string writes. |x| >= 1e15 is written as
// "%.6g", to bound the length
char* put_f6(double x, char* p);

char* put_uint(uint64 v, char* p);

// A byte buffer formatters write into through a raw pointer:
//   char* p = buf.reserve(maxlen);  ...write up to maxlen bytes...  buf.commit(p);
class TextBuf {
 public:
  char* reserve(size_t n) {
    if (len_ + n > cap_) grow(len_ + n);
    return b_.get() + len_;
  }
  void commit(const char* end) { len_ = (size_t)(end - b_.get()); }
  void append(const char* s, size_t n) {
    std::memcpy(reserve(n), s, n);
    len_ += n;
  }
  void append(const std::string& s) { append(s.data(), s.size()); }
  void clear() { len_ = 0; }
  size_t size() const { return len_; }
  const char* data() const { return b_.get(); }
  void write_to(std::ostream& os) const { os.write(b_.get(), (std::streamsize)len_); }

 private:
  void grow(size_t need) {  // not zero-filled, unlike a std::vector
    const size_t cap = std::max(need, 2 * cap_);
    std::unique_ptr<char[]> b(new char[cap]);
    if (len_) std::memcpy(b.get(), b_.get(), len_);
    b_.swap(b);
    cap_ = cap;
  }
  std::unique_ptr<char[]> b_;
  size_t len_ = 0, cap_ = 0;
};

// Write lines 0..n-1 to os, in order, each appended to a TextBuf by
// line(i, buf), '\n' included. `line` must not throw. The lines are formatted
// on every thread, a block of about TEXT_BLOCK_BYTES per thread at a time
// (bytes_per_line is a guess that sizes the blocks), and each block is written
// while the next one is formatted. The trick of EvalAdmix's write_matrix,
// generalised.
constexpr size_t TEXT_BLOCK_BYTES = size_t(4) << 20;  // fits in L2/L3, so the bytes are written from cache
template <typename LineFn>
void write_lines(std::ostream& os, Eigen::Index n, size_t bytes_per_line, LineFn&& line) {
  if (n <= 0) return;
  const int T = omp_in_parallel() ? 1 : std::max(1, omp_get_max_threads());
  const Eigen::Index want = (Eigen::Index)(T * TEXT_BLOCK_BYTES / std::max<size_t>(1, bytes_per_line));
  const Eigen::Index B = std::min<Eigen::Index>(n, std::max<Eigen::Index>(T, want));
  std::vector<TextBuf> cur(T), out(T);
  std::thread writer;
  for (Eigen::Index b0 = 0; b0 < n; b0 += B) {
    const Eigen::Index h = std::min(B, n - b0);
    for (auto& b : cur) b.clear();
#pragma omp parallel num_threads(T)
    {
      const Eigen::Index t = omp_get_thread_num(), nt = omp_get_num_threads();
      TextBuf& buf = cur[t];
      for (Eigen::Index i = b0 + h * t / nt; i < b0 + h * (t + 1) / nt; ++i) line(i, buf);
    }
    if (writer.joinable()) writer.join();
    cur.swap(out);
    writer = std::thread([&os, &out] {
      for (const auto& b : out) b.write_to(os);
    });
  }
  writer.join();
}

struct TextFormat {
  bool fixed = false;        // "%.6f"; else "%.6g", as an ostream and Eigen::IOFormat(6) write it
  const char* na = nullptr;  // written for a non-finite value; nullptr: nan, inf or -inf
  char sep = '\t';
  double scale = 1.0;      // every value is multiplied by this first
  bool by_column = false;  // line i is column i: for a symmetric matrix its row i, read contiguously
  // as Eigen's default `os << M` (with sep ' '): each value right-aligned to
  // the width of the widest
  bool align = false;
};

// One line per row of M (per column with by_column), the values separated by
// fmt.sep. With `rows`, line r is row rows[r] of M, or all fmt.na (nan) when
// rows[r] is negative. A vector is written one value per line. Every line
// ends with '\n', as `os << M.format(IOFormat(6, DontAlignCols, "\t", "\n"))
// << '\n'` wrote it.
void write_matrix(std::ostream& os,
                  const Eigen::Ref<const Mat2D>& M,
                  const TextFormat& fmt = {},
                  const std::vector<Eigen::Index>* rows = nullptr);
// the same into a new file, with an optional header line, or an error
void write_matrix(const std::string& path,
                  const Eigen::Ref<const Mat2D>& M,
                  const TextFormat& fmt = {},
                  const std::string& header = "",
                  const std::vector<Eigen::Index>* rows = nullptr);

// ---------------------------------------------------------------------------
// sample IDs
// ---------------------------------------------------------------------------

// FID and IID of each sample of the input, in its order: from the .fam; the
// .psam, whose #-header names its columns (without one it is laid out like a
// .fam); or the header line of a BEAGLE file.
struct SampleIds {
  String1D fid, iid;    // fid is "0" without an FID column, and the IID for BEAGLE
  bool has_fid = true;  // false for a .psam without an FID column, and for BEAGLE
  std::string error;    // why they could not be read; empty when they were
  size_t size() const { return iid.size(); }
};

SampleIds read_sample_ids(const Param& params);

// .eigvecs2, the PCs as plink2's .eigenvec: #FID IID PC1 .. PCk, then the FID,
// IID and PCs of each sample, tab-separated
void write_eigvecs2(const std::string& path, const Mat2D& U, const SampleIds& ids);

// ---------------------------------------------------------------------------
// the files of a reference PCA (-P/--USV)
// ---------------------------------------------------------------------------

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

// header is  #nsamples,nsnps[,key=value]...  then one singular value per line
void read_sigvals(const std::string& path, uint& N, uint& M, Mat1D& S, UsvTransform* transform = nullptr);

Mat1D read_eigvals(const std::string& path);

// the first k columns of the n rows of a .eigvecs or .loadings
Mat2D read_eigvecs(const std::string& path, int n, int k);

// a tab-separated matrix whose shape is not known
Mat2D read_usv(const std::string& path);

// the allele frequencies, column 7 of a .mbim
Mat1D read_frq(const std::string& path);

// the number of leading reference PCs a two-stage analysis (-P/--USV) uses: -k
// when it was given, else all `navail` PCs in `path`. refuses a -k above navail.
int ref_pcs(const Param& params, int navail, const std::string& path);

#endif  // PCAONE_TEXTIO_
