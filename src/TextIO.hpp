/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/TextIO.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 *
 * Text input and output: the files of a PCA (.eigvecs, .eigvals, .sigvals,
 * .loadings, .mbim), the sample IDs of the input, and the fast, parallel
 * writer and reader that every large text file goes through.
 ******************************************************************************/
#ifndef PCAONE_TEXTIO_
#define PCAONE_TEXTIO_

#include <omp.h>

#include <cstdio>
#include <cstring>
#include <limits>
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
// text to numbers
// ---------------------------------------------------------------------------
//
// Reading the M x k .loadings of a -P run with getline and strtod took ~90 ns
// a value on one thread: strtod goes through the locale, and the lines were
// read and parsed one by one. std::from_chars (C++17) converts without the
// locale, exactly and ~4x faster, but libc++ -- so Apple Clang -- has no
// from_chars for double (deleted up to LLVM 18 at least), and there an exact
// parser of the decimals PCAone writes takes its place. read_lines reads a
// file in blocks, the next block while this one is parsed, and parses the
// lines of a block on every thread.

// The number after any blanks at p, as strtod reads it, but within [p, end) and
// without the locale: its value in v and its end, or nullptr when there is no
// number. std::from_chars where the library has it; else up to 19 significant
// digits with a decimal exponent within +-22 -- all PCAone writes -- are one
// exact multiplication or division of two exact doubles, so correctly rounded
// (Clinger 1990). Out of range, nan, inf and longer numbers go to strtod. So
// every value is the double strtod gives. Hexadecimal is not read.
const char* parse_double(const char* p, const char* end, double& v);

struct LineSpan {
  const char* b;  // a line without its '\n'
  const char* e;
};

// The lines of a file, a block at a time, the next block being read on another
// thread while the caller parses this one. A line may be any length.
class LineReader {
 public:
  LineReader(const std::string& path, size_t block_bytes);
  ~LineReader();
  // the lines of the next block that are not blank (blanks are " \t\r\v\f"),
  // valid until the next call; false at the end of the file
  bool next_block(std::vector<LineSpan>& lines);

 private:
  void fill(int which, size_t offset);  // read into buf_[which] after the first `offset` bytes
  std::string path_;
  std::FILE* fp_ = nullptr;
  size_t block_;
  std::vector<char> buf_[2];
  size_t len_[2] = {0, 0};
  bool eof_[2] = {false, false};
  bool failed_ = false, done_ = false;
  int cur_ = 0;
  std::thread ahead_;
};

struct ReadLines {
  Eigen::Index lines = 0;  // the lines not blank, up to the end of the last block parsed
  Eigen::Index bad = -1;   // the first line whose parser returned a code, if any
  int code = 0;            // and that code
};

// Call line(i, b, e) for the i-th line [b, e) of the file that is not blank, on
// every thread, a block of about TEXT_BLOCK_BYTES per thread at a time, and
// grow(n) before the lines of a block, n being the lines up to its end, so the
// destination can be sized. line returns 0, or a code of the caller's for an
// error; the first line, in file order, with a code ends the reading. line must
// not throw.
template <typename LineFn, typename GrowFn>
ReadLines read_lines(const std::string& path, LineFn&& line, GrowFn&& grow) {
  const int T = omp_in_parallel() ? 1 : std::max(1, omp_get_max_threads());
  LineReader reader(path, T * TEXT_BLOCK_BYTES);
  ReadLines res;
  std::vector<LineSpan> spans;
  while (reader.next_block(spans)) {
    const Eigen::Index first = res.lines, nb = (Eigen::Index)spans.size();
    res.lines += nb;
    grow(res.lines);
    Eigen::Index bad = std::numeric_limits<Eigen::Index>::max();
    int code = 0;
#pragma omp parallel num_threads(T)
    {
      Eigen::Index my_bad = std::numeric_limits<Eigen::Index>::max();
      int my_code = 0;
#pragma omp for schedule(static)
      for (Eigen::Index i = 0; i < nb; ++i) {
        const int c = line(first + i, spans[i].b, spans[i].e);
        if (c && first + i < my_bad) my_bad = first + i, my_code = c;
      }
#pragma omp critical(pcaone_read_lines)
      if (my_bad < bad) bad = my_bad, code = my_code;
    }
    if (code) {
      res.bad = bad;
      res.code = code;
      break;
    }
  }
  return res;
}

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
