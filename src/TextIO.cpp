/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/TextIO.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "TextIO.hpp"

#include <charconv>
#include <cstdio>
#include <cstdlib>  // strtod
#include <fstream>
#include <sstream>

#include "Cmd.hpp"
#include "Utils.hpp"

// ---------------------------------------------------------------------------
// numbers to text
// ---------------------------------------------------------------------------

namespace {

constexpr double POW10[11] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10};  // all exact

// "00" .. "99"
struct Digits2 {
  char d[200];
  constexpr Digits2()
      : d() {
    for (int i = 0; i < 100; ++i) d[2 * i] = char('0' + i / 10), d[2 * i + 1] = char('0' + i % 10);
  }
};
constexpr Digits2 DIG2;

// ax * 10^d rounded to an integer as printf rounds: to nearest, ties to even,
// from the exact product. 10^d is exact for d <= 22, and fma gives the error of
// the rounded product, so ax * 10^d = y + e exactly. With y < 2^52, y - n and
// y - n - 0.5 are exact multiples of ulp(y), which is more than 2|e|, so s has
// the sign of the exact y + e - n - 0.5.
inline uint64 round_scaled(double ax, int d) {
  const double p = POW10[d];
  const double y = ax * p;
  const double e = std::fma(ax, p, -y);
  const double f = std::floor(y);
  uint64 n = (uint64)f;
  const double s = ((y - f) - 0.5) + e;
  if (s > 0 || (s == 0 && (n & 1))) ++n;
  return n;
}

// the 6 digits of 100000 <= n < 1000000
inline void six_digits(uint64 n, char* p) {
  const uint32_t a = (uint32_t)(n / 10000), b = (uint32_t)(n / 100 % 100), c = (uint32_t)(n % 100);
  std::memcpy(p, DIG2.d + 2 * a, 2);
  std::memcpy(p + 2, DIG2.d + 2 * b, 2);
  std::memcpy(p + 4, DIG2.d + 2 * c, 2);
}

// nan, inf, -inf
inline char* put_nonfinite(double x, char* p) {
  if (std::isnan(x)) {
    std::memcpy(p, "nan", 3);
    return p + 3;
  }
  if (x < 0) *p++ = '-';
  std::memcpy(p, "inf", 3);
  return p + 3;
}

// 10^(X+1) for X = -5 .. 5, near enough: a guess that the rounding corrects
constexpr double DECADE[11] = {1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6};

inline char* put_snprintf(const char* fmt, double x, char* p) {
  char tmp[64];
  const int n = std::snprintf(tmp, sizeof tmp, fmt, x);
  std::memcpy(p, tmp, (size_t)n);
  return p + n;
}

}  // namespace

char* put_uint(uint64 v, char* p) {
  char tmp[20];
  int t = 20;
  while (v >= 100) {
    std::memcpy(tmp + (t -= 2), DIG2.d + 2 * (v % 100), 2);
    v /= 100;
  }
  if (v >= 10)
    std::memcpy(tmp + (t -= 2), DIG2.d + 2 * v, 2);
  else
    tmp[--t] = char('0' + v);
  std::memcpy(p, tmp + t, 20 - t);
  return p + (20 - t);
}

char* put_f6(double x, char* p) {
  if (!std::isfinite(x)) return put_nonfinite(x, p);
  const double ax = std::fabs(x);
  // round_scaled needs ax * 1e6 < 2^52
  if (!(ax < 1e9)) return put_snprintf(ax < 1e15 ? "%.6f" : "%.6g", x, p);
  if (std::signbit(x)) *p++ = '-';  // "-0.000000" for -0 and what rounds to it, as printf
  const uint64 n = round_scaled(ax, 6);
  p = put_uint(n / 1000000, p);
  *p++ = '.';
  const uint64 f = n % 1000000;
  const uint32_t a = (uint32_t)(f / 10000), b = (uint32_t)(f / 100 % 100), c = (uint32_t)(f % 100);
  std::memcpy(p, DIG2.d + 2 * a, 2);
  std::memcpy(p + 2, DIG2.d + 2 * b, 2);
  std::memcpy(p + 4, DIG2.d + 2 * c, 2);
  return p + 6;
}

// %.6g: let X be the decimal exponent of x rounded to 6 significant digits.
// For -4 <= X < 6 it is x rounded to 5 - X decimals ("%f"), else 6 digits with
// an exponent ("%e"); then the trailing zeros, and a trailing point, go. The
// "%f" range is done here, which holds every PC, loading and frequency; the
// rest goes to snprintf.
char* put_g6(double x, char* p) {
  if (!std::isfinite(x)) return put_nonfinite(x, p);
  if (x == 0) {
    if (std::signbit(x)) *p++ = '-';
    *p++ = '0';
    return p;
  }
  const double ax = std::fabs(x);
  if (!(ax >= 1e-5 && ax < 1e6)) return put_snprintf("%.6g", x, p);
  // X is floor(log10(ax)), or one more when the rounding carries into the
  // next decade (9.999995 -> 10). floor(log2(ax)) * log10(2) is floor(log10(ax))
  // or one less, so X is found going up, at most twice: it is the first
  // exponent at which ax rounds to fewer than 7 digits. Up to 1e7 is rounded,
  // which round_scaled does exactly.
  uint64 bits;
  std::memcpy(&bits, &ax, sizeof bits);
  const int e2 = (int)((bits >> 52) & 0x7ff) - 1023;  // floor(log2(ax)), ax being normal here
  int X = std::max(-5, (e2 * 78913) >> 18);           // floor(e2 * log10(2)), also for e2 < 0
  X += ax >= DECADE[X + 5];                           // nearly always floor(log10(ax)) now, so one rounding
  uint64 n;
  while ((n = round_scaled(ax, 5 - X)) >= 1000000)
    if (++X > 5) return put_snprintf("%.6g", x, p);  // 999999.5 -> "1e+06"
  if (X < -4) return put_snprintf("%.6g", x, p);     // "9.99999e-05"
  *p = '-';
  p += std::signbit(x);
  // fixed-size copies, which compile to a move or two: up to 8 bytes past the
  // end that is returned are written, within TEXT_NUM_MAX
  char dig[16];
  six_digits(n, dig);
  int nd = 6;
  while (dig[nd - 1] == '0') --nd;  // the digits without the trailing zeros; dig[0] is not 0
  if (X >= 0) {                     // X + 1 digits before the point
    const int ni = X + 1;
    std::memcpy(p, dig, 8);
    p += ni;
    if (nd > ni) {
      *p++ = '.';
      std::memcpy(p, dig + ni, 8);
      p += nd - ni;
    }
    return p;
  }
  // 0.000ddd: -X - 1 zeros after the point
  std::memcpy(p, "0.000000", 8);
  p += 1 - X;
  std::memcpy(p, dig, 8);
  return p + nd;
}

// ---------------------------------------------------------------------------
// matrices to text
// ---------------------------------------------------------------------------

void write_matrix(std::ostream& os,
                  const Eigen::Ref<const Mat2D>& M,
                  const TextFormat& fmt,
                  const std::vector<Eigen::Index>* rows) {
  const Eigen::Index nl = rows ? (Eigen::Index)rows->size() : (fmt.by_column ? M.cols() : M.rows());
  const Eigen::Index nv = fmt.by_column ? M.rows() : M.cols();
  if (nl == 0 || nv == 0) return;
  const size_t nalen = fmt.na ? std::strlen(fmt.na) : 3;
  if (nalen > TEXT_NUM_MAX) cao.error("BUG: write_matrix: the NA string is too long");
  const bool fixed = fmt.fixed, by_column = fmt.by_column;
  const char sep = fmt.sep;
  const double scale = fmt.scale;

  // one value; non-finite ones as fmt.na
  auto put = [&](double v, char* p) -> char* {
    v *= scale;
    if (fmt.na && !std::isfinite(v)) {
      std::memcpy(p, fmt.na, nalen);
      return p + nalen;
    }
    return fixed ? put_f6(v, p) : put_g6(v, p);
  };

  // Eigen's alignment: every value padded on the left to the widest
  int width = 0;
  if (fmt.align) {
#pragma omp parallel for reduction(max : width) schedule(static)
    for (Eigen::Index j = 0; j < M.cols(); ++j) {
      char tmp[TEXT_NUM_MAX];
      for (Eigen::Index i = 0; i < M.rows(); ++i) width = std::max(width, (int)(put(M(i, j), tmp) - tmp));
    }
  }

  const size_t maxline = (size_t)nv * (std::max<size_t>(TEXT_NUM_MAX, width) + 1) + 1;
  write_lines(os, nl, (size_t)nv * 10 + 1, [&](Eigen::Index l, TextBuf& buf) {
    char* p = buf.reserve(maxline);
    const Eigen::Index r = rows ? (*rows)[l] : l;
    for (Eigen::Index j = 0; j < nv; ++j) {
      if (j) *p++ = sep;
      const double v = r < 0 ? std::numeric_limits<double>::quiet_NaN() : by_column ? M(j, r) : M(r, j);
      if (width) {
        char tmp[TEXT_NUM_MAX];
        const int len = (int)(put(v, tmp) - tmp);
        for (int s = len; s < width; ++s) *p++ = ' ';
        std::memcpy(p, tmp, len);
        p += len;
      } else {
        p = put(v, p);
      }
    }
    *p++ = '\n';
    buf.commit(p);
  });
}

void write_matrix(const std::string& path,
                  const Eigen::Ref<const Mat2D>& M,
                  const TextFormat& fmt,
                  const std::string& header,
                  const std::vector<Eigen::Index>* rows) {
  std::ofstream ofs(path, std::ios::binary);
  if (!ofs.is_open()) cao.error("can not open file for writing: " + path);
  ofs << header;
  write_matrix(ofs, M, fmt, rows);
  ofs.close();
  if (!ofs) cao.error("error writing " + path);
}

// ---------------------------------------------------------------------------
// text to numbers
// ---------------------------------------------------------------------------

namespace {

constexpr double EXACT_POW10[23] = {1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,  1e10, 1e11,
                                    1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};

inline bool is_blank(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f'; }
inline bool is_digit(char c) { return (unsigned)(c - '0') < 10u; }

// strtod of the number at p (no blanks). strtod needs a terminated string: a
// number holds no blank, so the token up to the next blank is enough
const char* parse_double_strtod(const char* p, const char* end, double& v) {
  const char* q = p;
  while (q < end && !is_blank(*q)) ++q;
  const std::string tok(p, q);
  char* stop = nullptr;
  v = std::strtod(tok.c_str(), &stop);
  return stop == tok.c_str() ? nullptr : p + (stop - tok.c_str());
}

}  // namespace

#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L

// libstdc++ has std::from_chars for double (GCC 11 on), as fast as the path
// below and exact. It reads no '+', and leaves v alone out of range, where
// strtod gives +-inf or a denormal: those few go to strtod.
const char* parse_double(const char* p, const char* end, double& v) {
  while (p < end && is_blank(*p)) ++p;
  const char* s = (p < end && *p == '+' && end - p > 1 && p[1] != '-') ? p + 1 : p;
  const auto r = std::from_chars(s, end, v);
  if (r.ec == std::errc()) return r.ptr;
  if (r.ec == std::errc::result_out_of_range) return parse_double_strtod(p, end, v);
  return nullptr;
}

#else

// libc++, so Apple Clang, has no std::from_chars for double: up to 19
// significant digits with a decimal exponent within +-22 are converted here,
// the rest by strtod
const char* parse_double(const char* p, const char* end, double& v) {
  while (p < end && is_blank(*p)) ++p;
  const char* const start = p;
  bool neg = false;
  if (p < end && (*p == '-' || *p == '+')) neg = *p++ == '-';
  uint64 m = 0;       // the significant digits
  int nd = 0, e = 0;  // how many, and the decimal exponent of the last
  bool any = false;
  for (; p < end && is_digit(*p); ++p) {
    any = true;
    if (m == 0 && *p == '0') continue;  // leading zeros
    if (++nd > 19) return parse_double_strtod(start, end, v);
    m = 10 * m + (uint64)(*p - '0');
  }
  if (p < end && *p == '.') {
    for (++p; p < end && is_digit(*p); ++p) {
      any = true;
      --e;
      if (m == 0 && *p == '0') continue;
      if (++nd > 19) return parse_double_strtod(start, end, v);
      m = 10 * m + (uint64)(*p - '0');
    }
  }
  if (!any) return parse_double_strtod(start, end, v);  // nan, inf, or no number at all
  if (p < end && (*p == 'e' || *p == 'E')) {
    const char* q = p + 1;
    bool eneg = false;
    if (q < end && (*q == '-' || *q == '+')) eneg = *q++ == '-';
    if (q < end && is_digit(*q)) {  // else "1e" is 1, and the 'e' is not read, as strtod does
      int x = 0;
      for (; q < end && is_digit(*q); ++q)
        if (x < 100000) x = 10 * x + (*q - '0');
      e += eneg ? -x : x;
      p = q;
    }
  }
  if (m == 0) {
    v = neg ? -0.0 : 0.0;
    return p;
  }
  if (m > (uint64(1) << 53) || e < -22 || e > 22) return parse_double_strtod(start, end, v);
  // m and 10^|e| are exact doubles, so one operation rounds correctly (Clinger 1990)
  const double d = e < 0 ? (double)m / EXACT_POW10[-e] : (double)m * EXACT_POW10[e];
  v = neg ? -d : d;
  return p;
}

#endif

LineReader::LineReader(const std::string& path, size_t block_bytes)
    : path_(path),
      block_(std::max<size_t>(block_bytes, 1 << 16)) {
  fp_ = std::fopen(path.c_str(), "rb");
  if (!fp_) cao.error("can not open ", path);
  for (auto& b : buf_) b.resize(block_ + 1);
  fill(0, 0);
}

LineReader::~LineReader() {
  if (ahead_.joinable()) ahead_.join();
  if (fp_) std::fclose(fp_);
}

void LineReader::fill(int which, size_t offset) {
  std::vector<char>& b = buf_[which];
  size_t len = offset;
  while (len < b.size() - 1) {
    const size_t got = std::fread(b.data() + len, 1, b.size() - 1 - len, fp_);
    len += got;
    if (got == 0) break;
  }
  if (std::ferror(fp_)) failed_ = true;
  len_[which] = len;
  eof_[which] = len < b.size() - 1;
}

bool LineReader::next_block(std::vector<LineSpan>& lines) {
  lines.clear();
  if (done_) return false;
  if (ahead_.joinable()) {
    ahead_.join();
    cur_ = 1 - cur_;
  }
  if (failed_) cao.error("error reading ", path_);
  std::vector<char>& b = buf_[cur_];
  size_t len = len_[cur_];
  if (len == 0) return false;
  // the whole lines in b; a line longer than the block grows it
  size_t cut = len;
  if (!eof_[cur_]) {
    for (;;) {
      size_t i = len;
      while (i > 0 && b[i - 1] != '\n') --i;
      if (i > 0) {
        cut = i;
        break;
      }
      b.resize(2 * b.size());
      fill(cur_, len);
      if (failed_) cao.error("error reading ", path_);
      len = len_[cur_];
      if (eof_[cur_]) {
        cut = len;
        break;
      }
    }
  }
  if (eof_[cur_]) {
    done_ = true;
  } else {  // the start of the next line goes first in the other buffer, read behind it on another thread
    const int other = 1 - cur_;
    const size_t tail = len - cut;
    if (buf_[other].size() < std::max(block_, 2 * tail) + 1) buf_[other].resize(std::max(block_, 2 * tail) + 1);
    std::memcpy(buf_[other].data(), b.data() + cut, tail);
    ahead_ = std::thread([this, other, tail] { fill(other, tail); });
  }
  const char* p = b.data();
  const char* const stop = b.data() + cut;
  while (p < stop) {
    const char* nl = static_cast<const char*>(std::memchr(p, '\n', stop - p));
    const char* e = nl ? nl : stop;
    const char* q = p;
    while (q < e && is_blank(*q)) ++q;
    if (q < e) lines.push_back({p, e});
    p = e + 1;
  }
  return true;
}

// ---------------------------------------------------------------------------
// sample IDs
// ---------------------------------------------------------------------------

// the IDs of a BEAGLE file: its header is marker allele1 allele2, then the
// sample name three times for each sample
static SampleIds read_beagle_ids(const std::string& path) {
  SampleIds ids;
  ids.has_fid = false;
  gzFile fp = gzopen(path.c_str(), "r");
  if (!fp) {
    ids.error = "can not open " + path;
    return ids;
  }
  std::string header;
  char chunk[1 << 16];
  while (gzgets(fp, chunk, sizeof chunk)) {
    header += chunk;
    if (!header.empty() && header.back() == '\n') break;
  }
  gzclose(fp);
  std::istringstream iss(header);
  std::string tok;
  int ncol = 0;
  while (iss >> tok)
    if (++ncol > 3 && (ncol - 1) % 3 == 0) ids.iid.push_back(tok);
  if (ncol % 3) {
    ids.iid.clear();
    ids.error = "the header of " + path + " has " + std::to_string(ncol) + " columns, not a multiple of 3";
  }
  ids.fid = ids.iid;
  return ids;
}

SampleIds read_sample_ids(const Param& params) {
  if (params.file_t == FileType::BEAGLE) return read_beagle_ids(params.filein);
  SampleIds ids;
  const std::string path = params.filein + (params.file_t == FileType::PGEN ? ".psam" : ".fam");
  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    ids.error = "can not open " + path;
    return ids;
  }
  int cfid = 0, ciid = 1;  // FID IID ..., as a .fam
  std::string line, tok;
  String1D t;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    const bool header = line[0] == '#';
    if (header && line.size() > 1 && line[1] == '#') continue;  // ## lines of a .psam
    std::istringstream iss(header ? line.substr(1) : line);
    t.clear();
    while (iss >> tok) t.push_back(tok);
    if (header) {  // #FID IID ... or #IID ...
      cfid = ciid = -1;
      for (int c = 0; c < (int)t.size(); ++c) {
        if (t[c] == "FID") cfid = c;
        if (t[c] == "IID") ciid = c;
      }
      ids.has_fid = cfid >= 0;
      if (ciid < 0) {
        ids.error = "the header line of " + path + " has no IID column";
        break;
      }
      continue;
    }
    if ((int)t.size() <= std::max(cfid, ciid)) {
      ids.error = "too few columns in " + path + " at: " + line;
      break;
    }
    ids.fid.push_back(cfid >= 0 ? t[cfid] : std::string("0"));
    ids.iid.push_back(t[ciid]);
  }
  if (!ids.error.empty()) ids.fid.clear(), ids.iid.clear();
  return ids;
}

void write_eigvecs2(const std::string& path, const Mat2D& U, const SampleIds& ids, const char* na) {
  if (!ids.error.empty() || (Eigen::Index)ids.size() != U.rows()) {
    cao.warn("no .eigvecs2 written:",
             ids.error.empty() ? std::to_string(ids.size()) + " sample IDs for " + std::to_string(U.rows()) + " samples"
                               : ids.error);
    return;
  }
  std::ofstream ofs(path, std::ios::binary);
  if (!ofs.is_open()) cao.error("can not open file for writing: " + path);
  ofs << "#FID\tIID";
  for (Eigen::Index k = 0; k < U.cols(); ++k) ofs << "\tPC" << k + 1;
  ofs << '\n';
  const Eigen::Index K = U.cols();
  const size_t nalen = na ? std::strlen(na) : 0;
  write_lines(ofs, U.rows(), 32 + 10 * K, [&](Eigen::Index i, TextBuf& buf) {
    buf.append(ids.fid[i]);
    buf.append("\t", 1);
    buf.append(ids.iid[i]);
    char* p = buf.reserve(K * (std::max(TEXT_NUM_MAX, nalen) + 1) + 1);
    for (Eigen::Index k = 0; k < K; ++k) {
      *p++ = '\t';
      if (na && !std::isfinite(U(i, k))) {
        std::memcpy(p, na, nalen);
        p += nalen;
      } else {
        p = put_g6(U(i, k), p);
      }
    }
    *p++ = '\n';
    buf.commit(p);
  });
  ofs.close();
  if (!ofs) cao.error("error writing " + path);
}

// ---------------------------------------------------------------------------
// the files of a reference PCA
// ---------------------------------------------------------------------------

// The key=value fields were added later; std::stoi stops at the next comma, so
// an older PCAone reads a newer file correctly and simply ignores them, and a
// newer PCAone reading an older file just finds none.
void read_sigvals(const std::string& path, uint& N, uint& M, Mat1D& S, UsvTransform* transform) {
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("can not open the singular values file\n => " + path);
  std::string line;
  getline(fin, line);
  // parse line #nsamples,nsnps
  size_t comma = line.find(',');
  N = std::stoi(line.substr(1, comma - 1));
  M = std::stoi(line.substr(comma + 1));
  if (transform != nullptr) {
    for (const auto& tok : split_string(line, ",")) {
      const size_t eq = tok.find('=');
      if (eq == std::string::npos) continue;
      const std::string key = tok.substr(0, eq), val_s = tok.substr(eq + 1);
      if (key == "scale")
        transform->scale = std::stoi(val_s);
      else if (key == "ploidy")
        transform->ploidy = std::stoi(val_s);
      else if (key == "gscale")
        transform->gscale = std::stoi(val_s);
      else
        continue;
      transform->known = true;
    }
  }
  Double1D V;
  while (getline(fin, line)) V.push_back(std::stod(line));
  S = Eigen::Map<Mat1D>(V.data(), V.size());
}

Mat1D read_eigvals(const std::string& path) {
  Double1D V;
  std::ifstream fin(path);
  std::string line;
  while (getline(fin, line)) V.push_back(std::stod(line));
  return Eigen::Map<Mat1D>(V.data(), V.size());
}

Mat2D read_eigvecs(const std::string& path, int n, int k) {
  cao.print(tick.date(), "read matrix from " + path);
  if (n <= 0 || k <= 0) cao.error("read_eigvecs: asked for a ", n, " x ", k, " matrix from ", path);
  // Every bound is checked before the write it guards: a .eigvecs from a
  // bigger cohort, or a short row, is an error, never a write past the end.
  // Blank lines are skipped, and the values after the k-th of a row ignored.
  Mat2D M(n, k);
  const ReadLines r = read_lines(
      path,
      [&](Eigen::Index j, const char* b, const char* e) -> int {
        if (j >= n) return -1;
        for (int c = 0; c < k; ++c) {
          double v;
          if (!(b = parse_double(b, e, v))) return c + 1;
          M(j, c) = v;
        }
        return 0;
      },
      [](Eigen::Index) {});
  if (r.code == -1)
    cao.error(path, " has more than the ", n, " rows this dataset needs. is it from a different cohort?");
  if (r.code > 0)
    cao.error(path, ": row ", r.bad + 1, " has only ", r.code - 1, " values but ", k,
              " are needed. rerun the reference with a larger -k, or use a smaller one here");
  if (r.lines != n) cao.error(path, " has ", r.lines, " rows but this dataset needs ", n);
  return M;
}

Mat2D read_usv(const std::string& path) {
  std::vector<Double1D> rows;
  // the tab-separated values of each row, every one a number
  const ReadLines r = read_lines(
      path,
      [&](Eigen::Index j, const char* b, const char* e) -> int {
        Double1D& row = rows[j];
        for (;;) {
          const char* t = static_cast<const char*>(std::memchr(b, '\t', e - b));
          const char* fe = t ? t : e;
          double v;
          const char* q = parse_double(b, fe, v);
          if (!q) return 1;
          while (q < fe && is_blank(*q)) ++q;
          if (q != fe) return 1;
          row.push_back(v);
          if (!t) return 0;
          b = t + 1;
        }
      },
      [&](Eigen::Index n) { rows.resize(n); });
  if (r.code) cao.error(path, ": row ", r.bad + 1, " has a value that is not a number");
  const Eigen::Index nrow = (Eigen::Index)rows.size(), ncol = nrow ? (Eigen::Index)rows[0].size() : 0;
  if (nrow < 1 || ncol < 1) cao.error("can not parse a matrix from file\n =>" + path);
  Mat2D M(nrow, ncol);
  for (Eigen::Index i = 0; i < nrow; ++i) {
    if ((Eigen::Index)rows[i].size() != ncol) cao.error("the columns are not aligned!\n =>" + path);
    for (Eigen::Index j = 0; j < ncol; ++j) M(i, j) = rows[i][j];
  }
  return M;
}

Mat1D read_frq(const std::string& path) {
  // without this the file simply reads as empty, F stays size 0, and the first
  // F(snp_idx) downstream segfaults with no message. -P/--USV points filebim at
  // <prefix>.mbim, which a PCA run does not write with --no-loadings (or, before
  // v0.8.0, without -V/--printv) or from input with no .bim/.pvar.
  if (!std::ifstream(path).is_open()) cao.error("can not open the allele frequency file\n => " + path);
  Mat1D F;
  // the 7th of exactly 7 fields, each blank (' ' or '\t') ending one, as split_string splits them
  const ReadLines r = read_lines(
      path,
      [&](Eigen::Index j, const char* b, const char* e) -> int {
        int nsep = 0;
        const char* f7 = nullptr;
        for (const char* q = b; q < e; ++q)
          if ((*q == ' ' || *q == '\t') && ++nsep == 6) f7 = q + 1;
        if (nsep != 6) return 1;
        double v;
        if (!parse_double(f7, e, v)) return 2;
        F(j) = v;
        return 0;
      },
      [&](Eigen::Index n) { F.conservativeResize(n); });
  if (r.code == 1) cao.error("the input file is not valid: line ", r.bad + 1, " does not have 7 columns\n => " + path);
  if (r.code == 2) cao.error("the input file is not valid: no allele frequency on line ", r.bad + 1, "\n => " + path);
  if (F.size() == 0) cao.error("no allele frequencies found in\n => " + path);
  return F;
}

int ref_pcs(const Param& params, int navail, const std::string& path) {
  if (navail <= 0) cao.error(path, "has no PCs");
  if (params.ref_k == 0) return navail;
  if ((int)params.ref_k > navail) cao.error("-k/--pc", params.ref_k, "is larger than the", navail, "PCs in", path);
  return params.ref_k;
}
