// The text I/O of TextIO.cpp against what it replaces.
//
// Output: put_g6, put_f6 and write_matrix format every value as printf("%.6g")
// / printf("%.6f") does, and a matrix comes out byte for byte as Eigen's
// IOFormat(6, DontAlignCols, "\t", "\n") and default `os << M` wrote it,
// whatever the number of threads. The values are random doubles over every
// magnitude, plus those where rounding is hardest: ties, the edges of each
// decade and of the "%f" range of "%g".
//
// Input: parse_double gives the bits strtod gives, and ends where it ends, on
// what PCAone writes and on what takes its slow path; LineReader returns every
// line in order across blocks, of any length; read_eigvecs, read_usv and
// read_frq read what getline + strtod read, and stop with an error where they
// did.
#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#define _DECLARE_TOOLBOX_HERE  // the logger and timer of the library, as Main.cpp
#include "../src/TextIO.hpp"
#include "../src/Utils.hpp"

static int failures = 0;
#define CHECK(cond, ...)                                                            \
  do {                                                                              \
    if (!(cond)) {                                                                  \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed: ", __FILE__, __LINE__, #cond); \
      std::fprintf(stderr, __VA_ARGS__);                                            \
      std::fprintf(stderr, "\n");                                                   \
      if (++failures > 20) std::exit(1);                                            \
    }                                                                               \
  } while (0)

static std::string printf_of(const char* fmt, double x) {
  char b[512];
  std::snprintf(b, sizeof b, fmt, x);
  std::string s(b);
  if (s == "-nan") s = "nan";  // put_* drop the sign of a NaN
  return s;
}

static void check_one(double x) {
  char b[64];
  const std::string g(b, put_g6(x, b) - b), want_g = printf_of("%.6g", x);
  CHECK(g == want_g, "put_g6(%.17g) = '%s', printf '%s'", x, g.c_str(), want_g.c_str());
  CHECK(g.size() <= TEXT_NUM_MAX, "put_g6(%.17g) is %zu bytes", x, g.size());
  if (std::fabs(x) < 1e15 || !std::isfinite(x)) {
    const std::string f(b, put_f6(x, b) - b), want_f = printf_of("%.6f", x);
    CHECK(f == want_f, "put_f6(%.17g) = '%s', printf '%s'", x, f.c_str(), want_f.c_str());
    CHECK(f.size() <= TEXT_NUM_MAX, "put_f6(%.17g) is %zu bytes", x, f.size());
  }
}

// parse_double against strtod: the same bits (any NaN for a NaN), and the
// same end
static void check_parse(const std::string& s) {
  char* stop = nullptr;
  const double want = std::strtod(s.c_str(), &stop);
  const long want_len = stop - s.c_str();
  double got = 12345.0;
  const char* q = parse_double(s.data(), s.data() + s.size(), got);
  if (want_len == 0) {
    CHECK(q == nullptr, "parse_double('%s') read a number, strtod none", s.c_str());
    return;
  }
  CHECK(q != nullptr, "parse_double('%s') read no number, strtod %.17g", s.c_str(), want);
  if (!q) return;
  const bool same = std::isnan(want) ? std::isnan(got) : std::memcmp(&want, &got, sizeof got) == 0;
  CHECK(same, "parse_double('%s') = %.17g, strtod %.17g", s.c_str(), got, want);
  CHECK(q - s.data() == want_len, "parse_double('%s') read %ld chars, strtod %ld", s.c_str(), (long)(q - s.data()),
        want_len);
}

static std::string tmp_path(const std::string& name) {
  return (std::filesystem::temp_directory_path() / ("pcaone_test_textio_" + name)).string();
}

static void write_file(const std::string& path, const std::string& text) {
  std::ofstream(path, std::ios::binary) << text;
}

// the old read_eigvecs: getline, then strtod for the first k values of each line that is not blank
static Mat2D read_ref(const std::string& path, int k) {
  std::ifstream in(path);
  std::string line;
  std::vector<std::vector<double>> rows;
  while (std::getline(in, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
    if (line.find_first_not_of(" \t") == std::string::npos) continue;
    const char* p = line.c_str();
    rows.emplace_back();
    for (int c = 0; c < k; ++c) {
      char* e = nullptr;
      rows.back().push_back(std::strtod(p, &e));
      p = e;
    }
  }
  Mat2D M(rows.size(), k);
  for (size_t i = 0; i < rows.size(); ++i)
    for (int c = 0; c < k; ++c) M(i, c) = rows[i][c];
  return M;
}

// the message cao.error() threw, its runs of spaces made one, or "" when it did not
template <typename F>
static std::string error_of(F&& f) {
  try {
    f();
  } catch (const std::runtime_error& e) {
    std::string m;
    for (const char* c = e.what(); *c; ++c)
      if (!(*c == ' ' && !m.empty() && m.back() == ' ')) m += *c;
    return m;
  }
  return "";
}

static bool bitwise_equal(const Mat2D& A, const Mat2D& B) {
  return A.rows() == B.rows() && A.cols() == B.cols() &&
         std::memcmp(A.data(), B.data(), sizeof(double) * A.size()) == 0;
}

static void test_input(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  long n = 0;
  char b[64];
  // what PCAone writes, and the same values at full precision
  for (int i = 0; i < 2000000; ++i, n += 3) {
    const double x = std::pow(10.0, -12.0 + 24.0 * u(rng)) * (u(rng) < 0.5 ? -1 : 1);
    check_parse(std::string(b, put_g6(x, b)));
    std::snprintf(b, sizeof b, "%.17g", x);
    check_parse(b);
    std::snprintf(b, sizeof b, i % 2 ? "%.6f" : "%.9e", x);
    check_parse(b);
  }
  // long mantissas, far exponents, 2^53 +- 1: the slow path
  for (int i = 0; i < 200000; ++i, ++n) {
    std::string s = u(rng) < 0.5 ? "-" : "";
    const int nd = 1 + (int)(u(rng) * 30);
    for (int d = 0; d < nd; ++d) s += char('0' + (int)(u(rng) * 10));
    if (u(rng) < 0.5) s.insert(s.size() - (size_t)(u(rng) * nd), ".");
    if (u(rng) < 0.7) s += "e" + std::to_string((int)(u(rng) * 700) - 350);
    check_parse(s);
  }
  const char* special[] = {"0", "-0", "+0", "0.0", "-0.000", "1", "-1", "+1.5", "1.", ".5", "-.5", ".", "-", "+", "+-1",
                           "-+1", "1e", "1e+", "1e-", "2e3", "2E-3", "1.5abc", "  \t 3.25", "1,5", "nan", "-nan",
                           "NaN", "nan(123)", "inf", "-inf", "+inf", "Infinity", "infin", "NA", "-NA", "", " ",
                           "9007199254740992", "9007199254740993", "9007199254740993e-5", "1e22", "1e23", "123e-22",
                           "123e-23", "1e400", "-1e400", "1e-400", "5e-324", "2.2250738585072011e-308",
                           "4.9406564584124654e-324", "1.7976931348623157e308", "1.7976931348623159e308",
                           "0.000000000000000000000000000001", "100000000000000000000000000000", "1e99999999999"};
  for (const char* s : special) check_parse(s), ++n;
  std::printf("parse_double: %ld strings, the bits and ends of strtod\n", n);

  // LineReader: every line in order, across blocks and a line longer than one
  {
    const std::string path = tmp_path("lines");
    std::vector<std::string> want;
    std::string text;
    for (int i = 0; i < 40000; ++i) {
      std::string line;
      const int len = i == 777 ? 300000 : (int)(u(rng) * 200);
      for (int c = 0; c < len; ++c) line += char('a' + (int)(u(rng) * 26));
      if (i % 97 == 0) line = " \t ";  // blank: skipped
      if (i % 89 == 0) line += '\r';
      text += line + "\n";
      if (line.find_first_not_of(" \t\r") != std::string::npos) want.push_back(line);
    }
    text += "last line without a newline";
    want.push_back("last line without a newline");
    write_file(path, text);
    LineReader rd(path, 1 << 16);
    std::vector<LineSpan> spans;
    size_t i = 0, blocks = 0;
    bool ok = true;
    while (rd.next_block(spans)) {
      ++blocks;
      for (const auto& sp : spans) ok = ok && i < want.size() && std::string(sp.b, sp.e) == want[i++];
    }
    CHECK(ok && i == want.size(), "LineReader: %zu of %zu lines", i, want.size());
    CHECK(blocks > 10, "LineReader: only %zu blocks", blocks);
    std::filesystem::remove(path);
    // an empty file
    write_file(path, "");
    LineReader empty(path, 1 << 16);
    CHECK(!empty.next_block(spans), "LineReader: an empty file has lines");
    std::filesystem::remove(path);
  }

  for (int threads : {1, 2, 8}) {
    omp_set_num_threads(threads);
    // read_eigvecs: what write_matrix wrote, as getline + strtod read it
    const std::string path = tmp_path("matrix");
    Mat2D V(20000, 7);
    for (Eigen::Index j = 0; j < V.cols(); ++j)
      for (Eigen::Index i = 0; i < V.rows(); ++i) V(i, j) = std::pow(10.0, -6 + 7 * u(rng)) * (u(rng) - 0.5);
    write_matrix(path, V);
    CHECK(bitwise_equal(read_eigvecs(path, 20000, 7), read_ref(path, 7)), "read_eigvecs, %d threads", threads);
    CHECK(bitwise_equal(read_eigvecs(path, 20000, 3), read_ref(path, 3)), "read_eigvecs -k 3, %d threads", threads);
    CHECK(bitwise_equal(read_usv(path), read_ref(path, 7)), "read_usv, %d threads", threads);
    // CRLF line ends and blank lines
    write_file(path, "1\t2\t3\r\n\r\n4\t5\t6\r\n \t\n7\t8\t9\r\n");
    Mat2D W(3, 3);
    W << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    CHECK(bitwise_equal(read_eigvecs(path, 3, 3), W), "read_eigvecs of CRLF lines");
    CHECK(bitwise_equal(read_usv(path), W), "read_usv of CRLF lines");
    // the errors of the old reader
    std::string err = error_of([&] { read_eigvecs(path, 2, 3); });
    CHECK(err.find("has more than the 2 rows") != std::string::npos, "too many rows: '%s'", err.c_str());
    err = error_of([&] { read_eigvecs(path, 4, 3); });
    CHECK(err.find("has 3 rows but this dataset needs 4") != std::string::npos, "too few rows: '%s'", err.c_str());
    err = error_of([&] { read_eigvecs(path, 3, 4); });
    CHECK(err.find("row 1 has only 3 values but 4") != std::string::npos, "short row: '%s'", err.c_str());
    write_file(path, "1\t2\n3\tNA\n");
    err = error_of([&] { read_eigvecs(path, 2, 2); });
    CHECK(err.find("row 2 has only 1 values") != std::string::npos, "NA: '%s'", err.c_str());
    err = error_of([&] { read_usv(path); });
    CHECK(err.find("row 2 has a value that is not a number") != std::string::npos, "read_usv NA: '%s'",
          err.c_str());
    write_file(path, "1\t2\n3\n");
    err = error_of([&] { read_usv(path); });
    CHECK(err.find("not aligned") != std::string::npos, "read_usv ragged: '%s'", err.c_str());
    // read_frq: the 7th column of a .mbim
    write_file(path, "1\trs1\t0\t100\tA\tG\t0.25\n1 rs2 0 200 C T 0.125\r\n1\trs3\t0\t300\tA\tC\t1e-3\n");
    const Mat1D F = read_frq(path);
    CHECK(F.size() == 3 && F(0) == 0.25 && F(1) == 0.125 && F(2) == 1e-3, "read_frq");
    write_file(path, "1\trs1\t0\t100\tA\tG\t0.25\n1\trs2\t0\t200\tC\tT\n");
    err = error_of([&] { read_frq(path); });
    CHECK(err.find("line 2 does not have 7 columns") != std::string::npos, "read_frq 6 columns: '%s'", err.c_str());
    err = error_of([&] { read_eigvecs(tmp_path("missing"), 1, 1); });
    CHECK(err.find("can not open") != std::string::npos, "missing file: '%s'", err.c_str());
    std::filesystem::remove(path);
  }
  std::printf("read_eigvecs, read_usv, read_frq: the values of getline + strtod and their errors, at 1, 2, 8 threads\n");
}

int main() {
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  long n = 0;

  // special values
  const double specials[] = {0.0, -0.0, 1.0, -1.0, 0.5, 1e-5, 1e-4, 9.99999e-5, 9.999995e-5, 0.0001, 0.00010000005,
                             999999.5, 999999.4, 999999.6, 100000.5, 100001.5, 1e6, 1e5, 123456.5, 2.5e-7, 5e-7,
                             1.0000005, 0.12345650000000001, 9.9999949999, 9.999995, 9.9999950000001, 1e300, -1e-300,
                             5e-324, std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN(), -std::numeric_limits<double>::quiet_NaN(),
                             1e9, 999999999.9999995, 1e15, 4503599627370496.0};
  for (double x : specials) check_one(x), check_one(-x), n += 2;

  // random magnitudes, 1e-12 .. 1e12, both signs
  for (int i = 0; i < 4000000; ++i, ++n) {
    const double x = std::pow(10.0, -12.0 + 24.0 * u(rng)) * (u(rng) < 0.5 ? -1 : 1);
    check_one(x);
  }
  // loadings, PCs and frequencies: [-1, 1] uniformly, and around 1/sqrt(M)
  for (int i = 0; i < 4000000; ++i, ++n) check_one(2.0 * u(rng) - 1.0), check_one(1e-3 * (u(rng) - 0.5));
  // ties and near-ties: k / 2^j has a short binary expansion, and m * 10^-d
  // +- 0.5 ulp sits on a rounding boundary of "%.6g" or "%.6f"
  for (int i = 0; i < 2000000; ++i, ++n) {
    const int e = (int)(u(rng) * 40) - 20;
    const double m = std::floor(u(rng) * 9e6 + 1e6) + 0.5;  // 7 digits and a half
    double x = std::ldexp(std::floor(u(rng) * (1 << 24)), -(int)(u(rng) * 30));
    check_one(x);
    x = m * std::pow(10.0, e);
    check_one(x), check_one(std::nextafter(x, 0.0)), check_one(std::nextafter(x, 1e300));
  }
  // every decade edge
  for (int e = -12; e <= 12; ++e) {
    const double x = std::pow(10.0, e);
    for (double y : {x, std::nextafter(x, 0.0), std::nextafter(x, 1e300), x * (1 - 5e-7), x * (1 - 4.9e-7)})
      check_one(y), check_one(-y), n += 2;
  }
  std::printf("put_g6/put_f6: %ld values checked against printf\n", n);

  // matrices, against Eigen's ostream output
  for (int threads : {1, 2, 3, 8}) {
    omp_set_num_threads(threads);
    for (auto shape : {std::pair<int, int>{1, 1}, {7, 1}, {1, 9}, {1000, 10}, {3001, 4}}) {
      Mat2D M(shape.first, shape.second);
      for (Eigen::Index j = 0; j < M.cols(); ++j)
        for (Eigen::Index i = 0; i < M.rows(); ++i) M(i, j) = std::pow(10.0, -6 + 8 * u(rng)) * (u(rng) - 0.5);
      std::ostringstream want, got;
      want << M.format(Eigen::IOFormat(6, Eigen::DontAlignCols, "\t", "\n")) << '\n';
      write_matrix(got, M);
      CHECK(got.str() == want.str(), "%d x %d matrix, %d threads", shape.first, shape.second, threads);
      // the vector of the eigenvalues
      std::ostringstream wv, gv;
      const Mat1D v = M.col(0);
      wv << v.format(Eigen::IOFormat(6, Eigen::DontAlignCols, "\t", "\n")) << '\n';
      write_matrix(gv, v);
      CHECK(gv.str() == wv.str(), "vector of %d, %d threads", shape.first, threads);
      // the .cov: Eigen's default format, right-aligned to the widest value
      std::ostringstream wa, ga;
      wa << M << '\n';
      TextFormat al;
      al.sep = ' ';
      al.align = true;
      write_matrix(ga, M, al);
      CHECK(ga.str() == wa.str(), "aligned %d x %d matrix, %d threads", shape.first, shape.second, threads);
    }
    // a symmetric matrix by columns, scaled and "%.6f", as EvalAdmix wrote it
    Mat2D S = Mat2D::Random(301, 301);
    S = (S + S.transpose()).eval() / 2;
    S(3, 5) = S(5, 3) = std::numeric_limits<double>::quiet_NaN();
    std::string want;
    char b[64];
    for (Eigen::Index i = 0; i < S.rows(); ++i) {
      for (Eigen::Index j = 0; j < S.cols(); ++j) {
        if (j) want += '\t';
        want += printf_of("%.6f", S(i, j) * 0.5);
      }
      want += '\n';
    }
    std::ostringstream got;
    TextFormat sym;
    sym.fixed = true;
    sym.scale = 0.5;
    sym.by_column = true;
    write_matrix(got, S, sym);
    CHECK(got.str() == want, "symmetric matrix by columns, %d threads", threads);
    // NA for what is not finite, and rows picked by index (-1: all NA)
    Mat2D R(3, 2);
    R << 1, std::numeric_limits<double>::infinity(), 0.25, -2, std::numeric_limits<double>::quiet_NaN(), 3;
    std::vector<Eigen::Index> rows{2, -1, 0, 1};
    std::ostringstream gn;
    TextFormat na;
    na.na = "NA";
    write_matrix(gn, R, na, &rows);
    CHECK(gn.str() == "NA\t3\nNA\tNA\n1\tNA\n0.25\t-2\n", "NA and row order: '%s'", gn.str().c_str());
    (void)b;
  }
  std::printf("write_matrix: same bytes as Eigen's output, at 1, 2, 3 and 8 threads\n");
  test_input(rng);
  if (failures) {
    std::fprintf(stderr, "%d failure(s)\n", failures);
    return 1;
  }
  return 0;
}
