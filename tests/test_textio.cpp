// put_g6, put_f6 and write_matrix (TextIO.cpp) against what they replace:
// every value is formatted as printf("%.6g") / printf("%.6f") formats it, and
// a matrix comes out byte for byte as Eigen's IOFormat(6, DontAlignCols, "\t",
// "\n") and default `os << M` wrote it, whatever the number of threads. The
// values are random doubles over every magnitude, plus those where rounding is
// hardest: ties, the edges of each decade and of the "%f" range of "%g".
#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <sstream>
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
  if (failures) {
    std::fprintf(stderr, "%d failure(s)\n", failures);
    return 1;
  }
  std::printf("write_matrix: same bytes as Eigen's output, at 1, 2, 3 and 8 threads\n");
  return 0;
}
