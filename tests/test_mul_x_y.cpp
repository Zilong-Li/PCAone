// mul_X_Y (Utils.cpp), the product H = X * G and H += X * G of the RSVD:
// the row-panel split gives the product's values (against a plain loop), for
// assignment and accumulation, for heights that are and are not multiples of
// the panel rounding, and the same bytes whatever the number of threads once
// there are two panels (192 rows). Without an external BLAS the split must
// actually be taken, or this test would only check Eigen's own product.
#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#define _DECLARE_TOOLBOX_HERE  // the logger and timer of the library, as Main.cpp
#include "../src/Utils.hpp"

static int failures = 0;
#define CHECK(cond, ...)                                              \
  do {                                                                \
    if (!(cond)) {                                                    \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed: ", __FILE__, __LINE__, #cond); \
      std::fprintf(stderr, __VA_ARGS__);                              \
      std::fprintf(stderr, "\n");                                     \
      ++failures;                                                     \
    }                                                                 \
  } while (0)

static Mat2D random(Eigen::Index r, Eigen::Index c, std::mt19937_64& rng) {
  std::normal_distribution<double> d;
  Mat2D m(r, c);
  for (Eigen::Index j = 0; j < c; ++j)
    for (Eigen::Index i = 0; i < r; ++i) m(i, j) = d(rng);
  return m;
}

// out (+)= X * Y by dot products in long double, the reference
static Mat2D reference(const Mat2D& X, const Mat2D& Y, const Mat2D* init) {
  Mat2D out(X.rows(), Y.cols());
  for (Eigen::Index j = 0; j < Y.cols(); ++j)
    for (Eigen::Index i = 0; i < X.rows(); ++i) {
      long double s = init ? (long double)(*init)(i, j) : 0.0L;
      for (Eigen::Index k = 0; k < X.cols(); ++k) s += (long double)X(i, k) * Y(k, j);
      out(i, j) = (double)s;
    }
  return out;
}

static bool same_bytes(const Mat2D& a, const Mat2D& b) {
  return a.rows() == b.rows() && a.cols() == b.cols() &&
         std::memcmp(a.data(), b.data(), sizeof(double) * a.size()) == 0;
}

int main() {
  std::mt19937_64 rng(2026);
  const int threads[] = {1, 2, 3, 4, 7};
  // rows: below the split, the smallest split, odd heights, a last panel of
  // one row with 7 threads (1153 = 6 x 192 + 1), joined to the one before
  const Eigen::Index rows[] = {150, 192, 193, 287, 481, 1001, 1153};
  // depth: one kernel step, and more than one depth block (kc is ~250 doubles)
  const Eigen::Index depths[] = {37, 700};
  // columns of out: -k 10 with the default oversampling, a ragged count, one
  const Eigen::Index cols[] = {20, 7, 1};
  int split_cases = 0;
  for (const Eigen::Index n : rows)
    for (const Eigen::Index b : depths)
      for (const Eigen::Index l : cols) {
        const Mat2D X = random(n, b, rng), Y = random(b, l, rng), H0 = random(n, l, rng);
        const Mat2D want_set = reference(X, Y, nullptr), want_add = reference(X, Y, &H0);
        // a generous bound on the rounding of b products: a wrong row or panel is off by O(1)
        const double tol = 64 * 2.3e-16 * b * (X.cwiseAbs().maxCoeff() * Y.cwiseAbs().maxCoeff() + H0.cwiseAbs().maxCoeff());
        Mat2D first_set, first_add;
        for (const int t : threads) {
          omp_set_num_threads(t);
          const Eigen::Index panels = mul_X_Y_panels(n);
#ifndef EIGEN_USE_BLAS
          // split from 192 rows with 2 threads or more, at most one panel per thread
          const bool split = t > 1 && n >= 192;
          CHECK(split ? panels >= 2 && panels <= t : panels == 1, "n=%ld t=%d panels=%ld", (long)n, t, (long)panels);
          if (panels > 1) ++split_cases;
#else
          CHECK(panels == 1, "external BLAS: n=%ld t=%d panels=%ld", (long)n, t, (long)panels);
#endif
          // assignment over a matrix holding garbage, into a block of a larger
          // matrix as Halko does (out is an Eigen::Ref)
          Mat2D big = Mat2D::Constant(n + 5, l + 2, 1e300);
          mul_X_Y(X, Y, big.block(3, 1, n, l));
          const Mat2D set = big.block(3, 1, n, l);
          CHECK((set - want_set).cwiseAbs().maxCoeff() <= tol, "assign n=%ld b=%ld l=%ld t=%d err=%g", (long)n,
                (long)b, (long)l, t, (set - want_set).cwiseAbs().maxCoeff());
          // the frame around the block is untouched
          Mat2D frame = big;
          frame.block(3, 1, n, l).setConstant(1e300);
          CHECK((frame.array() == 1e300).all(), "assign n=%ld t=%d wrote outside its block", (long)n, t);
          // accumulation
          Mat2D add = H0;
          mul_X_Y(X, Y, add, true);
          CHECK((add - want_add).cwiseAbs().maxCoeff() <= tol, "add n=%ld b=%ld l=%ld t=%d err=%g", (long)n, (long)b,
                (long)l, t, (add - want_add).cwiseAbs().maxCoeff());
          // with two panels or more, the bytes of one thread (Eigen's own GEMM
          // blocks the depth differently once it threads, so not below that)
          if (t == 1) {
            first_set = set;
            first_add = add;
          } else if (n >= 192) {
            CHECK(same_bytes(set, first_set), "assign n=%ld b=%ld l=%ld: %d threads differ from 1", (long)n, (long)b,
                  (long)l, t);
            CHECK(same_bytes(add, first_add), "add n=%ld b=%ld l=%ld: %d threads differ from 1", (long)n, (long)b,
                  (long)l, t);
          }
        }
      }
  // inside an OpenMP region, one panel: the caller's threads are taken
  {
    omp_set_num_threads(4);
    Eigen::Index inside = 0;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
      inside = mul_X_Y_panels(1000);
    }
    CHECK(inside == 1, "panels inside a parallel region: %ld", (long)inside);
  }
#ifndef EIGEN_USE_BLAS
  CHECK(split_cases > 0, "the row split was never taken");
#endif
  if (failures) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
#ifndef EIGEN_USE_BLAS
  std::printf("mul_X_Y: assignment and accumulation match, split into row panels in %d cases, "
              "the same bytes with 1 to 7 threads from 192 rows\n",
              split_cases);
#else
  std::printf("mul_X_Y: external BLAS, no split; assignment and accumulation match\n");
#endif
  return 0;
}
