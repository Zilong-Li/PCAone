// PackedGenotypes (PackedPca.cpp), the products of the in-core RSVD on the
// 2-bit calls: X' * Y and X * Y (and X * Y added) against the decoded matrix,
// for numbers of samples that do and do not fill the last byte of a record,
// missing calls imputed with 0 and with another value, widths that take one
// and several slabs of columns, sub-ranges of the sites, the records of a
// block given apart (out-of-core), and several numbers of threads.
#include <omp.h>

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <random>
#include <vector>

#define _DECLARE_TOOLBOX_HERE  // the logger and timer of the library, as Main.cpp
#include "../src/PackedPca.hpp"
#include "../src/Utils.hpp"

static int failures = 0;
#define CHECK(cond, ...)                                                            \
  do {                                                                              \
    if (!(cond)) {                                                                  \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed: ", __FILE__, __LINE__, #cond); \
      std::fprintf(stderr, __VA_ARGS__);                                            \
      std::fprintf(stderr, "\n");                                                   \
      ++failures;                                                                   \
    }                                                                               \
  } while (0)

static double maxrel(const Mat2D& a, const Mat2D& b) {
  return (a - b).cwiseAbs().maxCoeff() / std::max(1.0, b.cwiseAbs().maxCoeff());
}

int main() {
  std::mt19937_64 rng(2026);
  std::normal_distribution<double> nd;
  const double h[4] = {1, -9, 0.5, 0};  // BED2GENO, the missing code 1
  for (int N : {1, 3, 5, 64, 101, 1027, 4099}) {
    for (int M : {1, 7, 300, 2600}) {
      const uint64_t B = (N + 3) / 4;
      std::vector<uint8_t> rec(B * M);
      std::uniform_int_distribution<int> code(0, 3);
      std::bernoulli_distribution miss(0.05);
      Mat2D codes(N, M);
      for (int j = 0; j < M; ++j)
        for (int i = 0; i < N; ++i) {
          int c = code(rng);
          if (c == 1 && !miss(rng)) c = 0;  // few missing calls
          codes(i, j) = c;
          rec[j * B + i / 4] |= (uint8_t)(c << (2 * (i % 4)));
        }
      // the padding of the last byte is code 0, as in a BED
      auto buf = PCAone::PackedGenotypes::allocate(N, M);
      std::copy(rec.begin(), rec.end(), buf.get());
      PCAone::PackedGenotypes X(std::move(buf), N, M, h, 1);
      Mat1D a(M), b(M), m(M);
      for (int j = 0; j < M; ++j) a(j) = nd(rng), b(j) = nd(rng), m(j) = j % 3 ? 0.0 : nd(rng);
      X.set_sites(0, a, b, m);
      Mat2D D(N, M);
      for (int j = 0; j < M; ++j)
        for (int i = 0; i < N; ++i) {
          const int c = (int)codes(i, j);
          D(i, j) = c == 1 ? m(j) : a(j) + b(j) * h[c];
        }
      for (int K : {1, 6, 20, 33}) {
        Mat2D Y = Mat2D::NullaryExpr(N, K, [&]() { return nd(rng); });
        Mat2D Z = Mat2D::NullaryExpr(M, K, [&]() { return nd(rng); });
        for (int th : {1, 3, 8}) {
          omp_set_num_threads(th);
          Mat2D G(M, K);
          X.mul_Xt_Y(0, Y, G);
          Mat2D Gr = D.transpose() * Y;
          CHECK(maxrel(G, Gr) < 1e-12, "Xt_Y N=%d M=%d K=%d th=%d err=%g", N, M, K, th, maxrel(G, Gr));
          Mat2D H(N, K);
          X.mul_X_Y(0, Z, H);
          Mat2D Hr = D * Z;
          CHECK(maxrel(H, Hr) < 1e-12, "X_Y N=%d M=%d K=%d th=%d err=%g", N, M, K, th, maxrel(H, Hr));
          // a range of the sites, added
          const int j0 = M / 3, n = M - M / 3 - M / 5;
          if (n > 0) {
            Mat2D Gs(n, K);
            X.mul_Xt_Y(j0, Y, Gs);
            Mat2D Gsr = D.middleCols(j0, n).transpose() * Y;
            CHECK(maxrel(Gs, Gsr) < 1e-12, "Xt_Y range N=%d M=%d K=%d err=%g", N, M, K, maxrel(Gs, Gsr));
            Mat2D Ha = Hr;
            X.mul_X_Y(j0, Z.middleRows(j0, n), Ha, true);
            Mat2D Har = Hr + D.middleCols(j0, n) * Z.middleRows(j0, n);
            CHECK(maxrel(Ha, Har) < 1e-12, "X_Y add N=%d M=%d K=%d err=%g", N, M, K, maxrel(Ha, Har));
            // the same, the records of the sites given apart, as out-of-core
            PCAone::PackedGenotypes Xo(nullptr, N, M, h, 1);
            Xo.set_sites(0, a, b, m);
            const uint8_t* blk = rec.data() + (uint64_t)j0 * B;
            Mat2D Go(n, K), Ho = Hr;
            Xo.mul_Xt_Y(j0, Y, Go, blk);
            Xo.mul_X_Y(j0, Z.middleRows(j0, n), Ho, true, blk);
            CHECK(maxrel(Go, Gsr) < 1e-12 && maxrel(Ho, Har) < 1e-12, "block N=%d M=%d K=%d", N, M, K);
          }
        }
      }
    }
  }
  if (failures) std::fprintf(stderr, "%d failures\n", failures);
  else std::printf("test_packed_pca: all passed\n");
  return failures ? 1 : 0;
}
