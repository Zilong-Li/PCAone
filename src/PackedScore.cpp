/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PackedScore.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "PackedScore.hpp"

#include <algorithm>
#include <cstring>
#include <vector>

namespace PCAone {

namespace {

// SIMD vectors of doubles by the compiler's vector extension (GCC, Clang):
// the sums of a kernel as an array of them, fully unrolled, which both keep
// in registers. g++ left an array of doubles in memory, and was 5x slower.
#if defined(__AVX__)
constexpr int VW = 4;
#else
constexpr int VW = 2;  // NEON, SSE2
#endif
typedef double vd __attribute__((vector_size(VW * sizeof(double))));
inline vd loadu(const double* p) {
  vd v;
  std::memcpy(&v, p, sizeof v);
  return v;
}
inline void storeu(double* p, vd v) { std::memcpy(p, &v, sizeof v); }

// SPREAD[x]: the 4 codes of the byte x, one to a byte, the first in the lowest.
// The byte of codes of 4 sites for 4 samples is then
// SPREAD[r0] | SPREAD[r1] << 2 | SPREAD[r2] << 4 | SPREAD[r3] << 6, the r the
// bytes of the 4 sites: byte q of it is sample q's code of site s at bits 2 s.
struct Spread {
  uint32_t v[256];
  Spread() {
    for (uint32_t x = 0; x < 256; ++x)
      v[x] = (x & 3) | ((x >> 2) & 3) << 8 | ((x >> 4) & 3) << 16 | ((x >> 6) & 3) << 24;
  }
};
const Spread SPREAD;

// the most groups of 4 sites whose tables are used together: their tables take
// up to about 256 KB, and each sample adds their rows into its sums in one go
constexpr int kMaxTile = 16;
#ifndef PCAONE_TABLE_KB
#define PCAONE_TABLE_KB 256
#endif
constexpr Eigen::Index kTableBytes = PCAONE_TABLE_KB << 10;

// the table of one group: T[x] = sum over its sites s of t_s[(x >> 2 s) & 3],
// t_s[c] = value[4 s + c] * D(j_s, :), padded to kp with zeros. Sites past j1
// add nothing.
void build_table(const double* value,
                 const Eigen::Ref<const Mat2D>& D,
                 Eigen::Index j,
                 Eigen::Index j1,
                 int kp,
                 double* site,  // 4 x 4 x kp
                 double* pair,  // 2 x 16 x kp
                 double* T) {   // 256 x kp
  const int K = (int)D.cols();
  std::fill(site, site + 16 * kp, 0.0);
  for (int s = 0; s < 4 && j + s < j1; ++s)
    for (int c = 0; c < 4; ++c) {
      const double v = value[4 * (j + s) + c];
      double* t = site + (4 * s + c) * kp;
      for (int k = 0; k < K; ++k) t[k] = v * D(j + s, k);
    }
  for (int h = 0; h < 2; ++h)  // the pairs of sites (0, 1) and (2, 3)
    for (int x = 0; x < 16; ++x) {
      const double* a = site + (4 * (2 * h) + (x & 3)) * kp;
      const double* b = site + (4 * (2 * h + 1) + (x >> 2)) * kp;
      double* p = pair + (16 * h + x) * kp;
      for (int k = 0; k < kp; ++k) p[k] = a[k] + b[k];
    }
  for (int x = 0; x < 256; ++x) {
    const double* a = pair + (x & 15) * kp;
    const double* b = pair + (16 + (x >> 4)) * kp;
    double* t = T + x * kp;
    for (int k = 0; k < kp; ++k) t[k] = a[k] + b[k];
  }
}

template <int KP>
void product_add(const PackedBlock& blk,
                 Eigen::Index j0,
                 Eigen::Index j1,
                 const double* value,
                 const Eigen::Ref<const Mat2D>& D,
                 Eigen::Index i0,
                 Eigen::Index i1,
                 double* acc,
                 Eigen::Index ld,
                 int kp) {
  const Eigen::Index ngroups = (j1 - j0 + 3) / 4;
  const int tile = (int)std::max<Eigen::Index>(1, std::min<Eigen::Index>(kMaxTile, kTableBytes / (256 * 8 * kp)));
  std::vector<double> tables((size_t)tile * 256 * kp), site(16 * kp), pair(32 * kp);
  const uint64_t b0 = (uint64_t)i0 / 4, b1 = ((uint64_t)i1 + 3) / 4;
  for (Eigen::Index g0 = 0; g0 < ngroups; g0 += tile) {
    const int nt = (int)std::min<Eigen::Index>(tile, ngroups - g0);
    const uint8_t* r[kMaxTile][4];
    for (int t = 0; t < nt; ++t) {
      const Eigen::Index j = j0 + 4 * (g0 + t);
      build_table(value, D, j, j1, kp, site.data(), pair.data(), tables.data() + (size_t)t * 256 * kp);
      // a site past j1 has a table of zeros: any record will do
      for (int s = 0; s < 4; ++s) r[t][s] = blk.site(std::min(j + s, j1 - 1));
    }
    const double* tab = tables.data();
    for (uint64_t b = b0; b < b1; ++b) {
      uint32_t w[kMaxTile];
      for (int t = 0; t < nt; ++t)
        w[t] = SPREAD.v[r[t][0][b]] | SPREAD.v[r[t][1][b]] << 2 | SPREAD.v[r[t][2][b]] << 4 |
               SPREAD.v[r[t][3][b]] << 6;
      const int nq = (int)std::min<uint64_t>(4, (uint64_t)i1 - 4 * b);
      for (int q = 0; q < nq; ++q) {
        double* __restrict a = acc + (4 * b + q - (uint64_t)i0) * ld;
        if constexpr (KP > 0) {
          // the sums of the sample stay in registers across the tables
          constexpr int NV = KP / VW;
          vd sum[NV];
#pragma GCC unroll 16
          for (int v = 0; v < NV; ++v) sum[v] = loadu(a + v * VW);
          for (int t = 0; t < nt; ++t) {
            const double* e = tab + ((size_t)t * 256 + ((w[t] >> (8 * q)) & 255)) * KP;
#pragma GCC unroll 16
            for (int v = 0; v < NV; ++v) sum[v] += loadu(e + v * VW);
          }
#pragma GCC unroll 16
          for (int v = 0; v < NV; ++v) storeu(a + v * VW, sum[v]);
        } else {
          for (int t = 0; t < nt; ++t) {
            const double* __restrict e = tab + ((size_t)t * 256 + ((w[t] >> (8 * q)) & 255)) * kp;
            for (int k = 0; k < kp; ++k) a[k] += e[k];
          }
        }
      }
    }
  }
}

template <int KP>
void gram(const double* base, const uint32_t* idx, uint32_t n, int, double sign, double* out) {
  constexpr int NV = KP / VW;
  for (int a = 0; a < KP; a += 2) {
    vd s0[NV], s1[NV];
#pragma GCC unroll 16
    for (int v = 0; v < NV; ++v) s0[v] = s1[v] = vd{};
    for (uint32_t e = 0; e < n; ++e) {
      const double* v = base + (size_t)idx[e] * KP;
      const double x0 = v[a], x1 = v[a + 1];
      // only the vectors that reach the upper triangle: the caller reads no more
#pragma GCC unroll 16
      for (int u = 0; u < NV; ++u) {
        if ((u + 1) * VW <= a) continue;
        const vd r = loadu(v + u * VW);
        s0[u] += x0 * r, s1[u] += x1 * r;
      }
    }
#pragma GCC unroll 16
    for (int v = 0; v < NV; ++v) {
      storeu(out + a * KP + v * VW, loadu(out + a * KP + v * VW) + sign * s0[v]);
      storeu(out + (a + 1) * KP + v * VW, loadu(out + (a + 1) * KP + v * VW) + sign * s1[v]);
    }
  }
}

}  // namespace

void add_gram(const double* base, const uint32_t* idx, uint32_t n, int kp, double sign, double* out) {
  switch (kp) {
    case 4: return gram<4>(base, idx, n, kp, sign, out);
    case 8: return gram<8>(base, idx, n, kp, sign, out);
    case 12: return gram<12>(base, idx, n, kp, sign, out);
    case 16: return gram<16>(base, idx, n, kp, sign, out);
    case 20: return gram<20>(base, idx, n, kp, sign, out);
    case 24: return gram<24>(base, idx, n, kp, sign, out);
    case 28: return gram<28>(base, idx, n, kp, sign, out);
    case 32: return gram<32>(base, idx, n, kp, sign, out);
    default: break;
  }
  // wider than the registers: v v' one site at a time
  for (uint32_t e = 0; e < n; ++e) {
    const double* v = base + (size_t)idx[e] * kp;
    for (int a = 0; a < kp; ++a)
      for (int k = 0; k < kp; ++k) out[a * kp + k] += sign * v[a] * v[k];
  }
}

void packed_product_add(const PackedBlock& blk,
                        Eigen::Index j0,
                        Eigen::Index j1,
                        const double* value,
                        const Eigen::Ref<const Mat2D>& D,
                        Eigen::Index i0,
                        Eigen::Index i1,
                        double* acc,
                        Eigen::Index ld) {
  if (j1 <= j0 || i1 <= i0) return;
  const int K = (int)D.cols();
  const int kp = (K + 3) / 4 * 4;  // the rows of the tables, a whole number of SIMD registers
  if (ld < kp) {
    // the tables are padded to kp, which acc has no room for: go through a padded copy
    Mat2D pad = Mat2D::Zero(D.rows(), kp);
    pad.leftCols(K) = D;
    std::vector<double> a((size_t)(i1 - i0) * kp, 0.0);
    packed_product_add(blk, j0, j1, value, pad, i0, i1, a.data(), kp);
    for (Eigen::Index i = 0; i < i1 - i0; ++i)
      for (int k = 0; k < K; ++k) acc[i * ld + k] += a[(size_t)i * kp + k];
    return;
  }
  switch (kp) {
    case 4: return product_add<4>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 8: return product_add<8>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 12: return product_add<12>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 16: return product_add<16>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 20: return product_add<20>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 24: return product_add<24>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 28: return product_add<28>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    case 32: return product_add<32>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
    default: return product_add<0>(blk, j0, j1, value, D, i0, i1, acc, ld, kp);
  }
}

}  // namespace PCAone
