/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PackedPca.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "PackedPca.hpp"

#include <algorithm>
#include <cstring>

#include "PackedScore.hpp"

namespace PCAone {

namespace {
// the columns of a product are taken up to 32 at a time: the widest sums the
// kernels keep in registers
constexpr Eigen::Index kSlab = 32;
// The products are split in about this many tasks, whatever the number of
// threads: the partial sums are added in the same order, so the result does
// not depend on -n.
constexpr Eigen::Index kTasks = 64;
inline int pad4(Eigen::Index k) { return (int)((k + 3) / 4 * 4); }
inline Eigen::Index clampi(Eigen::Index x, Eigen::Index lo, Eigen::Index hi) { return std::max(lo, std::min(x, hi)); }
}  // namespace

PackedGenotypes::PackedGenotypes(
    std::unique_ptr<uint8_t[]> records, Index nsamples, Index nsites, const double h_[4], int missing)
    : rec(std::move(records)),
      N(nsamples),
      M(nsites),
      B(((uint64_t)nsamples + 3) / 4),
      miss(missing) {
  for (int c = 0; c < 4; ++c) h[c] = c == missing ? 0.0 : h_[c];
  sa = sb = sm = Mat1D::Zero(M);
  value.assign((size_t)4 * M, 0.0);
}

void PackedGenotypes::set_sites(Index j0, const Mat1D& a, const Mat1D& b, const Mat1D& m) {
  const Index n = a.size();
  sa.segment(j0, n) = a, sb.segment(j0, n) = b, sm.segment(j0, n) = m;
#pragma omp parallel for schedule(static)
  for (Index j = 0; j < n; ++j)
    for (int c = 0; c < 4; ++c) value[4 * (j0 + j) + c] = c == miss ? m(j) : a(j) + b(j) * h[c];
}

void PackedGenotypes::mul_Xt_Y(Index j0, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out, const uint8_t* block) {
  const Index n = out.rows(), K = Y.cols();
  if (n == 0 || K == 0) return;
  const PackedBlock blk{block ? block : record(j0), B, n, miss};
  // sites in chunks of about 2048: each tabulates the sums of its samples once
  // for them all, and the sums of its sites, which every table adds to, stay in
  // the cache. With fewer chunks than kTasks, the samples are split in panels
  // of at least 16 bytes as well, each panel adding into its own sums.
  const Index C = std::max<Index>(1, n / 2048), P = clampi((kTasks + C - 1) / C, 1, std::max<Index>(1, (Index)B / 16));
  for (Index c0 = 0; c0 < K; c0 += kSlab) {
    const Index kc = std::min(kSlab, K - c0);
    const int kp = pad4(kc);
    // the rows of Y, padded to kp, and zero past the last sample
    Yr.resize((size_t)4 * B * kp);
#pragma omp parallel for schedule(static)
    for (Index i = 0; i < 4 * (Index)B; ++i)
      for (Index k = 0; k < kp; ++k) Yr[(size_t)i * kp + k] = i < N && k < kc ? Y(i, c0 + k) : 0.0;
    ws.resize((size_t)P * n * kp);
#pragma omp parallel for schedule(dynamic, 1)
    for (Index t = 0; t < C * P; ++t) {
      const Index c = t / P, p = t % P;
      const Index s0 = n * c / C, s1 = n * (c + 1) / C;
      const uint64_t g0 = B * p / P, g1 = B * (p + 1) / P;
      double* acc = ws.data() + ((size_t)p * n + s0) * kp;
      std::fill(acc, acc + (s1 - s0) * kp, 0.0);
      packed_tproduct_add(blk, s0, s1, h, Yr.data() + g0 * 4 * kp, kp, g0, g1, acc);
    }
    // the sum of Y over the samples, for the centring
    Mat1D tot = Y.middleCols(c0, kc).colwise().sum().transpose();
#pragma omp parallel
    {
      Mat1D mis(kc);
#pragma omp for schedule(static)
      for (Index j = 0; j < n; ++j) {
        // the sum of Y over the samples whose call is missing at site j
        mis.setZero();
        bool any = false;
        for (uint64_t w = 0; w < (B + 7) / 8; ++w)
          for (uint64_t mw = missing_word(blk, j, w, (uint64_t)N); mw; mw &= mw - 1) {
            const Index i = (Index)(w * 32 + (__builtin_ctzll(mw) >> 1));
            mis += Y.row(i).segment(c0, kc).transpose();
            any = true;
          }
        const double a = sa(j0 + j), b = sb(j0 + j), m = sm(j0 + j);
        for (Index k = 0; k < kc; ++k) {
          double y = 0.0;
          for (Index p = 0; p < P; ++p) y += ws[((size_t)p * n + j) * kp + k];
          out(j, c0 + k) = a * tot(k) + b * y + (any ? (m - a) * mis(k) : 0.0);
        }
      }
    }
  }
}

void PackedGenotypes::mul_X_Y(
    Index j0, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out, bool add, const uint8_t* block) {
  const Index n = Y.rows(), K = Y.cols();
  if (!add) out.setZero();
  if (n == 0 || K == 0) return;
  const PackedBlock blk{block ? block : record(j0), B, n, miss};
  // samples in panels of at least 2048, as each tabulates the sums of its
  // sites once for them all, and the sites in chunks of at least 64 for about
  // kTasks tasks; each chunk adds into its own sums
  const Index P = std::max<Index>(1, N / 2048), C = clampi((kTasks + P - 1) / P, 1, std::max<Index>(1, n / 64));
  for (Index c0 = 0; c0 < K; c0 += kSlab) {
    const Index kc = std::min(kSlab, K - c0);
    const int kp = pad4(kc);
    // the rows of a panel start at a byte of the calls
    auto row0 = [&](Index p) { return (Index)(B * p / P) * 4; };
    ws.resize((size_t)C * 4 * B * kp);
#pragma omp parallel for schedule(dynamic, 1)
    for (Index t = 0; t < C * P; ++t) {
      const Index c = t / P, p = t % P;
      const Index s0 = n * c / C, s1 = n * (c + 1) / C;
      const Index i0 = row0(p), i1 = std::min(N, row0(p + 1));
      double* acc = ws.data() + ((size_t)c * 4 * B + i0) * kp;
      std::fill(acc, acc + (i1 - i0) * kp, 0.0);
      packed_product_add(blk, s0, s1, value.data() + 4 * j0, Y.middleCols(c0, kc), i0, i1, acc, kp);
    }
#pragma omp parallel for schedule(static)
    for (Index i = 0; i < N; ++i)
      for (Index k = 0; k < kc; ++k) {
        double s = 0.0;
        for (Index c = 0; c < C; ++c) s += ws[((size_t)c * 4 * B + i) * kp + k];
        out(i, c0 + k) += s;
      }
  }
}

}  // namespace PCAone
