/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PackedScore.hpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_PACKEDSCORE_HPP
#define PCAONE_PACKEDSCORE_HPP

#include <cstdint>
#include <cstring>

#include "Common.hpp"

namespace PCAone {

// A block of 2-bit calls: nsites records of `bytes` bytes, 4 samples a byte,
// the first sample in the low bits. That is the layout of a BED record and of
// PGEN's genovec, which differ only in what each code means.
struct PackedBlock {
  const uint8_t* rec = nullptr;
  uint64_t bytes = 0;      // per record, (N + 3) / 4
  Eigen::Index nsites = 0;
  int missing = 1;         // the code of a missing call: 1 in BED, 3 in PGEN
  const uint8_t* site(Eigen::Index j) const { return rec + (uint64_t)j * bytes; }
};

// acc(i - i0, :) += sum over the sites j in [j0, j1) of value[4 j + c_ij] * D(j, :)
// for the samples i in [i0, i1), with c_ij the code of sample i at site j. acc
// is row-major with a row of ld >= D.cols() doubles per sample, and i0 a
// multiple of 4.
//
// This is the product of the block, decoded and centred, with D, without
// decoding it: the sites are taken 4 at a time, and the 256 sums of their
// rows of D for the 4^4 codes the 4 sites can have are tabulated once, so each
// sample adds one K-vector per 4 sites, read by its byte of codes. That is a
// quarter of the arithmetic of the product and reads 2 bits per genotype
// instead of 8 bytes. value is 0 for a missing call.
void packed_product_add(const PackedBlock& blk,
                        Eigen::Index j0,
                        Eigen::Index j1,
                        const double* value,
                        const Eigen::Ref<const Mat2D>& D,
                        Eigen::Index i0,
                        Eigen::Index i1,
                        double* acc,
                        Eigen::Index ld);

// out (kp x kp, row-major) += sign * sum over e < n of v_e v_e', with v_e the
// row of kp doubles at base + idx[e] * kp, in its upper triangle (out(a, c),
// c >= a; the rest of out is not to be read). For the missing calls of
// --project 2: the sum of v_j v_j' over the sites j a sample misses. Two rows
// of out at a time stay in registers while the rows v_e go by.
void add_gram(const double* base, const uint32_t* idx, uint32_t n, int kp, double sign, double* out);

// One bit per sample (at the even bits) of the 32 samples from sample 32 * w
// on whose code at site j is blk.missing, up to sample n (exclusive).
inline uint64_t missing_word(const PackedBlock& blk, Eigen::Index j, uint64_t w, uint64_t n) {
  const uint8_t* r = blk.site(j) + w * 8;
  const uint64_t left = blk.bytes - w * 8;
  uint64_t x = 0;
  if (left >= 8)
    std::memcpy(&x, r, 8);  // one load
  else
    std::memcpy(&x, r, left);
  const uint64_t y = x ^ (0x5555555555555555ULL * (uint64_t)blk.missing);  // a missing code is 00
  uint64_t m = ~(y | (y >> 1)) & 0x5555555555555555ULL;
  const uint64_t ns = n - w * 32;  // samples of this word to keep
  if (ns < 32) m &= (1ULL << (2 * ns)) - 1;
  return m;
}

// f(i, j) for each missing call of the samples [i0, i1) at the sites
// [0, nsites), site by site; i0 a multiple of 32.
template <typename F>
void for_each_missing(const PackedBlock& blk, Eigen::Index i0, Eigen::Index i1, F&& f) {
  const uint64_t w0 = (uint64_t)i0 / 32, w1 = ((uint64_t)i1 + 31) / 32;
  for (Eigen::Index j = 0; j < blk.nsites; ++j)
    for (uint64_t w = w0; w < w1; ++w)
      for (uint64_t m = missing_word(blk, j, w, (uint64_t)i1); m; m &= m - 1)
        f((Eigen::Index)(w * 32 + (__builtin_ctzll(m) >> 1)), j);
}

// the code of sample i at site j
inline int packed_code(const PackedBlock& blk, Eigen::Index j, Eigen::Index i) {
  return (blk.site(j)[i >> 2] >> (2 * (i & 3))) & 3;
}

}  // namespace PCAone

#endif  // PCAONE_PACKEDSCORE_HPP
