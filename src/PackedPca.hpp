/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PackedPca.hpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 *
 * The genotype matrix of the RSVD (--svd 1 and 2) of a BED as its 2-bit calls:
 * a quarter of a byte a genotype where the decoded matrix takes 8 bytes, so
 * 32x less memory in-core, and the two products of every power iteration read
 * the calls straight into the sums by the tables of PackedScore.hpp, instead
 * of streaming a matrix of doubles from memory twice. Out-of-core, the blocks
 * of calls are not decoded either.
 ******************************************************************************/
#ifndef PCAONE_PACKEDPCA_HPP
#define PCAONE_PACKEDPCA_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "Common.hpp"

namespace PCAone {

class PackedGenotypes {
 public:
  using Index = Eigen::Index;

  // nsites records of (nsamples + 3) / 4 bytes, 4 samples a byte, the first in
  // the low bits; h[c] the genotype of code c (0 for the missing code), and
  // missing the code of a missing call. Without records (out-of-core), the
  // products are given the records of their sites.
  PackedGenotypes(std::unique_ptr<uint8_t[]> records, Index nsamples, Index nsites, const double h[4], int missing);
  // a buffer for nsites records, not initialized
  static std::unique_ptr<uint8_t[]> allocate(Index nsamples, Index nsites) {
    return std::unique_ptr<uint8_t[]>(new uint8_t[((uint64_t)nsamples + 3) / 4 * (uint64_t)nsites]);
  }

  Index rows() const { return N; }  // samples
  Index cols() const { return M; }  // sites
  uint64_t bytes() const { return B; }
  int missing_code() const { return miss; }
  const uint8_t* record(Index j) const { return rec.get() + (uint64_t)j * B; }

  // X(i, j) = a_j + b_j h[c_ij] for a called genotype, and m_j for a missing
  // call: the centred (a_j = -b_j f_j) and scaled (b_j) genotypes of the RSVD,
  // with m_j = 0 its missing calls imputed with the mean. For the sites
  // [j0, j0 + a.size()).
  void set_sites(Index j0, const Mat1D& a, const Mat1D& b, const Mat1D& m);

  // out = X' * Y for the sites [j0, j0 + out.rows()), Y nsamples x K; their
  // records at block, or else held
  void mul_Xt_Y(Index j0, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out, const uint8_t* block = nullptr);
  // out = X * Y, or out += X * Y with add, for the sites [j0, j0 + Y.rows())
  void mul_X_Y(Index j0,
               const Eigen::Ref<const Mat2D>& Y,
               Eigen::Ref<Mat2D> out,
               bool add = false,
               const uint8_t* block = nullptr);

 private:
  std::unique_ptr<uint8_t[]> rec;
  Index N, M;
  uint64_t B;
  double h[4];
  int miss;
  Mat1D sa, sb, sm;           // a_j, b_j, m_j of set_sites()
  std::vector<double> value;  // the value of each code at each site, 4 a site: of mul_X_Y
  std::vector<double> ws;     // the partial sums of the threads
  std::vector<double> Yr;     // the rows of Y of mul_Xt_Y, padded
};

}  // namespace PCAone

#endif  // PCAONE_PACKEDPCA_HPP
