/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileUSV.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileUSV.hpp"

#include <cmath>

using namespace std;

void FileUSV::read_all() {
  G = U * S.asDiagonal() * V.transpose();
  if (params.inbreed) {
// get \PI and store it in G
#pragma omp parallel for
    for (int i = 0; i < G.cols(); i++) {
      for (int j = 0; j < G.rows(); j++) {
        // FIX: U*S*V' is on the STANDARDISED scale -- Data::standardize_E
        // divides each column by sd = sqrt(f(1-f)) and multiplies by
        // sqrt(ploidy) -- so it must be returned to the genotype scale before
        // 2f is added.  Otherwise every deviation is inflated by
        // sqrt(ploidy)/sd, worst where sd is smallest, i.e. at rare variants.
        {
          const double sd = std::sqrt(F(i) * (1.0 - F(i)));
          const double g = (sd > 1e-9) ? G(j, i) * sd * std::sqrt(2.0) : 0.0;
          G(j, i) = (g + 2.0 * F(i)) * 0.5;
        }
        G(j, i) = fmin(fmax(G(j, i), 1e-4), 1.0 - 1e-4);
      }
    }
  }
}

// get a block
// NOTE: sanity check if blocks are continuous
void FileUSV::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = Mat2D::Zero(nsamples, actual_block_size);
  }
#pragma omp parallel for
  for (uint i = 0; i < actual_block_size; ++i) {
    uint64 snp_idx = start_idx + i;
    for (uint j = 0; j < nsamples; j++) {
      G(j, i) = 0.0;
      for (int k = 0; k < K; ++k) {
        G(j, i) += U(j, k) * S(k) * V(snp_idx, k);
      }
      if (params.inbreed) {
        //  map to domain -- same rescaling as read_all()
        const double sd = std::sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
        const double g = (sd > 1e-9) ? G(j, i) * sd * std::sqrt(2.0) : 0.0;
        G(j, i) = (g + 2.0 * F(snp_idx)) * 0.5;
        G(j, i) = fmin(fmax(G(j, i), 1e-4), 1.0 - 1e-4);
      }
    }
  }
}
