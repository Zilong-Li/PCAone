/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileUSV.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileUSV.hpp"

using namespace std;

// get \PI and store it in G
void FileUSV::read_all() {
  G = U * S.asDiagonal() * V.transpose();
  if (params.pcangsd) {
    for (int i = 0; i < F.size(); i++) {
      G.col(i) = (G.array().col(i) - F(i)) / 2.0;
    }
  }
}

// get a block
// NOTE: sanity check if blocks are continuous
void FileUSV::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < params.blocksize || (actual_block_size < params.blocksize)) {
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
      // NOTE:  map to domain(0,1)?
      if (params.pcangsd) {
        G(j, i) = (G(j, i) + 2.0 * F(snp_idx)) / 2.0;
        G(j, i) = fmin(fmax(G(j, i), 1e-4), 1.0 - 1e-4);
      }
    }
  }
}
