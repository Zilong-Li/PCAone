/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileUSV.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileUSV.hpp"

#include <cmath>

using namespace std;

// Undo the standardisation that the PCA applied, recovering the individual
// allele frequency pi on the 0..1 scale that F and BED2GENO use.
//
// U*S*V' reconstructs the matrix that was decomposed. For the default genetic
// path that is the STANDARDISED one: Data::standardize_E() scales column s by
// sqrt(ploidy)/sd with sd = sqrt(f(1-f)). The inverse is therefore
//
//     pi = U*S*V' * sd/sqrt(ploidy) + f
//
// Skipping the sd/sqrt(ploidy) factor inflates every deviation by
// sqrt(ploidy)/sd, worst where sd is smallest, i.e. at rare variants.
//
// standardize_E() leaves a column untouched when sd <= VAR_TOL, so the inverse
// leaves it untouched too -- hence the factor is 1, not 0.
static inline double usv_to_pi(double usv, double f, double inv_scale) {
  return fmin(fmax(usv * inv_scale + f, 1e-4), 1.0 - 1e-4);
}

// sd/sqrt(ploidy), or 1 where standardize_E() did not scale the column
static inline double usv_inv_scale(double f, double rploidy) {
  const double sd = std::sqrt(f * (1.0 - f));
  return (sd > VAR_TOL) ? sd / rploidy : 1.0;
}

void FileUSV::read_all() {
  G = U * S.asDiagonal() * V.transpose();
  if (params.inbreed) {
    // get \PI and store it in G. sd depends only on the site, so it is hoisted
    // out of the sample loop: gcc will not do it itself, because F and G are
    // unrelated pointers it cannot prove do not alias.
    const double rploidy = std::sqrt((double)params.ploidy);
#pragma omp parallel for
    for (int i = 0; i < G.cols(); i++) {
      const double f = F(i);
      const double inv_scale = usv_inv_scale(f, rploidy);
      for (int j = 0; j < G.rows(); j++) G(j, i) = usv_to_pi(G(j, i), f, inv_scale);
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
  const double rploidy = std::sqrt((double)params.ploidy);
#pragma omp parallel for
  for (uint i = 0; i < actual_block_size; ++i) {
    uint64 snp_idx = start_idx + i;
    // hoisted out of the sample loop, same reasoning as read_all()
    const double f = F(snp_idx);
    const double inv_scale = params.inbreed ? usv_inv_scale(f, rploidy) : 0.0;
    for (uint j = 0; j < nsamples; j++) {
      G(j, i) = 0.0;
      for (int k = 0; k < K; ++k) {
        G(j, i) += U(j, k) * S(k) * V(snp_idx, k);
      }
      //  map to domain -- same rescaling as read_all()
      if (params.inbreed) G(j, i) = usv_to_pi(G(j, i), f, inv_scale);
    }
  }
}
