/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileUSV.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileUSV.hpp"

#include <cmath>

using namespace std;

// Recover the individual allele frequency pi on the 0..1 scale that F and
// BED2GENO use, from the reconstruction U*S*V' of whatever matrix the reference
// PCA decomposed. All three cases are pi = U*S*V' * inv_scale(f) + f; only the
// factor differs, and FileUSV::inv_scale() picks it from the transform recorded
// in .sigvals:
//
//   standardized 0..1  sd/sqrt(ploidy)   the default genetic path; standardize_E()
//                                        scaled column s by sqrt(ploidy)/sd, and
//                                        skipping the inverse inflates every
//                                        deviation by sqrt(ploidy)/sd, worst
//                                        where sd is smallest i.e. at rare variants
//   centred 0..1       1                 --missme without --emu: never standardized
//   centred dosages    0.5               pcangsd: fit_with_pi() builds dosage - 2f
static inline double usv_to_pi(double usv, double f, double inv_scale) {
  return fmin(fmax(usv * inv_scale + f, 1e-4), 1.0 - 1e-4);
}

// Decide whether the reconstruction can be mapped back to allele frequencies at
// all, and say which assumption is being made when the reference run predates
// the recording of it.
void FileUSV::check_transform() {
  if (!usv.known) {
    // Written before .sigvals carried the transform. Guess from the input type:
    // a GL reference run is the pcangsd path, anything else took the defaults.
    usv.gscale = (params.file_t == FileType::BEAGLE) ? 2 : 1;
    usv.scale = SCALE_STANDARDIZE_GENETIC;
    usv.ploidy = params.ploidy;
    cao.warn(params.fileS,
             " predates the recording of the PCA scaling, so it is assumed to be the default (gscale=", usv.gscale,
             ", scale=-9, ploidy=", usv.ploidy,
             "). rerun the reference PCA to record it. the assumption is wrong if that run used -C/--scale, or "
             "--missme without --emu.");
    return;
  }
  if (usv.scale > 0)
    cao.error("the reference PCA used -C/--scale ", usv.scale,
              ", which does not map back to allele frequencies. --inbreed needs a reference run with the default "
              "genetic standardisation or none at all");
  if (usv.gscale == 2 && usv.scale == SCALE_STANDARDIZE_GENETIC)
    cao.error("the reference PCA standardised a dosage-scale matrix (gscale=2, scale=-9); --inbreed cannot invert "
              "that combination. rerun the reference PCA with --svd 1 or 2");
  cao.print(tick.date(), "USV transform: scale =", usv.scale, ", ploidy =", usv.ploidy, ", gscale =", usv.gscale);
}

void FileUSV::read_all() {
  // S holds every singular value of the reference, U and V only the first K
  // columns: with -k below the reference's k the product did not conform
  G = U * S.head(K).asDiagonal() * V.transpose();
  if (params.inbreed) {
    // get \PI and store it in G. the factor depends only on the site, so it is
    // hoisted out of the sample loop: gcc will not do it itself, because F and
    // G are unrelated pointers it cannot prove do not alias.
#pragma omp parallel for
    for (int i = 0; i < G.cols(); i++) {
      const double f = F(i);
      const double s = inv_scale(f);
      for (int j = 0; j < G.rows(); j++) G(j, i) = usv_to_pi(G(j, i), f, s);
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
  // U*S*V' of the block, a panel of sites per thread, each mapped to pi while it
  // is still in cache; it was a scalar triple loop, O(N K) per site on every
  // pass of every SQUAREM iteration
  const Mat2D US = U * S.head(K).asDiagonal();
  const Eigen::Index bs = std::max<Eigen::Index>(1, (Eigen::Index(1) << 15) / std::max<Eigen::Index>(1, nsamples));
  const Eigen::Index nb = ((Eigen::Index)actual_block_size + bs - 1) / bs;
#pragma omp parallel for schedule(static)
  for (Eigen::Index b = 0; b < nb; ++b) {
    const Eigen::Index c = b * bs, w = std::min<Eigen::Index>(bs, actual_block_size - c);
    G.middleCols(c, w).noalias() = US * V.middleRows(start_idx + c, w).transpose();
    if (!params.inbreed) continue;
    for (Eigen::Index i = c; i < c + w; ++i) {
      // hoisted out of the sample loop, same reasoning as read_all()
      const double f = F(start_idx + i);
      const double s = inv_scale(f);
      for (uint j = 0; j < nsamples; j++) G(j, i) = usv_to_pi(G(j, i), f, s);  // map to domain, as read_all()
    }
  }
}
