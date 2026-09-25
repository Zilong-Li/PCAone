/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Halko.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Halko.hpp"

#include "Common.hpp"
#include "RSVD.hpp"
#include "Utils.hpp"

using namespace std;

void RsvdOpData::initOmg() {
  PortableRng rng(data->params.seed);  // same test matrix on every platform
  Omg.resize(cols(), size());
  for (Index j = 0; j < Omg.cols(); ++j)
    for (Index i = 0; i < Omg.rows(); ++i) Omg(i, j) = data->params.rand ? rng.normal() : rng.uniform(-1.0, 1.0);
  Omg2 = Omg;
}

// Start an EM update from the previous PCs instead of a random matrix: they
// are the dominant subspace of the matrix that the update only nudged, so the
// power iteration starts almost converged. The remaining columns stay random,
// to pick up what the update added.
void RsvdOpData::warmOmg() {
  if (!update || U.rows() != Omg.rows() || U.cols() == 0 || U.cols() > Omg.cols()) return;
  Omg.leftCols(U.cols()) = U;
  Eigen::HouseholderQR<Mat2D> qr(Omg);
  Omg = qr.householderQ() * Mat2D::Identity(Omg.rows(), Omg.cols());
  Omg2 = Omg;
}

void RsvdOpData::computeUSV(int p, double tol) {
  const Index nk{ranks()};
  const Index nrow{rows()};  // nsnps
  const Index ncol{cols()};  // nsamples
  const int bands = data->params.bands;
  const bool winsvd = data->params.svd_t == SvdType::PCAoneAlg2;
  // a warm-started winSVD runs every epoch with the full band from the start:
  // each epoch then ends with an exact Rayleigh-Ritz step, so it can stop after
  // two, where the doubling schedule needs log2(bands) + 1 epochs to get there
  pi_offset = (winsvd && update && U.rows() == ncol) ? (int)std::ceil(std::log2((double)bands)) : 0;
  Mat2D Upre, Ucur, H(ncol, size()), G(nrow, size()), B(size(), ncol), R(size(), size());
  double diff = 1.0;
  for (int pi = 0; pi <= p; ++pi) {
    computeGandH(G, H, pi);
    // G = QR. B = Q'X' follows from H = XG = XQR as B = R^-T H', so only R is
    // needed here, and Q only when V is formed at the end. G is overwritten by
    // the Householder reflectors, and rewritten by the next computeGandH().
    // (A second QR of Q used to follow; Q is orthonormal to machine precision
    // already, so it cost as much as the first and changed nothing but signs.)
    Eigen::HouseholderQR<Eigen::Ref<Mat2D>> qr(G);
    R = qr.matrixQR().topRows(size()).triangularView<Eigen::Upper>();
    B.noalias() = R.transpose().fullPivHouseholderQr().solve(H.transpose());
    Eigen::JacobiSVD<Mat2D> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Ucur = svd.matrixV().leftCols(nk);
    if (pi > 0) {
      if (data->params.mev)
        diff = 1 - mev(Ucur, Upre);
      else
        diff = minSSE(Ucur, Upre).sum() / Upre.cols();
      if (data->params.verbose && !data->params.missme)
        cao.print(tick.date(), "running of epoch =", pi, ", diff =", diff);
      if (diff < tol || pi == p) {
        if (winsvd && std::pow(2, pi + pi_offset) < bands) {
          cao.print("PCAone winSVD converged but continues running to get S and V.");
          p = std::log2(bands);
        } else {
          U = Ucur;
          Mat2D W = Mat2D::Zero(nrow, nk);
          W.topRows(size()) = svd.matrixU().leftCols(nk);
          V.noalias() = qr.householderQ() * W;  // Q * W
          S = svd.singularValues().head(nk);
          if (data->params.verbose && !data->params.missme) cao.print(tick.date(), "stops at epoch =", pi + 1);
          if (!(diff < tol) && !data->params.missme)
            cao.warn("the RSVD reached --maxp " + std::to_string(p) + " epochs before converging (diff = " +
                     std::to_string(diff) + ", --tol-rsvd = " + std::to_string(tol) +
                     "). the trailing PCs may be inaccurate; consider a larger --maxp, or --svd 0");
          break;
        }
      } else {
        Upre = Ucur;
      }
    } else {
      Upre = Ucur;
    }
  }
  U = Ucur;
}

void NormalRsvdOpData::computeGandH(Mat2D& G, Mat2D& H, int pi) {
  if (!data->snpmajor) {
    cao.error("only work with snp major input data now.");
  }

  // reset omg to random, or to the previous PCs for an EM update
  if (pi == 0) {
    initOmg();
    warmOmg();
  }

  if (!data->params.out_of_core) {
    if (pi == 0) {
      if (update) {
        data->fit_with_pi(U, S, V.transpose());
      }
      if (standardize) {
        if (data->params.pcangsd) {
          data->pcangsd_standardize_E(U, S, V.transpose());
        } else {
          data->standardize_E();
        }
      }
    }
    if (pi > 0) {
      Eigen::HouseholderQR<Eigen::Ref<Mat2D>> qr(H);
      Omg.noalias() = qr.householderQ() * Mat2D::Identity(cols(), size);  // hold H in Omega
      PCAone::flipOmg(Omg2, Omg);
    }
    mul_Xt_Y(data->G, Omg, G);
    H.noalias() = data->G * G;
    return;
  }

  // for block version
  // data->G is always nsamples x nsnps;
  // for nsnps > nsamples
  if (pi > 0) {
    Eigen::HouseholderQR<Eigen::Ref<Mat2D>> qr(H);
    Omg.noalias() = qr.householderQ() * Mat2D::Identity(cols(), size);
  }
  H = Mat2D::Zero(cols(), size);
  data->check_file_offset_first_var();
  for (uint i = 0; i < data->nblocks; ++i) {
    start_idx = data->start[i];
    stop_idx = data->stop[i];
    actual_block_size = stop_idx - start_idx + 1;
    tick.clock();
    if (!update) {
      data->read_block_initial(start_idx, stop_idx, standardize);
    } else {
      data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
    }
    data->readtime += tick.reltime();
    mul_Xt_Y(data->G, Omg, G.middleRows(start_idx, actual_block_size));
    H.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
  }
}

void FancyRsvdOpData::computeGandH(Mat2D& G, Mat2D& H, int pi) {
  // check size of G and H first;
  if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
    cao.error("the size of G or H doesn't match with each other.");
  }
  if (pi == 0) {
    initOmg();
    warmOmg();
  }
  const int pe = pi + pi_offset;  // position in the band schedule
  if (std::pow(2, pe) >= data->params.bands) {
    // init H1, H2 to zero
    H1.setZero();
    H2.setZero();
  }
  if (!data->params.out_of_core) {
    if (pi == 0) {
      if (update) {
        data->fit_with_pi(U, S, V.transpose());
      }
      if (standardize) {
        if (data->params.pcangsd) {
          data->pcangsd_standardize_E(U, S, V.transpose());
        } else {
          data->standardize_E();
        }
      }
      bandsize = pi_offset > 0 ? data->params.bands : 1;
      // blocksize: how many snps in each block
      blocksize = (unsigned int)ceil((double)data->nsnps / data->params.bands);
      if (blocksize < data->params.bands)
        cao.warn("block size < window size. please consider the IRAM method with --svd 0");
      // Keep the same column order across EM iterations and final scaling.
      if (data->params.perm && !data->in_core_permuted) {
        cao.print(tick.date(), "permuting data matrix by columns in place");
        PCAone::permute_matrix(data->G, data->perm, data->params.seed);
        data->in_core_permuted = true;
      }
    }
    // bandsize: how many blocks in each band, 2, 4, 8, 16, 32, 64, ...
    bandsize = fmin(bandsize * 2, data->params.bands);
    // b: the index of current block
    for (uint b = 0, i = 1; b < data->params.bands; ++b, ++i) {
      // half-open [start_idx, start_idx + actual_block_size). With fewer than
      // bands^2 sites the trailing blocks are empty: ceil(M / bands) * b can
      // pass M, and the old inclusive stop_idx then underflowed the unsigned
      // block size, so the next line asked for ~2^64 rows (std::bad_alloc).
      // An empty block contributes nothing but keeps the band schedule intact.
      start_idx = std::min<uint64>((uint64)b * blocksize, data->nsnps);
      stop_idx = std::min<uint64>((uint64)(b + 1) * blocksize, data->nsnps);
      actual_block_size = stop_idx - start_idx;
      mul_Xt_Y(data->G.middleCols(start_idx, actual_block_size), Omg, G.middleRows(start_idx, actual_block_size));

      if (i <= bandsize / 2) {
        // continues to add in data based on current band
        H1.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
      } else {
        H2.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
      }

      // use the first quarter band of succesive iteration (H1)
      // for extra power iteration updates with the last used band (H2)
      const bool adjacent = (pe > 0 && (b + 1) == std::pow(2, pe - 1) && std::pow(2, pe) < data->params.bands);
      if ((b + 1) < bandsize && !adjacent) continue;

      // add up H and update Omg
      if (!((i == bandsize) || (i == bandsize / 2) || adjacent)) continue;
      H = H1 + H2;
      Eigen::HouseholderQR<Mat2D> qr(H);
      Omg.noalias() = qr.householderQ() * Mat2D::Identity(cols(), size);
      PCAone::flipOmg(Omg2, Omg);
      if (i == bandsize) {
        H1.setZero();
        i = 0;
      } else {
        H2.setZero();
      }
    }

    return;
  }

  // out-of-core implementation
  if (pi == 0) bandsize = pi_offset > 0 ? data->nblocks : data->bandFactor;
  data->check_file_offset_first_var();
  // band : 2, 4, 8, 16, 32, 64
  bandsize = fmin(bandsize * 2, data->nblocks);
  for (uint b = 0, i = 1; b < data->nblocks; ++b, ++i) {
    start_idx = data->start[b];
    stop_idx = data->stop[b];
    actual_block_size = stop_idx - start_idx + 1;
    tick.clock();
    if (!update) {
      data->read_block_initial(start_idx, stop_idx, standardize);
    } else {
      data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
    }
    data->readtime += tick.reltime();
    mul_Xt_Y(data->G, Omg, G.middleRows(start_idx, actual_block_size));

    if (i <= bandsize / 2) {
      H1.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
    } else {
      H2.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
    }

    const bool adjacent =
        (pe > 0 && (b + 1) == std::pow(2, pe - 1) * data->bandFactor && std::pow(2, pe) < data->params.bands);
    if ((b + 1) < bandsize && !adjacent) continue;

    // cao.print("i:", i, ",j:", j, ",bandsize:", bandsize, ",pi:", pi);
    // add up H and update Omg
    if (!((i == bandsize) || (i == bandsize / 2) || adjacent)) continue;
    H = H1 + H2;
    Eigen::HouseholderQR<Mat2D> qr(H);
    Omg.noalias() = qr.householderQ() * Mat2D::Identity(cols(), size);
    PCAone::flipOmg(Omg2, Omg);
    if (i == bandsize) {
      H1.setZero();
      i = 0;
    } else {
      H2.setZero();
    }
  }
}

void run_pca_with_halko(Data* data, const Param& params) {
  Mat2D Vpre;
  RsvdOpData* rsvd;
  if (params.svd_t == SvdType::PCAoneAlg2) {
    cao.print(tick.date(), "initialize window-based RSVD (winSVD) with",
              params.out_of_core ? "out-of-core" : "in-core");
    rsvd = new FancyRsvdOpData(data, params.k, params.oversamples);
  } else {
    cao.print(tick.date(), "initialize single-pass RSVD (sSVD) with", params.out_of_core ? "out-of-core" : "in-core");
    rsvd = new NormalRsvdOpData(data, params.k, params.oversamples);
  }
  if (!params.missme) {
    if (params.genetic) {
      rsvd->setFlags(false, true);
    } else {
      rsvd->setFlags(false, false);
    }
    rsvd->computeUSV(params.maxp, params.tol);
    // the sign of each PC is arbitrary; fix it as --svd 0 and 3 do, so that it
    // no longer depends on the internals of the QR and the SVD of B
    flip_UV(rsvd->U, rsvd->V);
  } else {
    if (data->p_miss == 0.0 && !params.out_of_core) cao.warn("there is no missing values");
    // for EM iteration
    rsvd->setFlags(false, false);
    rsvd->computeUSV(params.maxp, params.tol);
    flip_UV(rsvd->U, rsvd->V, false);
    double diff = 1.0;
    cao.print(tick.date(), "run EM-PCA. maxiter =", params.maxiter);
    for (uint i = 0; i < params.maxiter; ++i) {
      rsvd->setFlags(true, false);
      Vpre = rsvd->V;
      rsvd->computeUSV(params.maxp, params.tol);
      flip_UV(rsvd->U, rsvd->V, false);
      if (params.mev)
        diff = 1.0 - mev(rsvd->V, Vpre);
      else
        diff = minSSE(rsvd->V, Vpre).sum() / Vpre.cols();
      cao.print(tick.date(), "individual allele frequencies estimated iter =", i + 1, ", diff =", diff);
      if (diff < params.tolem) {
        cao.print(tick.date(), "come to convergence!");
        break;
      }
    }
    if (params.maxiter > 0 && !(diff < params.tolem)) warn_em_not_converged(params, diff);

    if (params.emu) {
      cao.print(tick.date(), "standardize the final matrix for EMU");
      rsvd->setFlags(true, true);
      rsvd->computeUSV(params.maxp, params.tol);
      flip_UV(rsvd->U, rsvd->V, false);
    }

    if (params.pcangsd && (params.file_t == FileType::BEAGLE)) {
      cao.print(tick.date(), "estimate GRM for pcangsd");
      data->pcangsd_standardize_E(rsvd->U, rsvd->S, rsvd->V.transpose());
      // TODO: use matrix-free method e.g Arnoldi to decompose the cov
      write_pcangsd_cov(data->G, data->Dc, data->nsnps, params);
    }
  }
  // output PI
  data->set_svd_transform(rsvd->standardize);
  data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->S, rsvd->U, rsvd->V);

  delete rsvd;

  cao.print(tick.date(), "PCAone - Randomized SVD done");
  return;
}
