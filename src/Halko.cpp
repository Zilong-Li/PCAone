/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Halko.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Halko.hpp"

#include "Common.hpp"
#include "Utils.hpp"

using namespace std;

void RsvdOpData::initOmg() {
  auto rng = std::default_random_engine{};
  if (data->params.rand)
    Omg = PCAone::StandardNormalRandom<Mat2D, std::default_random_engine>(cols(), size(), rng);
  else
    Omg = PCAone::UniformRandom<Mat2D, std::default_random_engine>(cols(), size(), rng);
  Omg2 = Omg;
}

Mat2D RsvdOpData::computeU(const Mat2D& G, const Mat2D& H) {
  const Index nrow{G.rows()};
  const Index nk{ranks()};
  Mat2D R(size(), size()), Rt(size(), size()), Gt(nrow, G.cols());
  Eigen::HouseholderQR<Mat2D> qr(G);
  R.noalias() = Mat2D::Identity(size(), nrow) * qr.matrixQR().triangularView<Eigen::Upper>();  // get R1
  Gt.noalias() = qr.householderQ() * Mat2D::Identity(nrow, size());
  {
    Eigen::HouseholderQR<Eigen::Ref<Mat2D>> qr(Gt);
    Rt.noalias() = Mat2D::Identity(size(), nrow) * qr.matrixQR().triangularView<Eigen::Upper>();  // get R2
    Gt.noalias() = qr.householderQ() * Mat2D::Identity(nrow, size());  // hold Q2 in Gt
  }
  R = Rt * R;  // get R = Rt * R
  // B is size x ncol
  // R.T * B = H.T
  Mat2D B = R.transpose().fullPivHouseholderQr().solve(H.transpose());
  Eigen::JacobiSVD<Mat2D> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // B is K x nsamples, thus we use matrixV() for PCs
  return svd.matrixV().leftCols(nk);
}

void RsvdOpData::computeUSV(int p, double tol) {
  const Index nk{ranks()};
  const Index nrow{rows()};  // nsnps
  const Index ncol{cols()};  // nsamples
  Mat2D Upre, H(ncol, size()), G(nrow, size()), B(size(), ncol), R(size(), size()), Rt(size(), size());
  double diff;
  for (int pi = 0; pi <= p; ++pi) {
    computeGandH(G, H, pi);
    // check if converged
    {
      Eigen::HouseholderQR<Eigen::Ref<Mat2D>> qr(G);
      R.noalias() = Mat2D::Identity(size(), nrow) * qr.matrixQR().triangularView<Eigen::Upper>();  // get R1
      G.noalias() = qr.householderQ() * Mat2D::Identity(nrow, size());  // hold Q1 in G
    }
    {
      Eigen::HouseholderQR<Eigen::Ref<Mat2D>> qr(G);
      Rt.noalias() = Mat2D::Identity(size(), nrow) * qr.matrixQR().triangularView<Eigen::Upper>();  // get R2
      G.noalias() = qr.householderQ() * Mat2D::Identity(nrow, size());  // hold Q2 in G
    }
    R = Rt * R;  // get R = R1R2;
    // B is size x ncol
    // R.T * B = H.T
    B.noalias() = R.transpose().fullPivHouseholderQr().solve(H.transpose());
    Eigen::JacobiSVD<Mat2D> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U = svd.matrixV().leftCols(nk);
    if (pi > 0) {
      if (data->params.mev)
        diff = 1 - mev(U, Upre);
      else
        diff = minSSE(U, Upre).sum() / Upre.cols();
      if (data->params.verbose) cao.print(tick.date(), "running of epoch =", pi, ", diff =", diff);
      if (diff < tol || pi == p) {
        if (data->params.svd_t == SvdType::PCAoneAlg2 && std::pow(2, pi) < data->params.bands) {
          cao.print("PCAone winSVD converged but continues running to get S and V.");
          p = std::log2(data->params.bands);
        } else {
          V.noalias() = G * svd.matrixU().leftCols(nk);
          S = svd.singularValues().head(nk);
          if (data->params.verbose) cao.print(tick.date(), "stops at epoch =", pi + 1);
          break;
        }
      } else {
        Upre = U;
      }
    } else {
      Upre = U;
    }
  }
}

void NormalRsvdOpData::computeGandH(Mat2D& G, Mat2D& H, int pi) {
  if (!data->snpmajor) {
    cao.error("only work with snp major input data now.");
  }

  // reset omg to random
  if (pi == 0 && reset) initOmg();

  if (!data->params.out_of_core) {
    if (pi == 0) {
      if (update) {
        data->update_batch_E(U, S, V.transpose());
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
    G.noalias() = data->G.transpose() * Omg;
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
    if (update) {
      data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
    } else {
      data->read_block_initial(start_idx, stop_idx, standardize);
    }
    data->readtime += tick.reltime();
    G.middleRows(start_idx, actual_block_size).noalias() = data->G.transpose() * Omg;
    H.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
  }
}

void FancyRsvdOpData::computeGandH(Mat2D& G, Mat2D& H, int pi) {
  // check size of G and H first;
  if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
    cao.error("the size of G or H doesn't match with each other.");
  }
  if (pi == 0 && reset) initOmg();
  if (std::pow(2, pi) >= data->params.bands) {
    // reset H1, H2 to zero
    H1.setZero();
    H2.setZero();
  }
  if (!data->params.out_of_core) {
    if (pi == 0) {
      if (update) {
        data->update_batch_E(U, S, V.transpose());
      }
      if (standardize) {
        if (data->params.pcangsd) {
          data->pcangsd_standardize_E(U, S, V.transpose());
        } else {
          data->standardize_E();
        }
      }
      bandsize = 1;
      // blocksize: how many snps in each block
      blocksize = (unsigned int)ceil((double)data->nsnps / data->params.bands);
      if (blocksize < data->params.bands)
        cao.warn("block size < window size. please consider the IRAM method with --svd 0");
      if (data->params.perm) {
        cao.print(tick.date(), "permuting data matrix by columns in place");
        PCAone::permute_matrix(data->G, data->perm);
      }
    }
    // bandsize: how many blocks in each band, 2, 4, 8, 16, 32, 64, ...
    bandsize = fmin(bandsize * 2, data->params.bands);
    // b: the index of current block
    for (uint b = 0, i = 1; b < data->params.bands; ++b, ++i) {
      start_idx = b * blocksize;
      stop_idx = (b + 1) * blocksize >= data->nsnps ? data->nsnps - 1 : (b + 1) * blocksize - 1;
      actual_block_size = stop_idx - start_idx + 1;
      G.middleRows(start_idx, actual_block_size).noalias() =
          data->G.middleCols(start_idx, actual_block_size).transpose() * Omg;

      if (i <= bandsize / 2) {
        // continues to add in data based on current band
        H1.noalias() +=
            data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
      } else {
        H2.noalias() +=
            data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
      }

      // use the first quarter band of succesive iteration (H1)
      // for extra power iteration updates with the last used band (H2)
      const bool adjacent =
          (pi > 0 && (b + 1) == std::pow(2, pi - 1) && std::pow(2, pi) < data->params.bands);
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
  if (pi == 0) bandsize = data->bandFactor;
  data->check_file_offset_first_var();
  // band : 2, 4, 8, 16, 32, 64
  bandsize = fmin(bandsize * 2, data->nblocks);
  for (uint b = 0, i = 1; b < data->nblocks; ++b, ++i) {
    start_idx = data->start[b];
    stop_idx = data->stop[b];
    actual_block_size = stop_idx - start_idx + 1;
    tick.clock();
    if (update) {
      data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
    } else {
      data->read_block_initial(start_idx, stop_idx, standardize);
    }
    data->readtime += tick.reltime();
    G.middleRows(start_idx, actual_block_size).noalias() = data->G.transpose() * Omg;

    if (i <= bandsize / 2) {
      H1.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
    } else {
      H2.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
    }

    const bool adjacent =
        (pi > 0 && (b + 1) == std::pow(2, pi - 1) * data->bandFactor && std::pow(2, pi) < data->params.bands);
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
  Mat1D S;
  RsvdOpData* rsvd;
  if (params.svd_t == SvdType::PCAoneAlg2) {
    cao.print(tick.date(), "initialize window-based RSVD (winSVD) with", params.out_of_core ? "out-of-core" : "in-core");
    rsvd = new FancyRsvdOpData(data, params.k, params.oversamples);
  } else {
    cao.print(tick.date(), "initialize single-pass RSVD (sSVD) with", params.out_of_core ? "out-of-core" : "in-core");
    rsvd = new NormalRsvdOpData(data, params.k, params.oversamples);
  }
  if (!params.impute) {
    if (params.file_t == FileType::PLINK || params.file_t == FileType::BGEN) {
      if (params.ld)
        rsvd->setFlags(false, false);
      else
        rsvd->setFlags(false, true);
    } else {
      rsvd->setFlags(false, false);
    }
    rsvd->computeUSV(params.maxp, params.tol);
  } else {
    // for EM iteration
    rsvd->setFlags(false, false, !params.fancyem);
    rsvd->computeUSV(params.maxp, params.tol);
    // flip_UV(rsvd->U, rsvd->V, false);
    double diff;
    cao.print(tick.date(), "run EM-PCA for data with uncertainty. maxiter =", params.maxiter);
    for (uint i = 0; i < params.maxiter; ++i) {
      rsvd->setFlags(true, false, !params.fancyem);
      Vpre = rsvd->V;
      rsvd->computeUSV(params.maxp, params.tol);
      // flip_UV(rsvd->U, rsvd->V, false);
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

    if (params.pcangsd) {
      cao.print(tick.date(), "estimate GRM for pcangsd");
      data->pcangsd_standardize_E(rsvd->U, rsvd->S, rsvd->V.transpose());
      // TODO: use matrix-free method e.g Arnoldi to decompose the cov
      Mat2D C = data->G * data->G.transpose();
      C.array() /= (double)data->nsnps;
      C.diagonal() = data->Dc.array() / (double)data->nsnps;
      std::ofstream fcov(params.fileout + ".cov");
      fcov << C << "\n";
      // Eigen::SelfAdjointEigenSolver<MyMatrix> eig(C);
      // use Eigen::JacobiSVD to get eigenvecs
      Eigen::JacobiSVD<Mat2D> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
      // output real eigenvectors of covariance in eigvecs2
      write_eigvecs2_beagle(svd.matrixU(), params.filein, params.fileout + ".eigvecs2");
    } else {
      cao.print(tick.date(), "standardize the final matrix for EMU");
      rsvd->setFlags(true, true, !params.fancyem);
      rsvd->computeUSV(params.maxp, params.tol);
    }
  }
  // output PI
  if (params.perm)
    data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->S, rsvd->U, data->perm * rsvd->V);
  else
    data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->S, rsvd->U, rsvd->V);
  if (params.ld) data->write_residuals(rsvd->S, rsvd->U, rsvd->V.transpose());

  delete rsvd;

  cao.print(tick.date(), "PCAone - Randomized SVD done");
  return;
}
