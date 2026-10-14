/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Halko.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Halko.hpp"

#include "Common.hpp"
#include "LinAlg.hpp"
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
    TallQR qr(G);
    R = qr.matrixR();
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
          V.resize(nrow, nk);
          qr.applyQ(svd.matrixU().leftCols(nk), V);  // Q * W
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

void RsvdOpData::add_block(uint64 start_idx, uint64 stop_idx, Mat2D& G, Mat2D& Hb) {
  const uint64 n = stop_idx - start_idx + 1;
  tick.clock();
  if (data->packed_pca()) {
    // the 2-bit calls of the block straight into the products (PackedPca.hpp)
    const uint8_t* rec = data->read_block_packed(start_idx, stop_idx);
    data->readtime += tick.reltime();
    data->set_packed_sites(standardize, start_idx, n);
    data->packed->mul_Xt_Y(start_idx, Omg, G.middleRows(start_idx, n), rec);
    data->packed->mul_X_Y(start_idx, G.middleRows(start_idx, n), Hb, true, rec);
    return;
  }
  if (winem) {
    data->read_block_initial(start_idx, stop_idx, false);  // missing calls are NaN
    data->readtime += tick.reltime();
    estep_block(start_idx, data->G.leftCols(n));
    mul_Xt_Y(data->G, Omg, G.middleRows(start_idx, n));
    mul_X_Y(data->G, G.middleRows(start_idx, n), Hb, true);
    return;
  }
  if (!update) {
    data->read_block_initial(start_idx, stop_idx, standardize);
  } else {
    data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
  }
  data->readtime += tick.reltime();
  mul_Xt_Y(data->G, Omg, G.middleRows(start_idx, n));
  mul_X_Y(data->G, G.middleRows(start_idx, n), Hb, true);
}

void RsvdOpData::update_em_basis(const TallQR& qr) {
  Eigen::JacobiSVD<Mat2D> svd(qr.matrixR(), Eigen::ComputeThinU);
  Uem.resize(cols(), ranks());
  qr.applyQ(svd.matrixU().leftCols(ranks()), Uem);
}

void RsvdOpData::estep_block(uint64 start_idx, Eigen::Ref<Mat2D> Z) {
  data->estep_reg(Uem, Z, start_idx);
  if (!standardize || data->params.scale != SCALE_STANDARDIZE_GENETIC) return;
#pragma omp parallel for
  for (Eigen::Index i = 0; i < Z.cols(); ++i) {
    const double f = data->F(start_idx + i), sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) Z.col(i) *= sqrt((double)data->params.ploidy) / sd;
  }
}

void RsvdOpData::runWinEM(uint maxpass, double tol) {
  winem = 1;
  setFlags(false, false);
  pi_offset = 0;
  Uem.resize(0, 0);
  // the doubling phase of the band, then two passes of the full band
  const uint minpass = (uint)std::ceil(std::log2((double)data->params.bands)) + 3;
  Mat2D G(rows(), size()), H(cols(), size()), Upre;
  double diff = 1.0;
  uint pi = 0;
  for (; pi < maxpass; ++pi) {
    computeGandH(G, H, pi);
    if (pi > 0) diff = 1.0 - mev(Uem, Upre);
    char d[32];  // the log may be in std::fixed by now, which shows 1e-7 as 0.000000
    snprintf(d, sizeof d, "%.3g", diff);
    cao.print(tick.date(), "winEM pass =", pi + 1, ", diff =", d);
    if (pi + 1 >= minpass && diff < tol) break;
    Upre = Uem;
  }
  if (pi < maxpass)
    cao.print(tick.date(), "winEM came to convergence after", pi + 1, "passes; E-step", data->esteptime, "seconds");
  else
    warn_em_not_converged(data->params, diff);
  U = Uem;
  winem = 2;
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
      if (data->packed) {
        data->set_packed_sites(standardize);  // no EM: packed_pca()
      } else {
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
    }
    if (pi > 0) {
      TallQR qr(H);
      qr.applyQ(Mat2D::Identity(size, size), Omg);  // hold H in Omega
      PCAone::flipOmg(Omg2, Omg);
    }
    if (data->packed) {
      data->packed->mul_Xt_Y(0, Omg, G);
      data->packed->mul_X_Y(0, G, H);
    } else {
      mul_Xt_Y(data->G, Omg, G);
      mul_X_Y(data->G, G, H);
    }
    return;
  }

  // for block version
  // data->G is always nsamples x nsnps;
  // for nsnps > nsamples
  if (pi > 0) {
    TallQR qr(H);
    qr.applyQ(Mat2D::Identity(size, size), Omg);
  }
  H = Mat2D::Zero(cols(), size);
  data->check_file_offset_first_var();
  for (uint i = 0; i < data->nblocks; ++i) add_block(data->start[i], data->stop[i], G, H);
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
    if (pi == 0 && !data->packed) {
      if (winem == 2) {
        data->estep(Uem, data->G, 0);
      } else if (update) {
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
    if (pi == 0) {
      bandsize = pi_offset > 0 ? data->params.bands : 1;
      // blocksize: how many snps in each block
      blocksize = (unsigned int)ceil((double)data->nsnps / data->params.bands);
      if (blocksize < data->params.bands)
        cao.warn("block size < window size. please consider the IRAM method with --svd 0");
      // Keep the same column order across EM iterations and final scaling.
      // (the 2-bit calls were read in the shuffled order: FileBed::read_all_packed())
      if (data->params.perm && !data->in_core_permuted) {
        cao.print(tick.date(), "permuting data matrix by columns in place");
        PCAone::permute_matrix(data->G, data->perm, data->params.seed);
        data->in_core_permuted = true;
      }
      if (data->packed) data->set_packed_sites(standardize);
    }
    // bandsize: how many blocks in each band, 2, 4, 8, 16, 32, 64, ...
    bandsize = fmin(bandsize * 2, data->params.bands);
    // The schedule of the band: each block b adds into H1 or H2, and Omg is
    // updated after some of them. A run of blocks with the same Omg that add
    // into the same H is one product over their sites, which are contiguous:
    // 2 products an epoch once the band is full, in place of one per block.
    struct Step {
      bool h1, update, full;
    };
    std::vector<Step> plan(data->params.bands);
    for (uint b = 0, i = 1; b < data->params.bands; ++b, ++i) {
      // continues to add in data based on current band
      plan[b].h1 = i <= bandsize / 2;
      // use the first quarter band of succesive iteration (H1)
      // for extra power iteration updates with the last used band (H2)
      const bool adjacent = (pe > 0 && (b + 1) == std::pow(2, pe - 1) && std::pow(2, pe) < data->params.bands);
      // add up H and update Omg
      plan[b].update = !((b + 1) < bandsize && !adjacent) && ((i == bandsize) || (i == bandsize / 2) || adjacent);
      plan[b].full = plan[b].update && i == bandsize;
      if (plan[b].full) i = 0;
    }
    for (uint b0 = 0, b1; b0 < data->params.bands; b0 = b1) {
      b1 = b0 + 1;
      while (b1 < data->params.bands && !plan[b1 - 1].update && plan[b1].h1 == plan[b0].h1) ++b1;
      // half-open [start_idx, stop_idx). With fewer than bands^2 sites the
      // trailing blocks are empty: ceil(M / bands) * b can pass M, and the old
      // inclusive stop_idx then underflowed the unsigned block size, so the
      // next line asked for ~2^64 rows (std::bad_alloc). An empty block
      // contributes nothing but keeps the band schedule intact.
      start_idx = std::min<uint64>((uint64)b0 * blocksize, data->nsnps);
      stop_idx = std::min<uint64>((uint64)b1 * blocksize, data->nsnps);
      actual_block_size = stop_idx - start_idx;
      Mat2D& Hb = plan[b0].h1 ? H1 : H2;
      // winEM: the E-step of the run, from the PCs of the last update
      if (winem == 1 && Uem.cols() > 0)
        data->estep(Uem, data->G.middleCols(start_idx, actual_block_size), start_idx);
      if (data->packed) {
        data->packed->mul_Xt_Y(start_idx, Omg, G.middleRows(start_idx, actual_block_size));
        data->packed->mul_X_Y(start_idx, G.middleRows(start_idx, actual_block_size), Hb, true);
      } else {
        mul_Xt_Y(data->G.middleCols(start_idx, actual_block_size), Omg, G.middleRows(start_idx, actual_block_size));
        mul_X_Y(data->G.middleCols(start_idx, actual_block_size), G.middleRows(start_idx, actual_block_size), Hb, true);
      }
      const Step& last = plan[b1 - 1];
      if (!last.update) continue;
      H = H1 + H2;
      Mat2D Hq = H;  // H is read again by computeUSV(): factor a copy
      TallQR qr(Hq);
      qr.applyQ(Mat2D::Identity(size, size), Omg);
      PCAone::flipOmg(Omg2, Omg);
      if (winem == 1) update_em_basis(qr);
      if (last.full)
        H1.setZero();
      else
        H2.setZero();
    }

    return;
  }

  // out-of-core implementation
  if (pi == 0) bandsize = pi_offset > 0 ? data->nblocks : data->bandFactor;
  data->check_file_offset_first_var();
  // band : 2, 4, 8, 16, 32, 64
  bandsize = fmin(bandsize * 2, data->nblocks);
  for (uint b = 0, i = 1; b < data->nblocks; ++b, ++i) {
    add_block(data->start[b], data->stop[b], G, i <= bandsize / 2 ? H1 : H2);

    const bool adjacent =
        (pe > 0 && (b + 1) == std::pow(2, pe - 1) * data->bandFactor && std::pow(2, pe) < data->params.bands);
    if ((b + 1) < bandsize && !adjacent) continue;

    // cao.print("i:", i, ",j:", j, ",bandsize:", bandsize, ",pi:", pi);
    // add up H and update Omg
    if (!((i == bandsize) || (i == bandsize / 2) || adjacent)) continue;
    H = H1 + H2;
    Mat2D Hq = H;  // H is read again by computeUSV(): factor a copy
    TallQR qr(Hq);
    qr.applyQ(Mat2D::Identity(size, size), Omg);
    PCAone::flipOmg(Omg2, Omg);
    if (winem == 1) update_em_basis(qr);
    if (i == bandsize) {
      H1.setZero();
      i = 0;
    } else {
      H2.setZero();
    }
  }
}

static RsvdOpData* new_rsvd(Data* data, const Param& params, uint k, uint os) {
  if (params.svd_t == SvdType::PCAoneAlg2) return new FancyRsvdOpData(data, k, os);
  return new NormalRsvdOpData(data, k, os);
}

void run_pca_with_halko(Data* data, const Param& params) {
  Mat2D Vpre;
  cao.print(tick.date(), "initialize",
            params.svd_t == SvdType::PCAoneAlg2 ? "window-based RSVD (winSVD) with" : "single-pass RSVD (sSVD) with",
            params.out_of_core ? "out-of-core" : "in-core", data->packed_pca() ? "on the 2-bit calls" : "");
  // EM-PCA fits the individual allele frequencies with em_k PCs, which is -k
  // unless --em-k says otherwise; the final matrix is decomposed for -k
  RsvdOpData* rsvd = params.missme ? new_rsvd(data, params, params.em_k, params.em_oversamples)
                                   : new_rsvd(data, params, params.k, params.oversamples);
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
  }
  if (params.missme && params.winem) {
    // winEM: the passes give the em_k PCs, and the final matrix is imputed from
    // them, standardized and decomposed for -k PCs from there, as below
    cao.print(tick.date(), "run winEM. maximum passes =", params.maxiter);
    if (params.out_of_core) {
      // the blocks hold the missing calls as NaN (BED only: Cmd.cpp), also if
      // the lookup of the centred calls was built before, with them at 0
      data->mark_missing = true;
      if (data->centered_geno_lookup.cols() > 0)
        data->centered_geno_lookup.row(1).setConstant(std::numeric_limits<double>::quiet_NaN());
    }
    rsvd->runWinEM(params.maxiter, params.tolem);
    if (params.em_k != params.k) {
      cao.print(tick.date(), "decompose the final matrix for", params.k, "PCs (-k), imputed from the", params.em_k,
                "PCs of winEM (--em-k)");
      RsvdOpData* out = new_rsvd(data, params, params.k, params.oversamples);
      out->U = rsvd->U;
      out->Uem = std::move(rsvd->Uem);
      out->winem = 2;
      delete rsvd;
      rsvd = out;
    }
    if (params.emu) cao.print(tick.date(), "standardize the final matrix for EMU");
    // update: warm-start from U. EMU decomposes its final matrix standardized,
    // PCAngsd its expected genotypes
    rsvd->setFlags(true, params.emu);
    rsvd->computeUSV(params.maxp, params.tol);
    flip_UV(rsvd->U, rsvd->V, false);
    if (params.pcangsd && params.file_t == FileType::BEAGLE) {
      // the individual allele frequencies of the em_k PCs: U U'E, with E the
      // expected genotypes the final decomposition imputed from them
      cao.print(tick.date(), "estimate GRM for pcangsd");
      const Mat2D& Uem = rsvd->Uem;
      const Mat2D Wt = Uem.transpose() * data->G;  // em_k x nsnps
      data->pcangsd_standardize_E(Uem, Mat1D::Ones(Uem.cols()), Wt);
      data->write_pcangsd_cov();
    }
  } else if (params.missme) {
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

    // PCAngsd's covariance, and the -k eigenvectors of it in .eigvecs2, come from
    // the individual allele frequencies of the --em-k PCs the EM fitted
    if (params.pcangsd && (params.file_t == FileType::BEAGLE)) {
      cao.print(tick.date(), "estimate GRM for pcangsd");
      data->pcangsd_standardize_E(rsvd->U, rsvd->S, rsvd->V.transpose());
      // TODO: use matrix-free method e.g Arnoldi to decompose the cov
      data->write_pcangsd_cov();
    }

    // EMU decomposes its final matrix standardized; PCAngsd writes the PCs of
    // the EM's matrix of expected genotypes, which needs one more decomposition
    // only when --em-k fitted another number of PCs than the -k written.
    if (params.emu || params.em_k != params.k) {
      if (params.emu) cao.print(tick.date(), "standardize the final matrix for EMU");
      if (params.em_k != params.k) {
        // The update below rebuilds the matrix from the U, S, V it is given, so
        // a solver for -k PCs takes over the fit of the --em-k ones: they give
        // the individual allele frequencies, and start its power iterations.
        cao.print(tick.date(), "decompose the final matrix for", params.k,
                  "PCs (-k), given the individual allele frequencies of the", params.em_k, "PCs of the EM (--em-k)");
        RsvdOpData* out = new_rsvd(data, params, params.k, params.oversamples);
        out->U = std::move(rsvd->U);
        out->S = std::move(rsvd->S);
        out->V = std::move(rsvd->V);
        delete rsvd;
        rsvd = out;
      }
      rsvd->setFlags(true, params.emu);
      rsvd->computeUSV(params.maxp, params.tol);
      flip_UV(rsvd->U, rsvd->V, false);
    }
  }
  // output PI
  data->set_svd_transform(rsvd->standardize);
  data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->S, rsvd->U, rsvd->V);

  delete rsvd;

  cao.print(tick.date(), "PCAone - Randomized SVD done");
  return;
}
