/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Projection.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Projection.hpp"

#include "Cmd.hpp"
#include "Common.hpp"
#include "Data.hpp"
#include "Utils.hpp"

static void solve_projection_scores(const Mat2D& V, const ArrBool& C, const Mat2D& G, Mat2D& U) {
  if (U.rows() == 0 || U.cols() == 0) return;
  double p_miss = C.size() ? (double)C.count() / (double)C.size() : 0.0;
  if (p_miss == 0.0) {
    Eigen::ColPivHouseholderQR<Mat2D> qr(V);
#pragma omp parallel for
    for (uint i = 0; i < (uint)U.rows(); i++) {
      U.row(i) = qr.solve(G.row(i).transpose());
    }
  } else {
#pragma omp parallel for
    for (uint i = 0; i < (uint)U.rows(); i++) {
      Int1D idx;
      for (int j = 0; j < V.rows(); j++) {
        if (!C(j * U.rows() + i)) idx.push_back(j);
      }
      U.row(i) = V(idx, Eigen::all).colPivHouseholderQr().solve(G(i, idx).transpose());
    }
  }
}

/**
 * options:
 * 1: simple, assume no missingness
 * 2: like smartPCA, solving g=Vx, can take missing genotypes
 * 3: iterative GL-aware projection (EM): alternates between updating individual allele frequencies
 *    from current PC scores and re-solving for PC scores with updated expected genotypes (BEAGLE only)
   // NOTE: we don't support out-of-core for projection.
 */
void run_projection(Data* data, const Param& params) {
  BimMatch match;
  if (params.file_t == FileType::BEAGLE) {
    match = match_beagle_to_mbim(params.filein, params.filebim);
    if (match.bim_indices.empty())
      cao.error("no overlapped SNPs found between " + params.filein + " and " + params.filebim);
  } else {
    match = match_bim_to_mbim(params.filein + ".bim", params.filebim);
    if (match.bim_indices.empty())
      cao.error("no overlapped SNPs found between " + params.filein + ".bim and " + params.filebim);
  }
  if (!match.identical) {
    data->keepSNPs = match.bim_indices;
    data->keepRefSNPs = match.mbim_indices;
    for (int k = 0; k < (int)match.flip.size(); ++k) {
      if (match.flip[k]) data->flipSNPs.push_back(k);
    }
    cao.warn("SNP info is not fully identical between input and reference .mbim");
    cao.print(tick.date(), "projection will use", match.bim_indices.size(), " overlapped sites. there are",
              data->flipSNPs.size(), " flipped alleles");
    if (!data->flipSNPs.empty())
      cao.warn(data->flipSNPs.size(), " SNPs have flipped ref/alt alleles and will be corrected");
  }
  cao.print(tick.date(), "run projection");
  data->prepare();  // read AF and resize F to matched size
  if (params.project != 3) data->standardize_E();
  cao.print(tick.date(), "start parsing V:", params.fileV, ", S:", params.fileS);
  uint nsamples, nsnps;
  Mat1D S;
  read_sigvals(params.fileS, nsamples, nsnps, S);
  // target number of PCs for getting individual allele frequency
  const int K = fmin(S.size(), params.k);
  Mat2D V = read_eigvecs(params.fileV, nsnps, K);
  if (!match.identical) {
    Mat2D V_overlap(match.mbim_indices.size(), K);
    for (int i = 0; i < (int)match.mbim_indices.size(); ++i) {
      V_overlap.row(i) = V.row(match.mbim_indices[i]);
      if (match.flip[i]) V_overlap.row(i) = -V_overlap.row(i);
    }
    V = V_overlap;
  }
  Mat2D U(data->nsamples, K);
  double p_miss = data->C.size() ? (double)data->C.count() / (double)data->C.size() : data->p_miss;

  if (params.project == 1) {
    if (p_miss > 0) cao.warn("there are missing genotypes. recommend using --project 2 or 3.");
    // get 1 / Singular = sqrt(Eigen * M)
    V = V * (S.array().inverse().matrix().asDiagonal());
    // G V = U D
    U = data->G * V;
  } else if (params.project == 2) {
    V = V * S.asDiagonal();
    if (p_miss == 0.0) cao.warn("there is no missing genotypes");
    solve_projection_scores(V, data->C, data->G, U);
  } else if (params.project == 3) {
    // project == 3: iterative GL-aware projection (EM)
    // E-step: update expected genotype G using per-sample allele frequencies from U
    // M-step: solve V*S*u_i = G_std_i (standardized with reference F)
    if (params.file_t != FileType::BEAGLE) cao.error("--project 3 requires BEAGLE genotype likelihood input");
    const bool filter = !data->keepSNPs.empty();

    V = V * S.asDiagonal();
    solve_projection_scores(V, data->C, data->G, U);

    // run EM
    for (uint iter = 0; iter < params.maxiter; ++iter) {
      Mat2D Uprev = U;

      // E-step: update G using individual allele frequencies
      // pred.noalias() = U * sdiag * V.transpose();  // nsamples x nsnps
#pragma omp parallel for
      for (uint j = 0; j < data->nsnps; ++j) {
        const double norm = sqrt(2.0 * data->F(j) * (1.0 - data->F(j)));
        uint s = filter ? data->keepSNPs[j] : (uint)j;
        for (uint i = 0; i < data->nsamples; ++i) {
          double pt = 0.0;
          for (int k = 0; k < K; ++k) {
            pt += U(i, k) * V(j, k);
          }
          if (params.scale == -9 && norm > VAR_TOL) pt *= norm;
          pt = fmin(fmax(pt + data->F(j), 1e-4), 1.0 - 1e-4);
          const double p0 = data->P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
          const double p1 = data->P(2 * i + 1, s) * 2.0 * pt * (1.0 - pt);
          const double p2 = (1.0 - data->P(2 * i + 0, s) - data->P(2 * i + 1, s)) * pt * pt;
          const double psum = p0 + p1 + p2;
          if (!std::isfinite(psum) || psum <= 0.0) {
            data->C[j * data->nsamples + i] = 1;
            data->G(i, j) = 0.0;
            continue;
          }
          data->C[j * data->nsamples + i] = 0;
          data->G(i, j) = (p1 + 2.0 * p2) / (2.0 * psum) - data->F(j);  // domain (0,1)
          if (params.scale == -9) {
            if (norm > VAR_TOL)
              data->G(i, j) /= norm;
            else
              data->G(i, j) = 0.0;
          }
        }
      }

      // M-step: standardize G and solve for U
      solve_projection_scores(V, data->C, data->G, U);
      // flip signs
      for (int k = 0; k < K; ++k) {
        if (U.col(k).dot(Uprev.col(k)) < 0.0) U.col(k) *= -1.0;
      }
      double denom = Uprev.norm();
      if (denom < 1e-12) denom = 1.0;
      double diff = (U - Uprev).norm() / denom;
      ;
      cao.print(tick.date(), "GL projection iter", iter + 1, ", diff =", diff);
      if (diff < params.tolem) break;
    }
  } else {
    cao.error("unsupported --project mode: " + std::to_string(params.project));
  }

  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  std::ofstream outu(params.fileout + ".eigvecs");
  if (outu.is_open()) outu << U.format(fmt) << '\n';
}
