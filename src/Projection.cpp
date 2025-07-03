/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Projection.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Projection.hpp"

#include "Utils.hpp"

/**
 * options:
 * 1: simple, assume no missingness
 * 2: like smartPCA, solving g=Vx, can take missing genotypes
 * 3: OADP, laser2, can take missing genotypes
 */
void run_projection(Data* data, const Param& params) {
  check_bim_vs_mbim(params.filein + ".bim", params.filebim);
  cao.print(tick.date(), "run projection");

  data->prepare();
  data->standardize_E();
  double p_miss = (double)data->C.count() / (double)data->C.size();
  // Singular = sqrt(Eigen * M)
  Mat1D S = (read_eigvals(params.fileS) * data->nsnps / params.ploidy).array().sqrt();
  const int pcs = S.size();
  Mat2D V = read_eigvecs(params.fileV, data->nsnps, pcs);  // M x K
  Mat2D U(data->nsamples, pcs);

  if (params.project == 1) {
    if (p_miss > 0) cao.warn("there are missing genotypes. recommend using --project 2 or 3.");
    // get 1 / Singular = sqrt(Eigen * M)
    V = V * (S.array().inverse().matrix().asDiagonal());
    // G V = U D
    U = data->G * V;
  } else if (params.project == 2) {
    V = V * S.asDiagonal();
    if (p_miss == 0.0) {
      cao.warn("there is no missing genotypes");
      Eigen::ColPivHouseholderQR<Mat2D> qr(V);
#pragma omp parallel for
      for (uint i = 0; i < data->nsamples; i++) {
        // Vx = g
        U.row(i) = qr.solve(data->G.row(i).transpose());
      }
    } else {
#pragma omp parallel for
      for (uint i = 0; i < data->nsamples; i++) {
        // find non-missing snps for sample i
        Int1D idx;
        for (int j = 0; j < V.rows(); j++) {
          if (!data->C(j * data->nsamples + i)) idx.push_back(j);
        }
        // Vx = g
        U.row(i) = V(idx, Eigen::all).colPivHouseholderQr().solve(data->G(i, idx).transpose());
      }
    }
  } else {
    cao.error("have not implemented yet");
  }

  data->write_eigs_files(S.array().square() / data->nsnps, S, U, V);
}
