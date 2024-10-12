/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/LD.cpp
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

  if (params.project == 1) {
    if (p_miss > 0)
      cao.warn(
          "there are missing genotypes. recommend the alternative method "
          "--project 2 or 3.");
    // get 1 / Singular = sqrt(Eigen * M)
    Mat1D S = 1.0 / (read_eigvals(params.fileS) * data->nsnps / params.ploidy).array().sqrt();
    int pcs = S.size();
    // G V = U D
    Mat2D V = read_eigvecs(params.fileV, data->nsnps, pcs) * S.asDiagonal();  // M x K
    Mat2D U = data->G * V;

    data->write_eigs_files(1.0 / S.array(), U, V);
  } else if (params.project == 2) {
    // get  Singular = sqrt(Eigen * M)
    Mat1D S = (read_eigvals(params.fileS) * data->nsnps / params.ploidy).array().sqrt();
    int pcs = S.size();
    Mat2D V = read_eigvecs(params.fileV, data->nsnps, pcs) * S.asDiagonal();  // M x K
    Mat2D U(data->nsamples, pcs);

    if (p_miss == 0.0) {
      cao.warn("there is no missing genotypes");
      Eigen::BDCSVD<Mat2D> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);

      for (uint i = 0; i < data->nsamples; i++) {
        // Vx = g
        U.row(i) = svd.solve(data->G.row(i).transpose());
      }
    } else {
      Int1D idx;
      for (uint i = 0; i < data->nsamples; i++) {
        // find non-missing snps for sample i
        idx.clear();
        for (int j = 0; j < V.rows(); j++) {
          if (!data->C(j * data->nsamples + i)) idx.push_back(j);
        }
        Eigen::BDCSVD<Mat2D> svd(V(idx, Eigen::all), Eigen::ComputeThinU | Eigen::ComputeThinV);
        Mat1D g = data->G(i, idx);
        // Vx = g
        U.row(i) = svd.solve(g);
      }
    }

    data->write_eigs_files(S, U, V);
  } else {
    cao.error("have not implemented yet");
  }
}
