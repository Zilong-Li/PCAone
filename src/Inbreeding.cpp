/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Inbreeding.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "Inbreeding.hpp"

#include "FilePlink.hpp"

// PI is N x M
// GL is (N x 2) x M
// Fd is M x 2. The 1st column is new F and the 2nd is the F_new - F
void calc_inbreed_coef(Mat2D& Fd, const Mat1D& F, const Mat2D& PI, const Mat2D& GL) {
  const int nsnps = PI.cols();
  const int nsamples = PI.rows();
#pragma omp parallel for
  for (int j = 0; j < nsnps; j++) {
    double obsH = 0.0, expH = 0.0;
    for (int i = 0; i < nsamples; i++) {
      // Fadj = (1-pi)*pi*F;
      double Fadj = (1.0 - PI(i, j)) * PI(i, j) * F(j);
      // (1-pi)^2 + pi(1-pi)*F
      double p0 = fmax(1e-4, (1 - PI(i, j) * (1 - PI(i, j)) + Fadj));
      // 2pi(1-pi)(1-F)
      double p1 = fmax(1e-4, 2 * PI(i, j) * (1 - PI(i, j)) - 2.0 * Fadj);
      // pi^2 + pi(1-pi)*F
      double p2 = fmax(1e-4, PI(i, j) * PI(i, j) + Fadj);
      // normalize
      double pSum = p0 + p1 + p2;
      p0 /= pSum;
      p1 /= pSum;
      p2 /= pSum;

      // posterior
      double pp0 = GL(2 * i + 0, j) * p0;
      double pp1 = GL(2 * i + 1, j) * p1;
      double pp2 = (1.0 - GL(2 * i + 0, j) - GL(2 * i + 1, j)) * p2;
      double ppSum = pp0 + pp1 + pp2;

      // sum over indiiduals
      obsH += pp1 / ppSum;
      // count heterozygotes
      expH += 2.0 * PI(i, j) * (1 - PI(i, j));
    }
    // ANGSD procedure
    obsH = fmax(1e-4, obsH / (double)nsamples);
    // Update the inbreeding coefficient
    Fd(j, 0) = 1.0 - (nsamples * obsH / expH);
    Fd(j, 0) = fmin(fmax(-1.0, Fd(j, 0)), 1.0);
    Fd(j, 1) = Fd(j, 0) - F(j);
  }
}

void run_inbreeding_em(Mat1D& F, const Mat2D& PI, const Mat2D& GL, const Param& params) {
  const int nsnps = PI.cols();
  Mat2D Fd(nsnps, 2);
  Mat1D f1(nsnps), d1(nsnps);
  double sr2, sv2, alpha, diff;
  AreClose areClose;
  for (uint it = 0; it < params.maxiter; it++) {
    calc_inbreed_coef(Fd, F, PI, GL);
    f1 = Fd.col(0);
    d1 = Fd.col(1);
    sr2 = d1.array().square().sum();
    calc_inbreed_coef(Fd, f1, PI, GL);
    f1 = F;  // f1 is F0
    sv2 = (Fd.col(1) - d1).array().square().sum();
    // safety break
    if (areClose(sv2, 0.0)) {
      F = Fd.col(0);  // copy
      cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, "RMSE = 0.0");
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    alpha = -fmax(1.0, sqrt(sr2 / sv2));
    F = F + 2 * alpha * d1 + alpha * alpha * (Fd.col(1) - d1);
    // map to domain [-1, 1]
    F = (F.array() < -1.0).select(-1.0, F);
    F = (F.array() > 1.0).select(1.0, F);
    // Stabilization step and convergence check
    calc_inbreed_coef(Fd, F, PI, GL);
    diff = rmse1d(f1, Fd.col(0));
    cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, "RMSE =", diff);
    if (diff < params.tolem) {
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    if (it == params.maxiter - 1) cao.warn("EM inbreeding coefficient not coverged!");
  }
}

void run_inbreeding(Data* data, const Param& params) {
  cao.print(tick.date(), "run inbreeding coefficient estimator");
  data->prepare();
  if (params.file_t == FileType::BEAGLE) {
    data->P = Mat2D::Zero(data->nsamples * 2, data->nsnps);  // genotype likelihood
    gzFile fp = gzopen(params.filein.c_str(), "r");
    parse_beagle_file(data->P, fp, data->nsamples, data->nsnps);
    gzclose(fp);
    // run EM-HWE
    data->F = Mat1D::Zero(data->nsnps);  // init inbreeding coef
    run_inbreeding_em(data->F, data->G, data->P, params);
  }
}
