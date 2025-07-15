/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Inbreeding.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "Inbreeding.hpp"

#include "Cmd.hpp"
#include "FileBeagle.hpp"
#include "FilePlink.hpp"
#include "Utils.hpp"

// type 1: Genotype input, {1, -9, 0.5, 0}, GL is N x M
// type 2: Genotype likelihood input, GL is (N x 2) x M
void calc_inbreed_coef(Mat1D& D, Mat1D& F, const Mat2D& PI, const Mat2D& GL, const int type, const int size,
                       const uint start) {
  const int nsnps = size;
  const int nsamples = PI.rows();
  if (type != 1 && type != 2) cao.error("type must be 1 or 2");
#pragma omp parallel for
  for (int j = 0; j < nsnps; j++) {
    double obsH = 0.0, expH = 0.0, Fadj, p0, p1, p2, pSum, pp0, pp1, pp2;
    int jj = j + start;
    for (int i = 0; i < nsamples; i++) {
      // Fadj = (1-pi)*pi*F;
      Fadj = (1.0 - PI(i, j)) * PI(i, j) * F(jj);
      // (1-pi)^2 + pi(1-pi)*F
      p0 = fmax(1e-4, (1.0 - PI(i, j)) * (1.0 - PI(i, j)) + Fadj);
      // 2pi(1-pi)(1-F)
      p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) * (1.0 - F(jj)));
      // pi^2 + pi(1-pi)*F
      p2 = fmax(1e-4, PI(i, j) * PI(i, j) + Fadj);
      // normalize
      pSum = 1.0 / (p0 + p1 + p2);
      p0 *= pSum;
      p1 *= pSum;
      p2 *= pSum;

      // posterior
      if (type == 1) {
        // should check genotype missingness
        if (GL(i, j) == 0.0) {
          pp0 = p0, pp1 = 0, pp2 = 0;
        } else if (GL(i, j) == 0.5) {
          pp0 = 0, pp1 = p1, pp2 = 0;
        } else if (GL(i, j) == 1.0) {
          pp0 = 0, pp1 = 0, pp2 = p2;
        } else {
          cao.error("missing genotypes found");
        }
      } else {
        pp0 = GL(2 * i + 0, j) * p0;
        pp1 = GL(2 * i + 1, j) * p1;
        pp2 = (1.0 - GL(2 * i + 0, j) - GL(2 * i + 1, j)) * p2;
      }
      pSum = pp0 + pp1 + pp2;

      // sum over indiiduals
      obsH += pp1 / pSum;
      // count heterozygotes
      expH += 2.0 * PI(i, j) * (1.0 - PI(i, j));
    }
    // Update the inbreeding coefficient
    double f = fmin(fmax(-1.0, 1.0 - (obsH / expH)), 1.0);
    D(jj) = f - F(jj);  // Fnew - Fold
    F(jj) = f;          // update with Fnew
  }
}

void calc_inbreed_site_lrt(Mat1D& T, const Mat1D& F, const Mat2D& PI, const Mat2D& GL, const int type,
                           const int size, const uint start) {
  if (type != 1 && type != 2) cao.error("type must be 1 or 2");
  const int nsnps = size;
  const int nsamples = PI.rows();
#pragma omp parallel for
  for (int j = 0; j < nsnps; j++) {
    double logAlt = 0.0, logNull = 0.0, Fadj, p0, p1, p2, pSum, l0, l1, l2;
    int jj = j + start;
    for (int i = 0; i < nsamples; i++) {
      // Fadj = (1-pi)*pi*F;
      Fadj = (1.0 - PI(i, j)) * PI(i, j) * F(jj);
      // (1-pi)^2 + pi(1-pi)*F
      p0 = fmax(1e-4, (1.0 - PI(i, j)) * (1.0 - PI(i, j)) + Fadj);
      // 2pi(1-pi)(1-F)
      p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) * (1.0 - F(jj)));
      // pi^2 + pi(1-pi)*F
      p2 = fmax(1e-4, PI(i, j) * PI(i, j) + Fadj);
      // normalize
      pSum = 1.0 / (p0 + p1 + p2);
      p0 *= pSum;
      p1 *= pSum;
      p2 *= pSum;

      // posterior // Likelihood*prior
      if (type == 1) {
        if (GL(i, j) == 0.0) {
          logAlt += log(p0);
          logNull += log((1.0 - PI(i, j)) * (1.0 - PI(i, j)));
        } else if (GL(i, j) == 0.5) {
          logAlt += log(p1);
          logNull += log(2.0 * PI(i, j) * (1.0 - PI(i, j)));
        } else if (GL(i, j) == 1.0) {
          logAlt += log(p2);
          logNull += log(PI(i, j) * PI(i, j));
        } else {
          cao.error("missing genotypes found");
        }
      } else {
        l0 = GL(2 * i + 0, j) * p0;
        l1 = GL(2 * i + 1, j) * p1;
        l2 = (1.0 - GL(2 * i + 0, j) - GL(2 * i + 1, j)) * p2;
        logAlt += log(l0 + l1 + l2);
        // Null model
        l0 = GL(2 * i + 0, j) * (1.0 - PI(i, j)) * (1.0 - PI(i, j));
        l1 = GL(2 * i + 1, j) * 2.0 * PI(i, j) * (1.0 - PI(i, j));
        l2 = (1.0 - GL(2 * i + 0, j) - GL(2 * i + 1, j)) * PI(i, j) * PI(i, j);
        logNull += log(l0 + l1 + l2);
      }
    }
    T(jj) = 2.0 * (logAlt - logNull);
  }
}

void write_hwe_per_site(const std::string& fout, const std::string& fbim, const Mat1D& hwe, const Mat1D& lrt,
                        const Mat1D& coef) {
  std::ifstream fin(fbim);
  if (!fin.is_open()) cao.error("can not open " + fbim);
  std::ofstream ohwe(fout);
  if (!ohwe.is_open()) cao.error("can not open " + fout);
  std::string line;
  const std::string sep{"\t"};
  int j = 0;
  ohwe << "#ID\tHWE_P\tLRT\tInbreeding_coefficient\n";
  while (getline(fin, line)) {
    auto id = split_string(line, sep)[1];
    ohwe << id << "\t" << hwe(j) << "\t" << lrt(j) << "\t" << coef(j) << "\n";
    j++;
  }
  assert(j == coef.size());
}

void run_inbreeding_em(int type, const Mat2D& GL, const Mat2D& PI, const Param& params) {
  if (type != 1 && type != 2) cao.error("type must be 1 or 2");
  const int nsnps = PI.cols();
  Mat1D F = Mat1D::Zero(nsnps), F0(nsnps);  // init inbreeding coef
  Mat1D D1(nsnps), D2(nsnps);               // store the diff between Fnew - Fold
  double sr2, sv2, alpha, diff;
  AreClose areClose;
  calc_inbreed_coef(D1, F, PI, GL, type, nsnps, 0);  // init F
  for (uint it = 0; it < params.maxiter; it++) {
    F0 = F;                                            // copy the initial F
    calc_inbreed_coef(D1, F, PI, GL, type, nsnps, 0);  // F is F1, D1 = F1 - F0
    sr2 = D1.array().square().sum();
    calc_inbreed_coef(D2, F, PI, GL, type, nsnps, 0);  // F is F2, D2 = F2 - F1
    sv2 = (D2 - D1).array().square().sum();
    // safety break
    if (areClose(sv2, 0.0)) {
      cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE = 0.0");
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    alpha = fmin(fmax(1.0, sqrt(sr2 / sv2)), 256.0);
    if (params.verbose > 1) cao.print("alpha:", alpha, ", sr2:", sr2, ", sv2:", sv2);
    F = F0 + 2 * alpha * D1 + alpha * alpha * (D2 - D1);  // F is F3
    // map to domain [-1, 1]
    F = (F.array() < -1.0).select(-1.0, F);
    F = (F.array() > 1.0).select(1.0, F);
    // Stabilization step and convergence check
    calc_inbreed_coef(D1, F, PI, GL, type, nsnps, 0);  // F is F4
    diff = rmse1d(F0, F);
    cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE =", diff);
    if (diff < params.tolem) {
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    if (it == params.maxiter - 1) cao.warn("EM inbreeding coefficient not coverged!");
  }
  cao.print(tick.date(), "compute the LRT test");
  calc_inbreed_site_lrt(D1, F, PI, GL, type, nsnps, 0);
  // #pragma omp parallel for
  for (int j = 0; j < F.size(); j++) {
    D2(j) = chisq1d(D1(j));
  }
  write_hwe_per_site(params.fileout + ".hwe", params.filebim, D2, D1, F);
}

void calc_inbreed_coef_outofcore(Mat1D& D1, Mat1D& F, Data* data, Data* Pi, const Param& params) {
  data->check_file_offset_first_var();
  for (uint b = 0; b < data->nblocks; b++) {
    data->read_block_initial(data->start[b], data->stop[b], false);
    Pi->read_block_initial(Pi->start[b], Pi->stop[b], false);
    if (params.file_t == FileType::PLINK) {
      calc_inbreed_coef(D1, F, Pi->G, data->G, 1, Pi->stop[b] - Pi->start[b] + 1,
                        Pi->start[b]);  // init F
    } else {
      calc_inbreed_coef(D1, F, Pi->G, data->P, 2, Pi->stop[b] - Pi->start[b] + 1,
                        Pi->start[b]);  // init F
    }
  }
}

void run_inbreeding(Data* Pi, const Param& params) {
  Pi->prepare();
  Data* data = nullptr;
  if (params.file_t == FileType::PLINK) {
    data = new FileBed(params);
  } else if (params.file_t == FileType::BEAGLE) {
    data = new FileBeagle(params);
  } else {
    cao.error("input file not supported for estmating inbreeding coefficient");
  }
  data->prepare();
  assert(data->blocksize == Pi->blocksize);

  cao.print(tick.date(), "run inbreeding coefficient estimator");
  if (!params.out_of_core) {
    if (params.file_t == FileType::PLINK) run_inbreeding_em(1, data->G, Pi->G, params);
    if (params.file_t == FileType::BEAGLE) run_inbreeding_em(2, data->P, Pi->G, params);
    delete data;
    return;
  }
  // out of core run
  int nsnps = Pi->nsnps;
  Mat1D F = Mat1D::Zero(nsnps);  // init inbreeding coef
  Mat1D F0(nsnps), D1(nsnps), D2(nsnps);
  double sr2, sv2, alpha, diff;
  AreClose areClose;
  calc_inbreed_coef_outofcore(D1, F, data, Pi, params);  // init F
  for (uint it = 0; it < params.maxiter; it++) {
    F0 = F;                                                // copy the initial F
    calc_inbreed_coef_outofcore(D1, F, data, Pi, params);  // F is F1
    sr2 = D1.array().square().sum();
    calc_inbreed_coef_outofcore(D2, F, data, Pi, params);  // F is F2
    sv2 = (D2 - D1).array().square().sum();
    // safety break
    if (areClose(sv2, 0.0)) {
      cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE = 0.0");
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    alpha = fmin(fmax(1.0, sqrt(sr2 / sv2)), 256.0);
    if (params.verbose > 1) cao.print("alpha:", alpha, ", sr2:", sr2, ", sv2:", sv2);
    F = F0 + 2 * alpha * D1 + alpha * alpha * (D2 - D1);  // F is F3
    // map to domain [-1, 1]
    F = (F.array() < -1.0).select(-1.0, F);
    F = (F.array() > 1.0).select(1.0, F);
    // Stabilization step and convergence check
    calc_inbreed_coef_outofcore(D1, F, data, Pi, params);  // F is F4
    diff = rmse1d(F0, F);
    cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE =", diff);
    if (diff < params.tolem) {
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    if (it == params.maxiter - 1) cao.warn("EM inbreeding coefficient not coverged!");
  }
  cao.print(tick.date(), "compute the LRT test");
  data->check_file_offset_first_var();
  for (uint b = 0; b < data->nblocks; b++) {
    data->read_block_initial(data->start[b], data->stop[b], false);
    Pi->read_block_initial(Pi->start[b], Pi->stop[b], false);
    if (params.file_t == FileType::PLINK) {
      calc_inbreed_site_lrt(D1, F, Pi->G, data->G, 1, Pi->stop[b] - Pi->start[b] + 1, Pi->start[b]);
    } else {
      calc_inbreed_site_lrt(D1, F, Pi->G, data->P, 2, Pi->stop[b] - Pi->start[b] + 1, Pi->start[b]);
    }
  }
#pragma omp parallel for
  for (int j = 0; j < F.size(); j++) {
    D2(j) = chisq1d(D1(j));
  }
  write_hwe_per_site(params.fileout + ".hwe", params.filebim, D2, D1, F);

  delete data;
}
