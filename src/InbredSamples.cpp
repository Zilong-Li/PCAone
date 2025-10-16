/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Inbreeding.cpp
 * @author      Zilong Li
 * Copyright (C) 2025. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "InbredSamples.hpp"

#include "Cmd.hpp"
#include "Common.hpp"
#include "FileBeagle.hpp"
#include "FilePlink.hpp"
#include "Utils.hpp"

// type 1: Genotype input, {1, -9, 0.5, 0}, GL is N x M
// type 2: Genotype likelihood input, GL is (N x 2) x M
void inbreed_coef_sample(Mat1D& D, Mat1D& F, const Mat2D& PI, const Mat2D& GL, const int type, const int size,
                         const uint start) {
  const int nsnps = size;
  const int nsamples = PI.rows();
  if (type != 1 && type != 2) cao.error("type must be 1 or 2");
  Double1D obsHall(nsamples, 0.0), expHall(nsamples, 0.0);
#pragma omp parallel
  for (int j = 0; j < nsnps; j++) {
    double Fadj, p0, p1, p2, pSum, pp0, pp1, pp2;
    Double1D obsH(nsamples, 0.0), expH(nsamples, 0.0);
    // int jj = j + start;
#pragma omp for
    for (int i = 0; i < nsamples; i++) {
      // Fadj = (1-pi)*pi*F;
      Fadj = (1.0 - PI(i, j)) * PI(i, j) * F(i);
      // (1-pi)^2 + pi(1-pi)*F
      p0 = fmax(1e-4, (1.0 - PI(i, j)) * (1.0 - PI(i, j)) + Fadj);
      // 2pi(1-pi)(1-F)
      p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) * (1.0 - F(i)));
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
        if (GL(i, j) == BED2GENO[3]) {
          pp0 = p0;
          pp1 = 0;
          pp2 = 0;
        } else if (GL(i, j) == BED2GENO[2]) {
          pp0 = 0;
          pp1 = p1;
          pp2 = 0;
        } else if (GL(i, j) == BED2GENO[0]) {
          pp0 = 0;
          pp1 = 0;
          pp2 = p2;
        } else {
          pp0 = 0.333333 * p0;
          pp1 = 0.333333 * p1;
          pp2 = 0.333333 * p2;
        }
      } else {
        pp0 = GL(2 * i + 0, j) * p0;
        pp1 = GL(2 * i + 1, j) * p1;
        pp2 = (1.0 - GL(2 * i + 0, j) - GL(2 * i + 1, j)) * p2;
      }
      pSum = pp0 + pp1 + pp2;

      // sum over indiiduals
      obsH[i] += pp1 / pSum;
      // count heterozygotes
      expH[i] += 2.0 * PI(i, j) * (1.0 - PI(i, j));
    }
#pragma omp critical
    {
      for (int i = 0; i < nsamples; i++) {
        obsHall[i] += obsH[i];
        expHall[i] += expH[i];
      }
    }
  }

  // Update the inbreeding coefficient
#pragma omp parallel for
  for (int i = 0; i < nsamples; i++) {
    double f = fmin(fmax(-1.0, 1.0 - (obsHall[i] / expHall[i])), 1.0);
    D(i) = f - F(i);  // Fnew - Fold
    F(i) = f;         // update with Fnew
  }
}

void inbreed_coef_sample_em(int type, const Mat2D& GL, const Mat2D& PI, const Param& params) {
  if (type != 1 && type != 2) cao.error("type must be 1 or 2");
  const int nsnps = PI.cols();
  const int nsamples = PI.rows();
  Mat1D F = Mat1D::Zero(nsamples), F0(nsamples);  // init inbreeding coef
  Mat1D D1(nsamples), D2(nsamples);               // store the diff between Fnew - Fold
  double sr2, sv2, alpha, diff;
  AreClose areClose;
  for (uint it = 0; it < params.maxiter; it++) {
    F0 = F;                                              // copy the initial F
    inbreed_coef_sample(D1, F, PI, GL, type, nsnps, 0);  // F is F1, D1 = F1 - F0
    sr2 = D1.array().square().sum();
    inbreed_coef_sample(D2, F, PI, GL, type, nsnps, 0);  // F is F2, D2 = F2 - F1
    sv2 = (D2 - D1).array().square().sum();
    // safety break
    if (areClose(sv2, 0.0)) {
      cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE = 0");
      cao.print(tick.date(), "EM inbreeding coefficient coverged!");
      break;
    }
    alpha = fmin(fmax(1.0, sqrt(sr2 / sv2)), 256.0);
    if (params.verbose > 1) cao.print("alpha:", alpha, ", sr2:", sr2, ", sv2:", sv2);
    F = F0 + 2.0 * alpha * D1 + alpha * alpha * (D2 - D1);  // F is F3
    // map to domain [-1, 1]
    F = (F.array() < -1.0).select(-1.0, F);
    F = (F.array() > 1.0).select(1.0, F);
    // Stabilization step and convergence check
    inbreed_coef_sample(D1, F, PI, GL, type, nsnps, 0);  // F is F4
    diff = rmse1d(F0, F);
    cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE =", diff);
    if (diff < params.tolem) {
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    if (it == params.maxiter - 1) cao.warn("EM inbreeding coefficient not coverged!");
  }
  write_coef_per_sample(params.fileout + ".inbred", params.filein + ".fam", F);
}

void calc_inbreed_coef_outofcore(Mat1D& D1, Mat1D& F, Data* data, Data* Pi, const Param& params) {
  data->check_file_offset_first_var();
  for (uint b = 0; b < data->nblocks; b++) {
    data->read_block_initial(data->start[b], data->stop[b], false);
    Pi->read_block_initial(Pi->start[b], Pi->stop[b], false);
    if (params.file_t == FileType::PLINK) {
      inbreed_coef_sample(D1, F, Pi->G, data->G, 1, Pi->stop[b] - Pi->start[b] + 1,
                          Pi->start[b]);  // init F
    } else {
      inbreed_coef_sample(D1, F, Pi->G, data->P, 2, Pi->stop[b] - Pi->start[b] + 1,
                          Pi->start[b]);  // init F
    }
  }
}

void run_inbreed_coef_sample(Data* Pi, const Param& params) {
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

  cao.print(tick.date(), "run inbreeding coefficient estimator per sample");
  if (!params.out_of_core) {
    if (params.file_t == FileType::PLINK) inbreed_coef_sample_em(1, data->G, Pi->G, params);
    if (params.file_t == FileType::BEAGLE) inbreed_coef_sample_em(2, data->P, Pi->G, params);
    delete data;
    return;
  }
  // out of core run
  int nsamples = Pi->nsamples;
  Mat1D F = Mat1D::Zero(nsamples);  // init inbreeding coef
  Mat1D F0(nsamples), D1(nsamples), D2(nsamples);
  double sr2, sv2, alpha, diff;
  AreClose areClose;
  for (uint it = 0; it < params.maxiter; it++) {
    F0 = F;                                                // copy the initial F
    calc_inbreed_coef_outofcore(D1, F, data, Pi, params);  // F is F1
    sr2 = D1.array().square().sum();
    calc_inbreed_coef_outofcore(D2, F, data, Pi, params);  // F is F2
    sv2 = (D2 - D1).array().square().sum();
    // safety break
    if (areClose(sv2, 0.0)) {
      cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE = 0");
      cao.print(tick.date(), "EM inbreeding coefficient coverged!");
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

  write_coef_per_sample(params.fileout + ".inbred", params.filein + ".fam", F);

  delete data;
}

void write_coef_per_sample(const std::string& fout, const std::string& fam, const Mat1D& coef) {
  std::ifstream fin(fam);
  if (!fin.is_open()) cao.error("can not open " + fam);
  std::ofstream ofs(fout);
  if (!ofs.is_open()) cao.error("can not open " + fout);
  ofs << "#ID\tInbreeding_coefficient\n";
  const std::string sep{" \t"};
  std::string line;
  int j = 0;
  while (getline(fin, line)) {
    auto id = split_string(line, sep)[1];
    ofs << id << "\t" << coef(j) << "\n";
    j++;
  }
  assert(j == coef.size());
}
