/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Inbreeding.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "Inbreeding.hpp"

#include "Cmd.hpp"
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
    double obsH = 0.0, expH = 0.0, Fadj, p0, p1, p2, pSum, pp0, pp1, pp2, ppSum;
    int jj = j + start;
    for (int i = 0; i < nsamples; i++) {
      // Fadj = (1-pi)*pi*F;
      Fadj = (1.0 - PI(i, j)) * PI(i, j) * F(jj);
      // (1-pi)^2 + pi(1-pi)*F
      p0 = fmax(1e-4, (1.0 - PI(i, j) * (1.0 - PI(i, j)) + Fadj));
      // 2pi(1-pi)(1-F)
      p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) - 2.0 * Fadj);
      // pi^2 + pi(1-pi)*F
      p2 = fmax(1e-4, PI(i, j) * PI(i, j) + Fadj);
      // normalize
      pSum = p0 + p1 + p2;
      p0 /= pSum;
      p1 /= pSum;
      p2 /= pSum;

      // posterior
      if (type == 1) {
        // should check genotype missingness
        if (GL(i, j) == 0.0) {
          pp0 = p0, pp1 = 0, pp2 = 0;
        }
        if (GL(i, j) == 0.5) {
          pp0 = 0, pp1 = p1, pp2 = 0;
        }
        if (GL(i, j) == 1.0) {
          pp0 = 0, pp1 = 0, pp2 = p2;
        }
      } else {
        pp0 = GL(2 * i + 0, j) * p0;
        pp1 = GL(2 * i + 1, j) * p1;
        pp2 = (1.0 - GL(2 * i + 0, j) - GL(2 * i + 1, j)) * p2;
      }
      ppSum = pp0 + pp1 + pp2;

      // sum over indiiduals
      obsH += pp1 / ppSum;
      // count heterozygotes
      expH += 2.0 * PI(i, j) * (1.0 - PI(i, j));
    }
    // ANGSD procedure
    obsH = fmax(1e-4, obsH / (double)nsamples);
    // Update the inbreeding coefficient
    double f = 1.0 - ((double)nsamples * obsH / expH);
    f = fmin(fmax(-1.0, f), 1.0);
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
      p0 = fmax(1e-4, (1.0 - PI(i, j) * (1.0 - PI(i, j)) + Fadj));
      // 2pi(1-pi)(1-F)
      p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) - 2.0 * Fadj);
      // pi^2 + pi(1-pi)*F
      p2 = fmax(1e-4, PI(i, j) * PI(i, j) + Fadj);
      // normalize
      pSum = p0 + p1 + p2;
      p0 /= pSum;
      p1 /= pSum;
      p2 /= pSum;

      // posterior // Likelihood*prior
      if (type == 1) {
        if (GL(i, j) == 0.0) {
          logAlt += log(p0);
          logNull += log((1.0 - PI(i, j)) * (1.0 - PI(i, j)));
        }
        if (GL(i, j) == 0.5) {
          logAlt += log(p1);
          logNull += log(2.0 * PI(i, j) * (1.0 - PI(i, j)));
        }
        if (GL(i, j) == 1.0) {
          logAlt += log(p2);
          logNull += log(PI(i, j) * PI(i, j));
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
  Mat1D F = Mat1D::Zero(nsnps);  // init inbreeding coef
  Mat1D F0(nsnps), D1(nsnps), D2(nsnps);
  double sr2, sv2, alpha, diff;
  AreClose areClose;
  for (uint it = 0; it < params.maxiter; it++) {
    F0 = F;                                            // copy the initial F
    calc_inbreed_coef(D1, F, PI, GL, type, nsnps, 0);  // F is F1
    sr2 = D1.array().square().sum();
    calc_inbreed_coef(D2, F, PI, GL, type, nsnps, 0);  // F is F2
    sv2 = (D2 - D1).array().square().sum();
    // safety break
    if (areClose(sv2, 0.0)) {
      cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE = 0.0");
      cao.print(tick.date(), "EM inbreeding coefficient coverged");
      break;
    }
    alpha = fmax(1.0, sqrt(sr2 / sv2));
    if (params.verbose) cao.print("alpha:", alpha, sr2, sv2);
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
#pragma omp parallel for
  for (int j = 0; j < F.size(); j++) {
    D2(j) = chisq1d(D1(j));
  }
  write_hwe_per_site(params.fileout + ".hwe", params.filebim, D2, D1, F);
}

void run_inbreeding(Data* data, const Param& params) {
  cao.print(tick.date(), "run inbreeding coefficient estimator");
  data->prepare();
  if (params.file_t == FileType::BEAGLE) {
    Mat2D G(data->nsamples * 2, data->nsnps);  // genotype likelihood
    gzFile fp = gzopen(params.filein.c_str(), "r");
    parse_beagle_file(G, fp, data->nsamples, data->nsnps);
    gzclose(fp);
    // run EM-HWE
    run_inbreeding_em(2, G, data->G, params);
  }
  if (params.file_t == FileType::PLINK) {
    FileBed* geno = new FileBed(params);
    geno->prepare();
    assert(geno->blocksize == data->blocksize);
    // run EM-HWE
    if (!params.out_of_core) {
      run_inbreeding_em(1, geno->G, data->G, params);
    } else {
      // out of core run
      int nsnps = data->nsnps;
      Mat1D F = Mat1D::Zero(nsnps);  // init inbreeding coef
      Mat1D F0(nsnps), D1(nsnps), D2(nsnps);
      double sr2, sv2, alpha, diff;
      AreClose areClose;
      for (uint it = 0; it < params.maxiter; it++) {
        F0 = F;  // copy the initial F
        geno->check_file_offset_first_var();
        for (uint b = 0; b < geno->nblocks; b++) {
          geno->read_block_initial(geno->start[b], geno->stop[b], false);
          data->read_block_initial(data->start[b], data->stop[b], false);
          calc_inbreed_coef(D1, F, data->G, geno->G, 1, data->stop[b] - data->start[b] + 1,
                            data->start[b]);  // F is F1
        }
        sr2 = D1.array().square().sum();
        geno->check_file_offset_first_var();
        for (uint b = 0; b < geno->nblocks; b++) {
          geno->read_block_initial(geno->start[b], geno->stop[b], false);
          data->read_block_initial(data->start[b], data->stop[b], false);
          calc_inbreed_coef(D2, F, data->G, geno->G, 1, data->stop[b] - data->start[b] + 1,
                            data->start[b]);  // F is F2
        }
        sv2 = (D2 - D1).array().square().sum();
        // safety break
        if (areClose(sv2, 0.0)) {
          cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, "RMSE = 0.0");
          cao.print(tick.date(), "EM inbreeding coefficient coverged");
          break;
        }
        alpha = fmax(1.0, sqrt(sr2 / sv2));
        if (params.verbose) cao.print("alpha:", alpha, sr2, sv2);
        F = F0 + 2 * alpha * D1 + alpha * alpha * (D2 - D1);  // F is F3
        // map to domain [-1, 1]
        F = (F.array() < -1.0).select(-1.0, F);
        F = (F.array() > 1.0).select(1.0, F);
        // Stabilization step and convergence check
        geno->check_file_offset_first_var();
        for (uint b = 0; b < geno->nblocks; b++) {
          geno->read_block_initial(geno->start[b], geno->stop[b], false);
          data->read_block_initial(data->start[b], data->stop[b], false);
          calc_inbreed_coef(D1, F, data->G, geno->G, 1, data->stop[b] - data->start[b] + 1,
                            data->start[b]);  // F is F4
        }
        diff = rmse1d(F0, F);
        cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE =", diff);
        if (diff < params.tolem) {
          cao.print(tick.date(), "EM inbreeding coefficient coverged");
          break;
        }
        if (it == params.maxiter - 1) cao.warn("EM inbreeding coefficient not coverged!");
      }
      cao.print(tick.date(), "compute the LRT test");
      geno->check_file_offset_first_var();
      for (uint b = 0; b < geno->nblocks; b++) {
        geno->read_block_initial(geno->start[b], geno->stop[b], false);
        data->read_block_initial(data->start[b], data->stop[b], false);
        calc_inbreed_site_lrt(D1, F, data->G, geno->G, 1, data->stop[b] - data->start[b] + 1, data->start[b]);
      }
#pragma omp parallel for
      for (int j = 0; j < F.size(); j++) {
        D2(j) = chisq1d(D1(j));
      }
      write_hwe_per_site(params.fileout + ".hwe", params.filebim, D2, D1, F);
    }
    delete geno;
  }
}
