/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileBeagle.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileBeagle.hpp"

using namespace std;

// read all data and estimate F
void FileBeagle::read_all() {
  P = Mat2D::Zero(nsamples * 2, nsnps);
  fp = gzopen(params.filein.c_str(), "r");
  parse_beagle_file(P, fp, nsamples, nsnps);
  if (!params.pca) return;
  cao.print(tick.date(), "begin to estimate allele frequencies");
  F = Mat1D::Constant(nsnps, 0.25);
  {  // out of scope: eigen object will be released;
    Mat1D Ft = Mat1D::Zero(nsnps);
    double scale = 1.0 / (2.0 * nsamples);
    double diff;
    // run EM to estimate allele frequencies
    for (uint it = 0; it < params.maxiter; it++) {
#pragma omp parallel for
      for (uint j = 0; j < nsnps; j++) {
        Ft(j) = F(j);
        double p0, p1, p2, pt = 0.0;
        for (uint i = 0; i < nsamples; i++) {
          p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
          p1 = P(2 * i + 1, j) * 2.0 * F(j) * (1.0 - F(j));
          p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
          pt += (p1 + 2.0 * p2) / (p0 + p1 + p2);
        }
        F(j) = pt * scale;
      }
      // calculate differences between iterations
      diff = sqrt((F - Ft).array().square().sum() / nsnps);
      // Check for convergence
      if (diff < params.tolmaf) {
        cao.print(tick.date(), "EM (MAF) converged at iteration:", it + 1);
        break;
      } else if (it == (params.maxiter - 1)) {
        cao.print(tick.date(), "EM (MAF) did not converge");
      }
    }
  }
  // filter snps and resize G;
  filter_snps_resize_F();
  // resize P, only keep columns matching the indecis in idx;
  // P = P(Eigen::all, idx).eval(); // aliasing issue!!!
  G = Mat2D::Zero(nsamples, nsnps);  // initial E which is G
#pragma omp parallel for
  for (uint j = 0; j < nsnps; j++) {
    double p0, p1, p2;
    uint s = params.keepsnp ? keepSNPs[j] : j;
    for (uint i = 0; i < nsamples; i++) {
      p0 = P(2 * i + 0, s) * (1.0 - F(j)) * (1.0 - F(j));
      p1 = P(2 * i + 1, s) * 2 * F(j) * (1.0 - F(j));
      p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * F(j) * F(j);
      G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
    }
  }
}

void FileBeagle::check_file_offset_first_var() {
  if (params.verbose) cao.print("reopen beagle file and read head line");
  fp = gzopen(params.filein.c_str(), "r");
  tgets(fp, &buffer, &bufsize);  // parse header line
  if (buffer != original) original = buffer;
}

void FileBeagle::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {
  if (params.pca) cao.error("doesn't support out-of-core PCAngsd algorithm");
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    P = Mat2D::Zero(nsamples * 2, actual_block_size);
  }
  char* tok;
  // read all GL data into P
  for (uint j = 0; j < actual_block_size; ++j) {
    tgets(fp, &buffer, &bufsize);  // get a line
    if (buffer != original) original = buffer;
    tok = strtok_r(buffer, delims, &buffer);
    tok = strtok_r(NULL, delims, &buffer);
    tok = strtok_r(NULL, delims, &buffer);
    for (uint i = 0; i < nsamples; i++) {
      tok = strtok_r(NULL, delims, &buffer);
      P(2 * i + 0, j) = strtod(tok, NULL);
      tok = strtok_r(NULL, delims, &buffer);
      P(2 * i + 1, j) = strtod(tok, NULL);
      tok = strtok_r(NULL, delims, &buffer);
    }
    buffer = original;
  }
}
