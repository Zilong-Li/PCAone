/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileBeagle.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileBeagle.hpp"

using namespace std;

// read all data and estimate F
void FileBeagle::read_all() {
  fp = gzopen(params.filein.c_str(), "r");
  parse_beagle_file(P, fp, nsamples, nsnps);
  gzclose(fp);
  cao.print(tick.date(), "begin to estimate allele frequencies");
  F = Mat1D::Constant(nsnps, 0.25);
  {  // out of scope: eigen object will be released;
    Mat1D Ft = Mat1D::Zero(nsnps);
    double diff;
    // run EM to estimate allele frequencies
    for (uint it = 0; it < params.maxiter; it++) {
#pragma omp parallel for
      for (uint j = 0; j < nsnps; j++) {
        Ft(j) = F(j);
        double p0, p1, p2, pt = 0.0;
        for (uint i = 0; i < nsamples; i++) {
          p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
          p1 = P(2 * i + 1, j) * 2 * F(j) * (1.0 - F(j));
          p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
          pt += (p1 + 2 * p2) / (2 * (p0 + p1 + p2));
        }
        F(j) = pt / (double)nsamples;
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
