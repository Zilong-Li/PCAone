/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileBeagle.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileBeagle.hpp"

using namespace std;

// inspired by Angsd
int tgets(gzFile gz, char **buf, uint64 *l) {
  int rlen = 0;
neverUseGoto:
  char *tok = gzgets(gz, *buf + rlen, *l - rlen);
  if (!tok) return rlen;
  int tmp = tok ? strlen(tok) : 0;
  if (tok[tmp - 1] != '\n') {
    rlen += tmp;
    *l *= 2;
    *buf = (char *)realloc(*buf, *l);
    goto neverUseGoto;
  }
  rlen += tmp;
  return rlen;
}

// read all data and estimate F
void FileBeagle::read_all() {
  fp = gzopen(params.filein.c_str(), "r");
  tgets(fp, &buffer, &bufsize);
  char *tok;
  uint i = 0, j = 0;
  // read all GP data into P
  while (tgets(fp, &buffer, &bufsize)) {
    if (buffer != original) original = buffer;
    tok = strtok_r(buffer, delims, &buffer);
    tok = strtok_r(NULL, delims, &buffer);
    tok = strtok_r(NULL, delims, &buffer);
    for (i = 0; i < nsamples; i++) {
      tok = strtok_r(NULL, delims, &buffer);
      P(2 * i + 0, j) = strtod(tok, NULL);
      tok = strtok_r(NULL, delims, &buffer);
      P(2 * i + 1, j) = strtod(tok, NULL);
      tok = strtok_r(NULL, delims, &buffer);
    }
    buffer = original;
    j++;
  }
  gzclose(fp);
  assert(j == nsnps);

  cao.print(tick.date(), "begin to estimate allele frequencies");
  F = MyVector::Constant(nsnps, 0.25);
  {  // out of scope: eigen object will be released;
    MyVector Ft = MyVector::Zero(nsnps);
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
  G = MyMatrix::Zero(nsamples, nsnps);  // initial E which is G
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
