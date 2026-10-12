/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileBeagle.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileBeagle.hpp"

#include "Data.hpp"
#include "Stats.hpp"

using namespace std;

// The GLs of a BEAGLE line (marker allele1 allele2, then AA AB BB of each
// sample) into column j of P: AA and AB of sample i in rows 2i and 2i + 1.
// BB is one minus the two.
static void parse_gl_line(char* line, Mat2D& P, uint j, uint nsamples) {
  const char* delims = "\t \n";
  char* save = nullptr;
  strtok_r(line, delims, &save);  // marker
  strtok_r(NULL, delims, &save);  // allele1
  strtok_r(NULL, delims, &save);  // allele2
  for (uint i = 0; i < nsamples; i++) {
    P(2 * i + 0, j) = strtod(strtok_r(NULL, delims, &save), NULL);
    P(2 * i + 1, j) = strtod(strtok_r(NULL, delims, &save), NULL);
    strtok_r(NULL, delims, &save);
  }
}

// read all data and estimate F
void FileBeagle::read_all() {
  uint total_nsnps = nsnps;
  P = Mat2D::Zero(nsamples * 2, total_nsnps);
  if (fp) gzclose(fp);  // the constructor's, left at the end of the file
  fp = gzopen(params.filein.c_str(), "r");
  if (!fp) cao.error("can not open " + params.filein);
  tgets(fp, &buffer, &bufsize);  // the header line
  uint j = 0;
  while (tgets(fp, &buffer, &bufsize)) {
    if (j >= total_nsnps)
      cao.error("more sites in " + params.filein + " than the " + to_string(total_nsnps) + " counted");
    parse_gl_line(buffer, P, j++, nsamples);
  }
  if (j != total_nsnps) cao.error("something wrong parsing beagle");
  original = buffer;  // tgets may have grown it
  // --inbreed pairs the GLs with the reference's individual allele frequencies
  // and needs nothing else: not the frequencies of this input, which the EM
  // below spent -maxiter passes on, nor the N x M expected genotypes
  if (params.inbreed != 0) return;
  if (!params.dopca && params.project > 0) {
    // projection mode: compute expected genotypes from GLs using reference F
    // F is read in Data::prepare() to include only overlapped sites
    const bool filter = !keepSNPs.empty();
    if (filter) nsnps = keepSNPs.size();
    G = Mat2D::Zero(nsamples, nsnps);
    C = ArrBool::Zero((uint64)nsnps * nsamples);

#pragma omp parallel for
    for (uint j = 0; j < nsnps; j++) {
      uint s = filter ? keepSNPs[j] : j;
      // const double norm = sqrt(2.0 * F(j) * (1.0 - F(j)));
      for (uint i = 0; i < nsamples; i++) {
        const double pt = fmin(fmax(F(j), 1e-4), 1.0 - 1e-4);
        const double p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        const double p1 = P(2 * i + 1, s) * 2.0 * pt * (1.0 - pt);
        const double p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        const double psum = p0 + p1 + p2;
        if (!std::isfinite(psum) || psum <= 0.0) {
          C[(uint64)j * nsamples + i] = 1;
          G(i, j) = 0.0;
          continue;
        }
        C[(uint64)j * nsamples + i] = 0;
        G(i, j) = (p1 + 2.0 * p2) / (2.0 * psum) - F(j);
        // if (params.scale == SCALE_STANDARDIZE_GENETIC && norm > VAR_TOL) G(i, j) /= norm;
      }
    }
    return;
  }
  cao.print(tick.date(), "begin to estimate allele frequencies");
  F = Mat1D::Constant(nsnps, 0.25);
  emMAF_with_GL(F, P, params.maxiter, params.tolmaf);
  // filter snps and resize G;
  filter_snps_resize_F();
  // resize P, only keep columns matching the indecis in idx;
  // P = P(Eigen::all, idx).eval(); // aliasing issue!!!
  G = Mat2D::Zero(nsamples, nsnps);  // initial E which is G
#pragma omp parallel for
  for (uint j = 0; j < nsnps; j++) {
    double p0, p1, p2;
    uint s = params.filterSNP ? keepSNPs[j] : j;
    for (uint i = 0; i < nsamples; i++) {
      p0 = P(2 * i + 0, s) * (1.0 - F(j)) * (1.0 - F(j));
      p1 = P(2 * i + 1, s) * 2 * F(j) * (1.0 - F(j));
      p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * F(j) * F(j);
      G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
    }
  }
}

void FileBeagle::check_file_offset_first_var() {
  if (params.verbose > 1) cao.print("reopen beagle file and read head line");
  // close the previous pass's handle: one leaked per pass, and the EM of
  // --inbreed makes three passes per iteration, past a 256 file limit
  if (fp) gzclose(fp);
  fp = gzopen(params.filein.c_str(), "r");
  if (!fp) cao.error("can not open " + params.filein);
  tgets(fp, &buffer, &bufsize);  // parse header line
  if (buffer != original) original = buffer;
}

void FileBeagle::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {
  if (params.dopca) cao.error("doesn't support out-of-core PCAngsd algorithm");
  uint actual_block_size = stop_idx - start_idx + 1;
  // P, not G: G is never allocated here, so P was reallocated for every block
  if (P.cols() != (Eigen::Index)actual_block_size) P = Mat2D::Zero(nsamples * 2, actual_block_size);
  for (uint j = 0; j < actual_block_size; ++j) {
    tgets(fp, &buffer, &bufsize);  // get a line
    parse_gl_line(buffer, P, j, nsamples);
  }
  original = buffer;  // tgets may have grown it
}
