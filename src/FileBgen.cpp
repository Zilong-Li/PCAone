/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileBgen.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileBgen.hpp"

#include <algorithm>
#include <numeric>

using namespace std;

void FileBgen::read_all() {
  uint i, j;
  cao.print(tick.date(), "start reading all data");
  if (!params.pcangsd) {
    F = Mat1D::Zero(nsnps);
    G = Mat2D(nsamples, nsnps);
    if (params.missme) C = ArrBool::Zero((uint64)nsnps * nsamples);
    // each thread decodes whole variants into their columns; the variants
    // that fail --maf are dropped afterwards
    std::vector<char> keep(nsnps, 0);
    std::string err;
#pragma omp parallel for private(i, j) schedule(dynamic, 8)
    for (j = 0; j < nsnps; j++) {
      const int thr = omp_get_thread_num();
      float* dose = thread_dosages[thr].data();
      try {
        reader->minor_dosage(thr, j, dose);
      } catch (const std::exception& e) {
#pragma omp critical
        if (err.empty()) err = e.what();
        continue;
      }
      uint64 gc = 0;
      double gs = 0.0;
      for (i = 0; i < nsamples; i++) {
        if (!std::isnan(dose[i])) {
          gs += dose[i] / 2.0;  // map to [0, 1];
          gc += 1;
        }
      }
      const double af = gc == 0 ? 0.0 : gs / gc;
      if (!(af > params.maf)) continue;
      keep[j] = 1;
      F(j) = af;
      // do centering and initialing
      for (i = 0; i < nsamples; i++) {
        if (std::isnan(dose[i])) {
          G(i, j) = 0;
          if (params.missme) C[(uint64)j * nsamples + i] = 1;
        } else {
          G(i, j) = dose[i] / 2.0 - af;  // map to [0, 1];
        }
      }
    }
    if (!err.empty()) cao.error(err);
    // move the kept variants to the front, in order
    uint64 k = 0;
    for (j = 0; j < nsnps; j++) {
      if (!keep[j]) continue;
      if (k != j) {
        G.col(k) = G.col(j);
        F(k) = F(j);
        if (params.missme)
          std::copy_n(C.data() + (uint64)j * nsamples, nsamples, C.data() + k * nsamples);
      }
      k++;
    }
    if (k == 0)
      cao.error("the number of SNPs after filtering is 0!");
    else
      cao.print(tick.date(), "number of SNPs after filtering by MAF >", params.maf, ":", k);
    // resize G, F, C;
    nsnps = k;  // resize nsnps;
    G.conservativeResize(Eigen::NoChange, nsnps);
    F.conservativeResize(nsnps);
    if (params.missme) C.conservativeResize((uint64)nsnps * nsamples);
  } else {
    // read all GP data into P;
    for (j = 0; j < nsnps; j++) {
      try {
        auto var = bg->next_var();
        probs1d.resize(nsamples * var.probs_per_sample());
        var.probs_1d(probs1d.data());
#pragma omp parallel for
        for (i = 0; i < nsamples; i++) {
          P(i * 2 + 0, j) = probs1d[i * 3 + 0];
          P(i * 2 + 1, j) = probs1d[i * 3 + 1];
          // no need to parse probs1d[i * 3 + 2]
        }
      } catch (const std::out_of_range& e) {
        throw e.what();
      }
    }
    assert(j == nsnps);
    cao.print(tick.date(), "begin to estimate allele frequencies using GP");
    F = Mat1D::Constant(nsnps, 0.25);
    emMAF_with_GL(F, P, params.maxiter, params.tolmaf);
    filter_snps_resize_F();
    // initial E which is G
    G = Mat2D::Zero(nsamples, nsnps);
#pragma omp parallel for
    for (j = 0; j < nsnps; j++) {
      double p0, p1, p2;
      uint s = params.filterSNP ? keepSNPs[j] : j;
      for (i = 0; i < nsamples; i++) {
        p0 = P(2 * i + 0, s) * (1.0 - F(j)) * (1.0 - F(j));
        p1 = P(2 * i + 1, s) * 2 * F(j) * (1.0 - F(j));
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * F(j) * F(j);
        G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
      }
    }
  }
}

void FileBgen::variants_of(uint64 first, uint64 count, std::vector<uint32_t>& out) const {
  out.resize(count);
  for (uint64 j = 0; j < count; ++j) out[j] = params.perm ? (uint32_t)perm.indices()(first + j) : (uint32_t)(first + j);
  if (params.perm) std::sort(out.begin(), out.end());
}

void FileBgen::begin_block(uint64 start_idx, uint64 stop_idx) {
  const uint64 count = stop_idx - start_idx + 1;
  requests.resize(count);
  for (uint64 j = 0; j < count; ++j)
    requests[j] = {params.perm ? (uint32_t)perm.indices()(start_idx + j) : (uint32_t)(start_idx + j), (uint32_t)j};
  if (params.perm)
    std::sort(requests.begin(), requests.end(),
              [](const ReadRequest& a, const ReadRequest& b) { return a.variant < b.variant; });
  if (params.noprefetch) return;
  // request the next block of a pass, or the first block after the last
  const bool sequential = start_idx == 0 || start_idx == last_end;
  last_end = stop_idx + 1;
  if (start_idx == 0) first_count = count;
  if (sequential && count < nsnps) {
    uint64 next = stop_idx + 1, n = std::min<uint64>(count, nsnps - std::min<uint64>(next, nsnps));
    if (next >= nsnps) next = 0, n = std::min<uint64>(first_count ? first_count : count, nsnps);
    variants_of(next, n, next_variants);
    try {
      reader->prefetch(next_variants);
    } catch (const std::exception&) {
      // only a hint
    }
  }
}

void FileBgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {
  const uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = Mat2D::Zero(nsamples, actual_block_size);
  }
  begin_block(start_idx, stop_idx);
  const bool genetic_sd = standardize && params.scale == SCALE_STANDARDIZE_GENETIC;
  const bool estimate = !frequency_was_estimated;
  std::string err;
  uint r, j;
#pragma omp parallel for private(r, j) schedule(dynamic, 4)
  for (r = 0; r < actual_block_size; ++r) {
    const uint i = requests[r].column;
    const uint64 snp_idx = start_idx + i;
    const int thr = omp_get_thread_num();
    float* dose = thread_dosages[thr].data();
    try {
      reader->minor_dosage(thr, requests[r].variant, dose);
    } catch (const std::exception& e) {
#pragma omp critical
      if (err.empty()) err = e.what();
      continue;
    }
    if (estimate) {
      uint64 gc = 0;
      double gs = 0.0;
      for (j = 0; j < nsamples; j++) {
        if (!std::isnan(dose[j])) {
          gs += dose[j] / 2.0;
          gc += 1;
        }
      }
      if (gc == 0) {
#pragma omp critical
        if (err.empty()) err = "the allele frequency should not be 0. do filtering first";
        continue;
      }
      F(snp_idx) = gs / gc;
    }
    const double f = F(snp_idx);
    double scale = 1.0;
    if (genetic_sd) {
      const double sd = sqrt(f * (1 - f));
      if (sd > VAR_TOL) scale = sqrt((double)params.ploidy) / sd;
    }
    double* g = &G(0, i);
    for (j = 0; j < nsamples; j++) g[j] = std::isnan(dose[j]) ? 0.0 : (dose[j] / 2.0 - f) * scale;
  }
  if (!err.empty()) cao.error(err);
  if (estimate && stop_idx + 1 == nsnps) frequency_was_estimated = true;
}

PermMat compute_bgen_perm(uint nsnps, int seed) {
  std::vector<int> order(nsnps);
  std::iota(order.begin(), order.end(), 0);
  PortableRng rng(seed);  // --seed
  portable_shuffle(order.begin(), order.end(), rng);
  PermMat P;
  P.indices() = Eigen::Map<Eigen::VectorXi>(order.data(), order.size());
  return P;
}
