/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FilePgen.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FilePgen.hpp"

#include <omp.h>

#include "Common.hpp"
#include "Utils.hpp"

using namespace std;

// ReadHardcalls(allele_idx=1) returns counts of ALT allele: 0.0=HomRef,
// 1.0=Het, 2.0=HomAlt, -3.0=missing.  Divide by 2 to get [0,1] dosage
// matching the BED convention (BED2GENO = {1,-9,0.5,0} for A1=ALT).
static constexpr double PGEN_MISSING = -3.0;
static inline double pgen2dosage(double v) { return (v == PGEN_MISSING) ? BED_MISSING_VALUE : v / 2.0; }

static inline double centered_pgen_value(double v, double af) {
  double dosage = pgen2dosage(v);
  return (dosage == BED_MISSING_VALUE) ? 0.0 : dosage - af;
}

void FilePgen::read_all() {
  uint i, j;
  if (params.dopca && !frequency_was_estimated && !params.filterSNP) {
    // One pass: decode each variant into G, take F from that column, centre it.
    // The two-pass path below decodes every variant twice (once for F, once for
    // G), which a --maf filter needs; without one it is ~1.4x the read time.
    F = Mat1D::Zero(nsnps);
    G = Mat2D::Zero(nsamples, nsnps);
    if (params.missme) C = ArrBool::Zero((uint64)nsnps * nsamples);
    uint64 nmono = 0;
#pragma omp parallel for private(i, j) schedule(static) reduction(+ : nmono)
    for (i = 0; i < nsnps; ++i) {
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      if (dosage_mode) {
        reader.Read(buf, nsamples, thr, i, 1);
      } else {
        reader.ReadHardcalls(buf, nsamples, thr, i, 1);
      }
      uint64 c = 0;
      double sum = 0.0;
      for (j = 0; j < nsamples; ++j) {
        if (buf[j] != PGEN_MISSING) {
          sum += buf[j] / 2.0;
          ++c;
        }
      }
      F(i) = (c > 0) ? sum / c : 0.0;
      if (F(i) == 0.0 || F(i) == 1.0) ++nmono;
      for (j = 0; j < nsamples; ++j) {
        G(j, i) = pgen2dosage(buf[j]);
        if (params.missme && G(j, i) == BED_MISSING_VALUE) C[(uint64)i * nsamples + j] = 1;
        if (params.center) G(j, i) = (G(j, i) == BED_MISSING_VALUE) ? 0.0 : G(j, i) - F(i);
      }
    }
    warn_monomorphic(nmono);
    if (params.missme) {
      p_miss = (double)C.count() / (double)C.size();
      cao.print(tick.date(), "the proportion of missingness  is", p_miss);
    }
    return;
  }
  if (params.dopca && !frequency_was_estimated) {
    F = Mat1D::Zero(nsnps);
    uint64 nmono = 0;
#pragma omp parallel for private(i, j) schedule(static) reduction(+ : nmono)
    for (i = 0; i < nsnps; ++i) {
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      if (dosage_mode) {
        reader.Read(buf, nsamples, thr, i, 1);
      } else {
        reader.ReadHardcalls(buf, nsamples, thr, i, 1);
      }
      uint64 c = 0;
      double sum = 0.0;
      for (j = 0; j < nsamples; ++j) {
        if (buf[j] != PGEN_MISSING) {
          sum += buf[j] / 2.0;
          ++c;
        }
      }
      F(i) = (c > 0) ? sum / c : 0.0;
      if (F(i) == 0.0 || F(i) == 1.0) ++nmono;  // counted: cao is not thread-safe
    }
    warn_monomorphic(nmono);
    filter_snps_resize_F();
  }
  const bool filter = !keepSNPs.empty();
  if (filter) nsnps = keepSNPs.size();
  G = Mat2D::Zero(nsamples, nsnps);

  if (params.missme) C = ArrBool::Zero((uint64)nsnps * nsamples);

#pragma omp parallel for private(i, j) schedule(static)
  for (i = 0; i < nsnps; ++i) {
    int thr = omp_get_thread_num();
    double* buf = thread_bufs[thr].data();
    uint s = filter ? keepSNPs[i] : i;
    if (dosage_mode) {
      reader.Read(buf, nsamples, thr, s, 1);
    } else {
      reader.ReadHardcalls(buf, nsamples, thr, s, 1);
    }
    for (j = 0; j < nsamples; ++j) {
      G(j, i) = pgen2dosage(buf[j]);
      if ((params.missme && G(j, i) == BED_MISSING_VALUE)) C[(uint64)i * nsamples + j] = 1;
    }
    if (params.center) {
      for (j = 0; j < nsamples; ++j) {
        if (G(j, i) == BED_MISSING_VALUE)
          G(j, i) = 0.0;  // impute to mean
        else
          G(j, i) -= F(i);
      }
    }
  }

  if (params.missme) {
    p_miss = (double)C.count() / (double)C.size();
    cao.print(tick.date(), "the proportion of missingness  is", p_miss);
  }
}

PermMat compute_pgen_perm(uint nsnps, uint nbatches, uint blocksize, uint nthreads, int seed) {
  if (nsnps == 0) return PermMat(0);

  nbatches = std::max<uint>(1, std::min<uint>(nbatches, nsnps));
  blocksize = std::max<uint>(1, blocksize);
  nthreads = std::max<uint>(1, std::min<uint>(nthreads, std::min<uint>(blocksize, nsnps)));

  std::vector<uint64> read_count_by_thread(nthreads, 0);
  for (uint batch = 0; batch < nbatches; ++batch) {
    uint64 batch_start = ((uint64)batch * nsnps) / nbatches;
    uint64 batch_stop = ((uint64)(batch + 1) * nsnps) / nbatches;
    for (uint64 block_start = batch_start; block_start < batch_stop; block_start += blocksize) {
      uint block_len = (uint)std::min<uint64>(blocksize, batch_stop - block_start);
      uint active_threads = std::min<uint>(nthreads, block_len);
      uint base = block_len / active_threads;
      uint extra = block_len % active_threads;
      for (uint t = 0; t < active_threads; ++t) {
        read_count_by_thread[t] += base + (t < extra ? 1u : 0u);
      }
    }
  }

  std::vector<std::vector<int>> source_by_thread(nthreads);
  uint64 source_start = 0;
  for (uint t = 0; t < nthreads; ++t) {
    uint64 source_stop = source_start + read_count_by_thread[t];
    source_by_thread[t].reserve(source_stop - source_start);
    for (uint64 snp_idx = source_start; snp_idx < source_stop; ++snp_idx) {
      source_by_thread[t].push_back((int)snp_idx);
    }
    source_start = source_stop;
  }

  PortableRng rng(seed);
  for (auto& source : source_by_thread) portable_shuffle(source.begin(), source.end(), rng);

  std::vector<uint64> next_by_thread(nthreads, 0);
  Eigen::VectorXi indices(nsnps);

  for (uint batch = 0; batch < nbatches; ++batch) {
    uint64 batch_start = ((uint64)batch * nsnps) / nbatches;
    uint64 batch_stop = ((uint64)(batch + 1) * nsnps) / nbatches;
    for (uint64 block_start = batch_start; block_start < batch_stop; block_start += blocksize) {
      uint block_len = (uint)std::min<uint64>(blocksize, batch_stop - block_start);
      uint active_threads = std::min<uint>(nthreads, block_len);
      uint base = block_len / active_threads;
      uint extra = block_len % active_threads;
      uint offset = 0;
      for (uint t = 0; t < active_threads; ++t) {
        uint len = base + (t < extra ? 1u : 0u);
        for (uint j = 0; j < len; ++j) {
          if (next_by_thread[t] >= source_by_thread[t].size()) cao.error("BUG: exhausted PGEN permutation partition.");
          indices((Eigen::Index)(block_start + offset + j)) = source_by_thread[t][next_by_thread[t]++];
        }
        offset += len;
      }
    }
  }

  return PermMat(indices);
}

FilePgen::~FilePgen() {
  if (block_reader) {
    block_reader->cancel();
    if (params.verbose > 1)
      cao.print(tick.date(), "PGEN blocks predicted:", block_reader->predicted(),
                ", requested from the disk:", block_reader->requested(), ", page cache probe:",
                block_reader->probe_name());
  }
}

// the file indices of the logical SNPs [first, first + count), sorted
void FilePgen::variants_of(uint64 first, uint64 count, std::vector<uint32_t>& out) const {
  out.resize(count);
  for (uint64 j = 0; j < count; ++j) out[j] = params.perm ? (uint32_t)perm.indices()(first + j) : (uint32_t)(first + j);
  if (params.perm) std::sort(out.begin(), out.end());
}

void FilePgen::begin_block(uint64 start_idx, uint64 stop_idx) {
  const uint64 count = stop_idx - start_idx + 1;
  requests.resize(count);
  for (uint64 j = 0; j < count; ++j)
    requests[j] = {params.perm ? (uint32_t)perm.indices()(start_idx + j) : (uint32_t)(start_idx + j), (uint32_t)j};
  if (params.perm)
    std::sort(requests.begin(), requests.end(),
              [](const ReadRequest& a, const ReadRequest& b) { return a.variant < b.variant; });

  if (!block_reader_tried) {
    block_reader_tried = true;
    if (!params.noprefetch) {
      try {
        block_reader = std::make_unique<PCAone::PgenBlockReader>(params.filein + ".pgen", nsamples, reader_threads);
      } catch (const std::exception& e) {
        cao.warn("reading the PGEN variant by variant:", e.what());
      }
    }
  }
  if (!block_reader) return;
  block_variants.resize(count);
  for (uint64 j = 0; j < count; ++j) block_variants[j] = requests[j].variant;
  block_reader->load(block_variants);
  // read ahead the next block of a pass, or the first block after the last
  const bool sequential = start_idx == 0 || start_idx == last_end;
  last_end = stop_idx + 1;
  if (start_idx == 0) first_count = count;
  if (sequential && count < nsnps) {
    uint64 next = stop_idx + 1, n = std::min<uint64>(count, nsnps - std::min<uint64>(next, nsnps));
    if (next >= nsnps) next = 0, n = std::min<uint64>(first_count ? first_count : count, nsnps);
    variants_of(next, n, next_variants);
    block_reader->prefetch(next_variants);
  }
}

void FilePgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize) G = Mat2D::Zero(nsamples, actual_block_size);

  begin_block(start_idx, stop_idx);
  uint i, j, r;
  uint64 snp_idx;

  if (!params.dopca) frequency_was_estimated = true;
  if (frequency_was_estimated) {
    // centred: the 2-bit calls through the SNP's 4 values, straight into G, then
    // the dosages over them, with the arithmetic of the loop below
    const bool lookup = block_reader && params.center;
#pragma omp parallel for private(i, j, r, snp_idx) schedule(static)
    for (r = 0; r < actual_block_size; ++r) {
      i = requests[r].column;
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      snp_idx = start_idx + i;
      uint64 pgen_idx = requests[r].variant;
      const double f = F(snp_idx);
      double scale_factor = 1.0;
      if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
        double sd = sqrt(f * (1.0 - f));
        if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
      }
      if (lookup) {
        alignas(16) double table[32];
        double* g = &G(0, i);
        if (!dosage_mode) {
          for (int c = 0; c < 4; ++c) table[2 * c] = centered_geno_lookup(c, snp_idx) * scale_factor;
          plink2::InitLookup16x8bx2(table);
          plink2::GenoarrLookup16x8bx2(block_reader->genovec(thr, (uint32_t)pgen_idx), table, nsamples, g);
          continue;
        }
        const uintptr_t* present;
        const uint16_t* dmain;
        uint32_t dosage_ct = 0;
        const uintptr_t* geno = block_reader->genovec_dosages(thr, (uint32_t)pgen_idx, &present, &dmain, &dosage_ct);
        for (int c = 0; c < 4; ++c) table[2 * c] = centered_pgen_value(c == 3 ? PGEN_MISSING : c, f) * scale_factor;
        plink2::InitLookup16x8bx2(table);
        plink2::GenoarrLookup16x8bx2(geno, table, nsamples, g);
        // as Dosage16ToDoubles: the d-th dosage belongs to the d-th set bit
        for (uint32_t w = 0, d = 0; d < dosage_ct; ++w) {
          for (uintptr_t bits = present[w]; bits; bits &= bits - 1) {
            const uint32_t sample = w * plink2::kBitsPerWord + plink2::ctzw(bits);
            g[sample] = centered_pgen_value(S_CAST(double, dmain[d++]) * 0.00006103515625, f) * scale_factor;
          }
        }
        continue;
      }
      read_variant(thr, pgen_idx, buf);

      for (j = 0; j < nsamples; ++j) {
        if (!params.center) {
          G(j, i) = pgen2dosage(buf[j]);
        } else if (dosage_mode) {
          G(j, i) = centered_pgen_value(buf[j], f);
        } else {
          G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        }
        G(j, i) *= scale_factor;
      }
    }
  } else {
    if (start_idx == 0) nmono_seen = 0;  // a pass over the blocks restarted before F was complete
    uint64 nmono = 0;
#pragma omp parallel for private(i, j, r, snp_idx) schedule(static) reduction(+ : nmono)
    for (r = 0; r < actual_block_size; ++r) {
      i = requests[r].column;
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      snp_idx = start_idx + i;
      read_variant(thr, requests[r].variant, buf);

      uint64 c = 0;
      double sum = 0.0;
      for (j = 0; j < nsamples; ++j) {
        double dosage = pgen2dosage(buf[j]);
        if (dosage != BED_MISSING_VALUE) {
          sum += dosage;
          ++c;
        }
      }
      F(snp_idx) = (c > 0) ? sum / c : 0.0;
      if (F(snp_idx) == 0.0 || F(snp_idx) == 1.0) ++nmono;
      if (!dosage_mode) {
        // centered_geno_lookup: rows 0=HomRef, 1=Het, 2=HomAlt, 3=missing
        centered_geno_lookup(3, snp_idx) = 0.0;               // missing: impute to mean
        centered_geno_lookup(0, snp_idx) = 0.0 - F(snp_idx);  // HomRef
        centered_geno_lookup(1, snp_idx) = 0.5 - F(snp_idx);  // Het
        centered_geno_lookup(2, snp_idx) = 1.0 - F(snp_idx);  // HomAlt
      }
      const double f = F(snp_idx);
      double scale_factor = 1.0;
      if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
        double sd = sqrt(f * (1.0 - f));
        if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
      }
      for (j = 0; j < nsamples; ++j) {
        if (dosage_mode) {
          G(j, i) = centered_pgen_value(buf[j], f);
        } else {
          G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        }
        G(j, i) *= scale_factor;
      }
    }
    nmono_seen += nmono;
  }

  if (stop_idx + 1 == nsnps && !frequency_was_estimated) {
    frequency_was_estimated = true;
    warn_monomorphic(nmono_seen);
  }
}

void FilePgen::read_block_update(
    uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals, const Mat2D& VT, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize) G = Mat2D::Zero(nsamples, actual_block_size);

  begin_block(start_idx, stop_idx);
  uint i, j, k, r;
  uint64 snp_idx;
  uint ks = svals.rows();

#pragma omp parallel for private(i, j, k, r, snp_idx) schedule(static)
  for (r = 0; r < actual_block_size; ++r) {
    i = requests[r].column;
    int thr = omp_get_thread_num();
    double* buf = thread_bufs[thr].data();
    snp_idx = start_idx + i;
    read_variant(thr, requests[r].variant, buf);
    const double f = F(snp_idx);
    double scale_factor = 1.0;
    if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
      double sd = sqrt(f * (1.0 - f));
      if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
    }

    for (j = 0; j < nsamples; ++j) {
      bool is_missing = (buf[j] == PGEN_MISSING);
      if (dosage_mode) {
        G(j, i) = centered_pgen_value(buf[j], f);
      } else {
        G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
      }
      if (params.emu && is_missing) {  // missing: predict via EMU
        G(j, i) = 0.0;
        for (k = 0; k < ks; ++k) G(j, i) += U(j, k) * svals(k) * VT(k, snp_idx);
        G(j, i) = fmin(fmax(G(j, i), 0.0 - f), 1.0 - f);
      }
      G(j, i) *= scale_factor;
    }
  }
}
