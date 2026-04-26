/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FilePgen.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FilePgen.hpp"

#include "Common.hpp"
#include "Utils.hpp"
#include <omp.h>

using namespace std;

// ReadHardcalls(allele_idx=1) returns counts of ALT allele: 0.0=HomRef,
// 1.0=Het, 2.0=HomAlt, -3.0=missing.  Divide by 2 to get [0,1] dosage
// matching the BED convention (BED2GENO = {1,-9,0.5,0} for A1=ALT).
static constexpr double PGEN_MISSING = -3.0;
static inline double pgen2dosage(double v) {
  return (v == PGEN_MISSING) ? BED_MISSING_VALUE : v / 2.0;
}

static inline double centered_pgen_value(double v, double af) {
  double dosage = pgen2dosage(v);
  return (dosage == BED_MISSING_VALUE) ? 0.0 : dosage - af;
}

void FilePgen::read_all() {
  uint i, j;
  if (params.dopca && !frequency_was_estimated) {
    F = Mat1D::Zero(nsnps);
 #pragma omp parallel for private(i, j) schedule(static)
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
      if (F(i) == 0.0 || F(i) == 1.0) cao.warn("sites with MAF=0 found! remove them first! SNP index:", i);
      if (params.ld && params.verbose > 1 && F(i) == 0.5)
        cao.warn("sites with MAF=0.5 found in LD estimation. NaN values expected! SNP index:", i);
    }
    filter_snps_resize_F();
  }
  const bool filter = !keepSNPs.empty();
  if (filter) nsnps = keepSNPs.size();
  G = Mat2D::Zero(nsamples, nsnps);

  if(params.missme) C = ArrBool::Zero(nsnps * nsamples);

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
      if((params.missme && G(j, i)== BED_MISSING_VALUE))
        C[i * nsamples + j] = 1;
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

  auto rng = std::default_random_engine{};
  rng.seed(seed);
  for (auto& source : source_by_thread) std::shuffle(source.begin(), source.end(), rng);

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

static std::string trim_pgen_list_line(const std::string& line) {
  const auto first = line.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = line.find_last_not_of(" \t\r\n");
  return line.substr(first, last - first + 1);
}

static std::vector<std::string> read_pgen_prefix_list(const std::string& path) {
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("Cannot open --pgen-list file: " + path);
  std::vector<std::string> prefixes;
  std::string line;
  while (std::getline(fin, line)) {
    line = trim_pgen_list_line(line);
    if (line.empty() || line[0] == '#') continue;
    prefixes.push_back(line);
  }
  if (prefixes.empty()) cao.error("--pgen-list is empty: " + path);
  return prefixes;
}

static std::string read_text_file(const std::string& path) {
  std::ifstream fin(path, std::ios::binary);
  if (!fin.is_open()) cao.error("Cannot open file: " + path);
  return std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
}

static uint count_psam_samples(const std::string& path) {
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("Cannot open psam file: " + path);
  uint n = 0;
  std::string line;
  while (std::getline(fin, line)) {
    if (!line.empty() && line[0] != '#') ++n;
  }
  return n;
}

static void read_pvar_lines(const std::string& path, std::vector<std::string>* headers,
                            std::vector<std::string>& variants) {
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("Cannot open pvar file: " + path);
  std::string line;
  while (std::getline(fin, line)) {
    if (!line.empty() && line[0] == '#') {
      if (headers) headers->push_back(line);
    } else if (!line.empty()) {
      variants.push_back(line);
    }
  }
}

static void write_u32_le(std::ofstream& out, uint32_t value) {
  unsigned char bytes[4] = {static_cast<unsigned char>(value & 0xff),
                            static_cast<unsigned char>((value >> 8) & 0xff),
                            static_cast<unsigned char>((value >> 16) & 0xff),
                            static_cast<unsigned char>((value >> 24) & 0xff)};
  out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

static unsigned char pgen_hardcall_code(double value) {
  if (value == PGEN_MISSING) return 3;
  if (value <= 0.5) return 0;
  if (value <= 1.5) return 1;
  return 2;
}

void merge_permute_pgen(const Param& params, uint nthreads) {
  struct PgenShard {
    std::string prefix;
    uint nsnps = 0;
    bool dosage_mode = false;
    std::unique_ptr<PgenReader> reader;
  };

  const auto prefixes = read_pgen_prefix_list(params.pgen_list);
  const uint reader_threads = std::max(1u, nthreads);
  const uint nsamples = count_psam_samples(prefixes.front() + ".psam");
  const std::string first_psam = read_text_file(prefixes.front() + ".psam");

  std::vector<PgenShard> shards;
  std::vector<uint64> offsets;
  std::vector<std::string> pvar_headers;
  std::vector<std::string> pvar_lines;
  offsets.push_back(0);

  cao.print(tick.date(), "merge and permute PGEN shards. inputs:", prefixes.size(), ", samples:", nsamples);
  for (size_t i = 0; i < prefixes.size(); ++i) {
    const std::string& prefix = prefixes[i];
    const uint shard_samples = count_psam_samples(prefix + ".psam");
    if (shard_samples != nsamples)
      cao.error("sample count mismatch in " + prefix + ".psam");
    if (i > 0 && read_text_file(prefix + ".psam") != first_psam)
      cao.error("all .psam files must be identical and in the same sample order. mismatch: " + prefix + ".psam");

    std::vector<std::string> shard_pvars;
    read_pvar_lines(prefix + ".pvar", i == 0 ? &pvar_headers : nullptr, shard_pvars);

    auto reader = std::make_unique<PgenReader>();
    reader->Load(prefix + ".pgen", nsamples, {}, reader_threads);
    const uint nsnps = reader->GetVariantCt();
    const bool dosage_mode = reader->DosagePresent();
    if (shard_pvars.size() != nsnps)
      cao.error("variant count mismatch between .pgen and .pvar for prefix: " + prefix);

    pvar_lines.insert(pvar_lines.end(), shard_pvars.begin(), shard_pvars.end());
    shards.push_back(PgenShard{prefix, nsnps, dosage_mode, std::move(reader)});
    offsets.push_back(offsets.back() + nsnps);
    cao.print(tick.date(), "shard", i + 1, ", variants:", nsnps, ", dosage_mode:", dosage_mode, ", prefix:", prefix);
  }

  const uint64 total_snps = offsets.back();
  if (total_snps > 0x7ffffffdULL)
    cao.error("basic PGEN writer currently supports at most 2^31-3 variants.");
  cao.print(tick.date(), "total variants before permutation:", total_snps);

  std::vector<uint32_t> perm(total_snps);
  std::iota(perm.begin(), perm.end(), 0);
  std::default_random_engine rng;
  rng.seed(params.seed);
  std::shuffle(perm.begin(), perm.end(), rng);

  {
    std::ofstream psam_out(params.fileout + ".psam", std::ios::binary);
    if (!psam_out.is_open()) cao.error("Cannot write psam file: " + params.fileout + ".psam");
    psam_out << first_psam;
  }

  std::ofstream pvar_out(params.fileout + ".pvar");
  if (!pvar_out.is_open()) cao.error("Cannot write pvar file: " + params.fileout + ".pvar");
  for (const auto& header : pvar_headers) pvar_out << header << '\n';

  std::ofstream pgen_out(params.fileout + ".pgen", std::ios::binary);
  if (!pgen_out.is_open()) cao.error("Cannot write pgen file: " + params.fileout + ".pgen");
  const unsigned char magic[3] = {0x6c, 0x1b, 0x02};
  pgen_out.write(reinterpret_cast<const char*>(magic), sizeof(magic));
  write_u32_le(pgen_out, static_cast<uint32_t>(total_snps));
  write_u32_le(pgen_out, nsamples);
  pgen_out.put('\0');

  const uint bytes_per_variant = (nsamples + 3) >> 2;
  const uint64 target_block_bytes = std::max<uint64>(1, params.buffer) * 1073741824ULL;
  const uint64 variants_per_block = std::max<uint64>(1, target_block_bytes / bytes_per_variant);
  std::vector<std::vector<double>> thread_genotypes(reader_threads, std::vector<double>(nsamples));

  struct PgenMergeRequest {
    uint32_t src_global = 0;
    uint64 out_offset = 0;
    size_t shard_idx = 0;
    uint local_idx = 0;
  };

  cao.print(tick.date(), "physical PGEN permutation block size: ", variants_per_block,
            " variants (~", (variants_per_block * bytes_per_variant) / 1048576, " MiB)");

  for (uint64 block_start = 0; block_start < total_snps; block_start += variants_per_block) {
    const uint64 block_len = std::min<uint64>(variants_per_block, total_snps - block_start);
    std::vector<PgenMergeRequest> requests(block_len);
    for (uint64 j = 0; j < block_len; ++j) {
      const uint32_t src_global = perm[block_start + j];
      auto shard_it = std::upper_bound(offsets.begin(), offsets.end(), src_global);
      const size_t shard_idx = static_cast<size_t>(std::distance(offsets.begin(), shard_it) - 1);
      requests[j] = PgenMergeRequest{src_global, j, shard_idx, static_cast<uint>(src_global - offsets[shard_idx])};
    }

    std::sort(requests.begin(), requests.end(), [](const PgenMergeRequest& a, const PgenMergeRequest& b) {
      if (a.shard_idx != b.shard_idx) return a.shard_idx < b.shard_idx;
      return a.local_idx < b.local_idx;
    });

    std::vector<unsigned char> out_block(block_len * bytes_per_variant);

#pragma omp parallel for schedule(static) num_threads(reader_threads)
    for (uint64 req_idx = 0; req_idx < block_len; ++req_idx) {
      const auto& req = requests[req_idx];
      const int thr = omp_get_thread_num();
      double* genotypes = thread_genotypes[thr].data();
      if (shards[req.shard_idx].dosage_mode) {
        shards[req.shard_idx].reader->Read(genotypes, nsamples, thr, req.local_idx, 1);
      } else {
        shards[req.shard_idx].reader->ReadHardcalls(genotypes, nsamples, thr, req.local_idx, 1);
      }

      unsigned char* packed = out_block.data() + req.out_offset * bytes_per_variant;
      for (uint sample_idx = 0; sample_idx < nsamples; ++sample_idx) {
        packed[sample_idx >> 2] |= static_cast<unsigned char>(pgen_hardcall_code(genotypes[sample_idx])
                                                              << ((sample_idx & 3) * 2));
      }
    }

    pgen_out.write(reinterpret_cast<const char*>(out_block.data()), out_block.size());
    for (uint64 j = 0; j < block_len; ++j) pvar_out << pvar_lines[perm[block_start + j]] << '\n';

    if (params.verbose > 1)
      cao.print(tick.date(), "  written variants:", block_start + block_len, "/", total_snps);
  }
  cao.print(tick.date(), "wrote permuted PGEN prefix:", params.fileout,
            "(.pgen/.pvar/.psam), variants:", total_snps);
}

void FilePgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize)
    G = Mat2D::Zero(nsamples, actual_block_size);

  uint i, j;
  uint64 snp_idx;

  if (!params.dopca) frequency_was_estimated = true;
  if (frequency_was_estimated) {
 #pragma omp parallel for private(i, j, snp_idx) schedule(static)
    for (i = 0; i < actual_block_size; ++i) {
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      snp_idx = start_idx + i;
      uint64 pgen_idx = params.perm ? (uint64)perm.indices()(snp_idx) : snp_idx;
      if (dosage_mode) {        
        reader.Read(buf, nsamples, thr, pgen_idx, 1);
      } else {
        reader.ReadHardcalls(buf, nsamples, thr, pgen_idx, 1);
      }
      const double f = F(snp_idx);
      double scale_factor = 1.0;
      if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
        double sd = sqrt(f * (1.0 - f));
        if (sd > VAR_TOL) scale_factor = sqrt((double)params.ploidy) / sd;
      }

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
 #pragma omp parallel for private(i, j, snp_idx) schedule(static)
    for (i = 0; i < actual_block_size; ++i) {
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      snp_idx = start_idx + i;
      uint64 pgen_idx = params.perm ? (uint64)perm.indices()(snp_idx) : snp_idx;
      if (dosage_mode) {        
        reader.Read(buf, nsamples, thr, pgen_idx, 1);
      } else {
        reader.ReadHardcalls(buf, nsamples, thr, pgen_idx, 1);
      }

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
      if (F(snp_idx) == 0.0 || F(snp_idx) == 1.0) cao.warn("sites with MAF=0 found! remove them first!");
      if (params.ld && params.verbose > 1 && F(snp_idx) == 0.5)
        cao.warn("MAF for site ", snp_idx, " is 0.5. NaN values expected in calculating LD R2.");
      if (!dosage_mode) {
        // centered_geno_lookup: rows 0=HomRef, 1=Het, 2=HomAlt, 3=missing
        centered_geno_lookup(3, snp_idx) = 0.0;                    // missing: impute to mean
        centered_geno_lookup(0, snp_idx) = 0.0 - F(snp_idx);       // HomRef
        centered_geno_lookup(1, snp_idx) = 0.5 - F(snp_idx);       // Het
        centered_geno_lookup(2, snp_idx) = 1.0 - F(snp_idx);       // HomAlt
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
  }

  if (stop_idx + 1 == nsnps) frequency_was_estimated = true;
}

void FilePgen::read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals,
                                 const Mat2D& VT, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize)
    G = Mat2D::Zero(nsamples, actual_block_size);

  uint i, j, k;
  uint64 snp_idx;
  uint ks = svals.rows();

 #pragma omp parallel for private(i, j, k, snp_idx) schedule(static)
  for (i = 0; i < actual_block_size; ++i) {
    int thr = omp_get_thread_num();
    double* buf = thread_bufs[thr].data();
    snp_idx = start_idx + i;
    uint64 pgen_idx = params.perm ? (uint64)perm.indices()(snp_idx) : snp_idx;
    if (dosage_mode) {        
      reader.Read(buf, nsamples, thr, pgen_idx, 1);
    } else {
      reader.ReadHardcalls(buf, nsamples, thr, pgen_idx, 1);
    }
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
