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
  if (params.dopca && !frequency_was_estimated) {
    F = Mat1D::Zero(nsnps);
#pragma omp parallel for private(i, j) schedule(dynamic)
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

  if (params.missme) C = ArrBool::Zero(nsnps * nsamples);

#pragma omp parallel for private(i, j) schedule(dynamic)
  for (i = 0; i < nsnps; ++i) {
    int thr = omp_get_thread_num();
    double* buf = thread_bufs[thr].data();
    uint s = params.filterSNP ? keepSNPs[i] : i;
    if (dosage_mode) {
      reader.Read(buf, nsamples, thr, s, 1);
    } else {
      reader.ReadHardcalls(buf, nsamples, thr, s, 1);
    }
    for (j = 0; j < nsamples; ++j) {
      G(j, i) = pgen2dosage(buf[j]);
      if ((params.missme && G(j, i) == BED_MISSING_VALUE)) C[i * nsamples + j] = 1;
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

static PermMat compute_pgen_perm_indices(uint nsnps, uint nbands) {
  if (nsnps == 0) cao.error("Cannot permute an empty PGEN file");
  if (nbands == 0) cao.error("Cannot permute PGEN file with zero bands");
  if (nbands > nsnps) nbands = nsnps;
  uint bufsize = nsnps / nbands;
  uint twoGB_snps = bufsize * nbands;
  uint nblocks = (nsnps + twoGB_snps - 1) / twoGB_snps;
  uint modr2 = nsnps % twoGB_snps;
  uint64 bufidx = bufsize;

  uint modr = nsnps % nbands;
  uint bandsize = (nsnps + nbands - 1) / nbands;
  std::vector<uint64> bandidx(nbands);
  if (modr == 0) {
    for (uint i = 0; i < nbands; ++i) bandidx[i] = (uint64)i * bandsize;
  } else {
    for (uint i = 0; i < nbands; ++i) {
      if (i < modr)
        bandidx[i] = (uint64)i * bandsize;
      else
        bandidx[i] = (uint64)modr * bandsize + (uint64)(bandsize - 1) * (i - modr);
    }
  }

  Eigen::VectorXi indices(nsnps);
  for (uint i = 0; i < nblocks; ++i) {
    uint64 cur_bufsize = bufsize;
    uint cur_modr2 = modr2;
    if (i == nblocks - 1 && modr2 != 0) {
      uint64 twoGB_snps2 = nsnps - (uint64)(nblocks - 1) * twoGB_snps;
      cur_bufsize = (twoGB_snps2 + nbands - 1) / nbands;
      cur_modr2 = twoGB_snps2 % nbands;
    }
    for (uint b = 0; b < nbands; ++b) {
      for (uint64 j = 0; j < cur_bufsize - 1; ++j) {
        uint64 ia = (uint64)i * twoGB_snps + j * nbands + b;
        uint64 ib = (uint64)i * bufidx + bandidx[b] + j;
        indices(ib) = (int)ia;
      }
      uint64 j = cur_bufsize - 1;
      if (i != nblocks - 1 || cur_modr2 == 0 || b < cur_modr2) {
        uint64 ia = (uint64)i * twoGB_snps + j * nbands + b;
        uint64 ib = (uint64)i * bufidx + bandidx[b] + j;
        indices(ib) = (int)ia;
      }
    }
  }
  return PermMat(indices);
}

static void write_u32_le(std::ofstream& out, uint32_t v) {
  uchar bytes[4] = {static_cast<uchar>(v & 0xff), static_cast<uchar>((v >> 8) & 0xff),
                    static_cast<uchar>((v >> 16) & 0xff), static_cast<uchar>((v >> 24) & 0xff)};
  out.write(reinterpret_cast<const char*>(bytes), 4);
}

static void write_u16_le(std::ofstream& out, uint16_t v) {
  uchar bytes[2] = {static_cast<uchar>(v & 0xff), static_cast<uchar>((v >> 8) & 0xff)};
  out.write(reinterpret_cast<const char*>(bytes), 2);
}

static void pack_pgen_hardcalls(const std::vector<double>& hardcalls, std::vector<uchar>& packed) {
  std::fill(packed.begin(), packed.end(), 0);
  for (uint s = 0; s < hardcalls.size(); ++s) {
    uint code = (hardcalls[s] == PGEN_MISSING) ? 3u : static_cast<uint>(hardcalls[s]);
    packed[s >> 2] |= static_cast<uchar>((code & 3u) << (2 * (s & 3u)));
  }
}

static uint16_t pgen_dosage16(double v) {
  if (v == PGEN_MISSING) return 65535;
  v = fmin(fmax(v, 0.0), 2.0);
  return static_cast<uint16_t>(llround(v * 16384.0));
}

PermMat permute_plink_pgen(std::string& fin, const std::string& fout, uint nbands) {
  std::string fpsam = fin + ".psam";
  std::string fpgen = fin + ".pgen";
  std::string fpvar = fin + ".pvar";

  uint nsamples = 0;
  {
    std::ifstream psam(fpsam);
    if (!psam.is_open()) cao.error("Cannot open psam file: " + fpsam);
    std::string line;
    while (std::getline(psam, line))
      if (!line.empty() && line[0] != '#') ++nsamples;
  }

  PgenReader reader;
  reader.Load(fpgen, nsamples, {}, 1);
  uint nsnps = reader.GetVariantCt();
  bool dosage_mode = reader.DosagePresent();
  cao.print(tick.date(), "permute pgen files. nsnps:", nsnps, ", nsamples:", nsamples,
            ". dosage_present:", dosage_mode);

  PermMat perm = compute_pgen_perm_indices(nsnps, nbands);

  std::ifstream in_pvar(fpvar);
  if (!in_pvar.is_open()) cao.error("Cannot open pvar file: " + fpvar);
  std::ofstream out_pvar(fout + ".perm.pvar");
  if (!out_pvar.is_open()) cao.error("Cannot write pvar file: " + fout + ".perm.pvar");
  std::vector<std::string> pvar_header;
  std::vector<std::string> pvar_records;
  std::string line;
  while (std::getline(in_pvar, line)) {
    if (!line.empty() && line[0] == '#')
      pvar_header.push_back(line);
    else
      pvar_records.push_back(line);
  }
  if (pvar_records.size() != nsnps) cao.error("pvar variant count does not match pgen variant count");
  for (const auto& header_line : pvar_header) out_pvar << header_line << "\n";
  for (int i = 0; i < perm.indices().size(); ++i) out_pvar << pvar_records[perm.indices()(i)] << "\n";

  std::ifstream in_psam(fpsam);
  std::ofstream out_psam(fout + ".perm.psam");
  if (!out_psam.is_open()) cao.error("Cannot write psam file: " + fout + ".perm.psam");
  out_psam << in_psam.rdbuf();

  std::ofstream out_pgen(fout + ".perm.pgen", std::ios::binary);
  if (!out_pgen.is_open()) cao.error("Cannot write pgen file: " + fout + ".perm.pgen");
  uchar magic[3] = {0x6c, 0x1b, static_cast<uchar>(dosage_mode ? 0x03 : 0x02)};
  out_pgen.write(reinterpret_cast<char*>(&magic[0]), 3);
  write_u32_le(out_pgen, nsnps);
  write_u32_le(out_pgen, nsamples);
  out_pgen.put(0);

  std::vector<double> hardcalls(nsamples);
  std::vector<double> dosages(nsamples);
  std::vector<uchar> packed((nsamples + 3) >> 2);
  for (int i = 0; i < perm.indices().size(); ++i) {
    uint src_idx = static_cast<uint>(perm.indices()(i));
    reader.ReadHardcalls(hardcalls.data(), nsamples, 0, src_idx, 1);
    pack_pgen_hardcalls(hardcalls, packed);
    out_pgen.write(reinterpret_cast<char*>(packed.data()), packed.size());
    if (dosage_mode) {
      reader.Read(dosages.data(), nsamples, 0, src_idx, 1);
      for (uint s = 0; s < nsamples; ++s) write_u16_le(out_pgen, pgen_dosage16(dosages[s]));
    }
  }

  reader.Close();
  fin = fout + ".perm";
  return perm;
}

void FilePgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize) G = Mat2D::Zero(nsamples, actual_block_size);

  uint i, j;
  uint64 snp_idx;

  if (!params.dopca) frequency_was_estimated = true;
  if (frequency_was_estimated) {
#pragma omp parallel for private(i, j, snp_idx) schedule(dynamic)
    for (i = 0; i < actual_block_size; ++i) {
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      snp_idx = start_idx + i;
      if (dosage_mode) {
        reader.Read(buf, nsamples, thr, snp_idx, 1);
      } else {
        reader.ReadHardcalls(buf, nsamples, thr, snp_idx, 1);
      }
      for (j = 0; j < nsamples; ++j) {
        if (dosage_mode) {
          G(j, i) = centered_pgen_value(buf[j], F(snp_idx));
        } else {
          G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        }
        if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
          double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
          if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
        }
      }
    }
  } else {
#pragma omp parallel for private(i, j, snp_idx) schedule(dynamic)
    for (i = 0; i < actual_block_size; ++i) {
      int thr = omp_get_thread_num();
      double* buf = thread_bufs[thr].data();
      snp_idx = start_idx + i;
      if (dosage_mode) {
        reader.Read(buf, nsamples, thr, snp_idx, 1);
      } else {
        reader.ReadHardcalls(buf, nsamples, thr, snp_idx, 1);
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
        centered_geno_lookup(3, snp_idx) = 0.0;               // missing: impute to mean
        centered_geno_lookup(0, snp_idx) = 0.0 - F(snp_idx);  // HomRef
        centered_geno_lookup(1, snp_idx) = 0.5 - F(snp_idx);  // Het
        centered_geno_lookup(2, snp_idx) = 1.0 - F(snp_idx);  // HomAlt
      }
      for (j = 0; j < nsamples; ++j) {
        if (dosage_mode) {
          G(j, i) = centered_pgen_value(buf[j], F(snp_idx));
        } else {
          G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        }
        if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
          double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
          if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
        }
      }
    }
  }

  if (stop_idx + 1 == nsnps) frequency_was_estimated = true;
}

void FilePgen::read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals,
                                 const Mat2D& VT, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize) G = Mat2D::Zero(nsamples, actual_block_size);

  uint i, j, k;
  uint64 snp_idx;
  uint ks = svals.rows();

#pragma omp parallel for private(i, j, k, snp_idx) schedule(dynamic)
  for (i = 0; i < actual_block_size; ++i) {
    int thr = omp_get_thread_num();
    double* buf = thread_bufs[thr].data();
    snp_idx = start_idx + i;
    if (dosage_mode) {
      reader.Read(buf, nsamples, thr, snp_idx, 1);
    } else {
      reader.ReadHardcalls(buf, nsamples, thr, snp_idx, 1);
    }
    for (j = 0; j < nsamples; ++j) {
      bool is_missing = (buf[j] == PGEN_MISSING);
      if (dosage_mode) {
        G(j, i) = centered_pgen_value(buf[j], F(snp_idx));
      } else {
        G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
      }
      if (params.emu && is_missing) {  // missing: predict via EMU
        G(j, i) = 0.0;
        for (k = 0; k < ks; ++k) G(j, i) += U(j, k) * svals(k) * VT(k, snp_idx);
        G(j, i) = fmin(fmax(G(j, i), -F(snp_idx)), 1.0 - F(snp_idx));
      }
      if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
        double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
        if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
      }
    }
  }
}
