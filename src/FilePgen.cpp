/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FilePgen.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FilePgen.hpp"

#include "Common.hpp"
#include "Utils.hpp"

using namespace std;

// ReadHardcalls(allele_idx=1) returns counts of ALT allele: 0.0=HomRef,
// 1.0=Het, 2.0=HomAlt, -3.0=missing.  Divide by 2 to get [0,1] dosage
// matching the BED convention (BED2GENO = {1,-9,0.5,0} for A1=ALT).
static constexpr double PGEN_MISSING = -3.0;
static inline double pgen2dosage(double v) {
  return (v == PGEN_MISSING) ? BED_MISSING_VALUE : v / 2.0;
}

void FilePgen::read_all() {
  uint i, j;
  if (params.estaf) {
    F = Mat1D::Zero(nsnps);
    C = ArrBool::Zero(nsnps * nsamples);
    for (i = 0; i < nsnps; ++i) {
      if (dosage_mode) {        
        reader.Read(buf.data(), nsamples, 0, i, 1);
      } else {
        reader.ReadHardcalls(buf.data(), nsamples, 0, i, 1);
      }
      uint64 c = 0;
      double sum = 0.0;
      for (j = 0; j < nsamples; ++j) {
        if (buf[j] != PGEN_MISSING) {
          sum += buf[j] / 2.0;
          ++c;
        } else {
          C[i * nsamples + j] = 1;
        }
      }
      F(i) = (c > 0) ? sum / c : 0.0;
      if (F(i) == 0.0 || F(i) == 1.0) cao.warn("sites with MAF=0 found! remove them first! SNP index:", i);
      if (params.ld && params.verbose > 1 && F(i) == 0.5)
        cao.warn("sites with MAF=0.5 found in LD estimation. NaN values expected! SNP index:", i);
    }
    if (C.count() > 0) {
      impute = true;
      double p_miss = (double)C.count() / (double)C.size();
      cao.print(tick.date(), "the proportion of missingness  is", p_miss);
    } else {
      C.resize(0); // clear 
      filter_snps_resize_F();
    }
  }

  G = Mat2D::Zero(nsamples, nsnps);

  for (i = 0; i < nsnps; ++i) {
    uint s = params.filterSNP ? keepSNPs[i] : i;
    if (dosage_mode) {        
      reader.Read(buf.data(), nsamples, 0, s, 1);
    } else {
      reader.ReadHardcalls(buf.data(), nsamples, 0, s, 1);
    }
    for (j = 0; j < nsamples; ++j) {
      G(j, i) = pgen2dosage(buf[j]);
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

}

PermMat compute_pgen_perm(uint nsnps, uint nbands) {
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

void FilePgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || actual_block_size < blocksize)
    G = Mat2D::Zero(nsamples, actual_block_size);

  uint i, j;
  uint64 snp_idx;

  if (!params.estaf) frequency_was_estimated = true;
  if (frequency_was_estimated) {
    for (i = 0; i < actual_block_size; ++i) {
      snp_idx = start_idx + i;
      uint64 pgen_idx = params.perm ? (uint64)perm.indices()(snp_idx) : snp_idx;
      if (dosage_mode) {        
        reader.Read(buf.data(), nsamples, 0, pgen_idx, 1);
      } else {
        reader.ReadHardcalls(buf.data(), nsamples, 0, pgen_idx, 1);
      }

#pragma omp parallel for private(j)
      for (j = 0; j < nsamples; ++j) {
        G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        if (standardize && params.scale == -9) {
          double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
          if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
        }
      }
    }
  } else {
    for (i = 0; i < actual_block_size; ++i) {
      snp_idx = start_idx + i;
      uint64 pgen_idx = params.perm ? (uint64)perm.indices()(snp_idx) : snp_idx;
      if (dosage_mode) {        
        reader.Read(buf.data(), nsamples, 0, pgen_idx, 1);
      } else {
        reader.ReadHardcalls(buf.data(), nsamples, 0, pgen_idx, 1);
      }

      uint64 c = 0;
      double sum = 0.0;
      for (j = 0; j < nsamples; ++j) {
        if (buf[j] != PGEN_MISSING) {
          sum += buf[j] / 2.0;
          ++c;
        }
      }
      F(snp_idx) = (c > 0) ? sum / c : 0.0;
      if (F(snp_idx) == 0.0 || F(snp_idx) == 1.0) cao.warn("sites with MAF=0 found! remove them first!");
      if (params.ld && params.verbose > 1 && F(snp_idx) == 0.5)
        cao.warn("MAF for site ", snp_idx, " is 0.5. NaN values expected in calculating LD R2.");
      // centered_geno_lookup: rows 0=HomRef, 1=Het, 2=HomAlt, 3=missing
      centered_geno_lookup(3, snp_idx) = 0.0;                    // missing: impute to mean
      centered_geno_lookup(0, snp_idx) = 0.0 - F(snp_idx);       // HomRef
      centered_geno_lookup(1, snp_idx) = 0.5 - F(snp_idx);       // Het
      centered_geno_lookup(2, snp_idx) = 1.0 - F(snp_idx);       // HomAlt
#pragma omp parallel for private(j)
      for (j = 0; j < nsamples; ++j) {
        G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        if (standardize && params.scale == -9) {
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
  if (G.cols() < blocksize || actual_block_size < blocksize)
    G = Mat2D::Zero(nsamples, actual_block_size);

  uint i, j, k;
  uint64 snp_idx;
  uint ks = svals.rows();

  for (i = 0; i < actual_block_size; ++i) {
    snp_idx = start_idx + i;
    uint64 pgen_idx = params.perm ? (uint64)perm.indices()(snp_idx) : snp_idx;
    if (dosage_mode) {        
      reader.Read(buf.data(), nsamples, 0, pgen_idx, 1);
    } else {
      reader.ReadHardcalls(buf.data(), nsamples, 0, pgen_idx, 1);
    }

#pragma omp parallel for private(j, k)
    for (j = 0; j < nsamples; ++j) {
      int code = pgen_code(buf[j]);
      G(j, i) = centered_geno_lookup(code, snp_idx);
      if (params.emu && code == 3) {  // missing: predict via EMU
        G(j, i) = 0.0;
        for (k = 0; k < ks; ++k) G(j, i) += U(j, k) * svals(k) * VT(k, snp_idx);
        G(j, i) = fmin(fmax(G(j, i), -F(snp_idx)), 1.0 - F(snp_idx));
      }
      if (standardize && params.scale == -9) {
        double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
        if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
      }
    }
  }
}
