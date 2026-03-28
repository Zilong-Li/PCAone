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
  const bool filter = keepSNPs.size() > 0 ? true : false;
  if (filter) nsnps = keepSNPs.size();
  G = Mat2D::Zero(nsamples, nsnps);

  if(params.miss) C = ArrBool::Zero(nsnps * nsamples);

 #pragma omp parallel for private(i, j) schedule(static)
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
      if((params.miss && G(j, i)== BED_MISSING_VALUE))
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

  if (params.miss) {
      impute = true;
      p_miss = (double)C.count() / (double)C.size();
      cao.print(tick.date(), "the proportion of missingness  is", p_miss);
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
      for (j = 0; j < nsamples; ++j) {
        if (dosage_mode) {
          G(j, i) = centered_pgen_value(buf[j], F(snp_idx));
        } else {
          G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        }
        if (standardize && params.scale == -9) {
          double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
          if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
        }
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
      for (j = 0; j < nsamples; ++j) {
        if (dosage_mode) {
          G(j, i) = centered_pgen_value(buf[j], F(snp_idx));
        } else {
          G(j, i) = centered_geno_lookup(pgen_code(buf[j]), snp_idx);
        }
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
      if (standardize && params.scale == -9) {
        double sd = sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
        if (sd > VAR_TOL) G(j, i) = G(j, i) * sqrt((double)params.ploidy) / sd;
      }
    }
  }
}
