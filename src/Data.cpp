/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Data.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Data.hpp"

#include "Eigen/src/Core/util/Meta.h"
#include "LD.hpp"
#include "Utils.hpp"

using namespace std;

void Data::prepare() {
  if (nsamples > nsnps) nsamples_ge_nsnps = true;

  if (!params.out_of_core) {
    tick.clock();
    read_all();
    readtime += tick.reltime();
    cao.print(tick.date(), "done reading all data");
    return;
  }

  // some common settings for out-of-core
  if (params.pca) {
    if (params.svd_t == SvdType::IRAM) {
      // ram of arnoldi = n * b * 8 / 1024 kb
      blocksize = (uint)ceil((double)params.memory * 134217728 / nsamples);
    } else {
      // ram of halko = (3*n*l + 2*m*l + 5*m + n*b)*8/1024 Kb
      uint l = params.k + params.oversamples;
      double m = (double)(3 * nsamples * l + 2 * nsnps * l + 5 * nsnps) / 134217728;
      if (params.memory > 1.1 * m)
        m = 0;
      else
        cao.warn("minimum RAM required is ", m, " GB. trying to allocate more RAM.");
      blocksize = (unsigned int)ceil(
          (double)((m + params.memory) * 134217728 - 3 * nsamples * l - 2 * nsnps * l - 5 * nsnps) /
          nsamples);
    }
  } else {
    // ram of non-pca run
    blocksize = (uint)ceil((double)params.memory * 134217728 / nsamples);
  }

  nblocks = (unsigned int)ceil((double)nsnps / blocksize);
  cao.print(tick.date(), "initial setting by -m/--memory: blocksize =", blocksize, ", nblocks =", nblocks,
            ", factor =", bandFactor);
  if (nblocks == 1) cao.error("only one block exists. please remove -m option");
  if (params.pca && params.svd_t == SvdType::PCAoneAlg2) {
    // decrease blocksize to fit the fancy halko
    if (nblocks < params.bands) {
      blocksize = (unsigned int)ceil((double)nsnps / params.bands);
    } else {
      bandFactor = (unsigned int)ceil((double)nblocks / params.bands);
      blocksize = (unsigned int)ceil((double)nsnps / (params.bands * bandFactor));
    }
    nblocks = (unsigned int)ceil((double)nsnps / blocksize);
    cao.print(tick.date(), "after adjustment by PCAone: -w =", params.bands, ", blocksize =", blocksize,
              ", nblocks =", nblocks, ", factor =", bandFactor);
  }
  start.resize(nblocks);
  stop.resize(nblocks);
  for (uint i = 0; i < nblocks; i++) {
    start[i] = i * blocksize;
    stop[i] = start[i] + blocksize - 1;
    stop[i] = stop[i] >= nsnps ? nsnps - 1 : stop[i];
  }
}

void Data::filter_snps_resize_F() {
  if (params.keepsnp && params.maf > 0 &&
      params.maf <= 0.5) {  // filter snps, update keepSNPs, reassign nsnps;
    Mat1D Fnew(F.size());   // make a temp F
    int i, j;
    for (i = 0, j = 0; j < (int)F.size(); j++) {
      if (MAF(F(j)) > params.maf) {
        keepSNPs.push_back(j);  // keep track of index of element > maf
        Fnew(i++) = F(j);
      }
    }
    nsnps = keepSNPs.size();  // new number of SNPs
    cao.print(tick.date(), "number of SNPs after filtering by MAF >", params.maf, ":", nsnps);
    if (nsnps < 1) cao.error("no SNPs left after filtering!");
    // resize F
    F.noalias() = Fnew.head(nsnps);
  }
}

// initially only works with plink inputs
// but can work with beagle file as long as there is beagle.gz.bim file
void Data::save_snps_in_bim() {
  cao.print(tick.date(), "save matched sites in .mbim file and params.perm is", params.perm);
  // could be permuted
  std::ifstream ifs_bim(params.filein + ".bim");
  if (!ifs_bim.is_open()) cao.error(params.filein + ".bim not exists!");
  std::ofstream ofs_bim(params.fileout + ".mbim");
  std::string line;
  int i, j;
  if (params.perm && params.out_of_core) {
    vector<std::string> bims2;
    bims2.resize(nsnps);
    Mat1D Ftmp(F.size());
    j = 0;
    while (getline(ifs_bim, line)) {
      i = perm.indices()[j];
      bims2[i] = line;
      Ftmp(i) = F(j);
      j++;
    }
    for (j = 0; j < (int)bims2.size(); j++) {
      ofs_bim << bims2[j] << "\t" << Ftmp(j) << "\n";
    }
  } else {  // plink.bim is not permuted
    j = 0, i = 0;
    while (getline(ifs_bim, line)) {
      if (params.keepsnp) {
        if (keepSNPs[i] == j) {
          ofs_bim << line << "\t" << F(i) << "\n";
          i++;
        }
      } else {
        ofs_bim << line << "\t" << F(j) << "\n";
      }
      j++;
    }
  }
  ofs_bim.close();
}

/**
  T = U'/s
  VT = (U'/s) * G = T * G
  V = G' * (U/s) // calculate V is not a good idea
 **/
void Data::calcu_vt_initial(const Mat2D &T, Mat2D &VT, bool standardize) {
  if (nblocks == 1) {
    cao.error("only one block exists. please use in-memory mode instead");
  }
  uint actual_block_size;
  check_file_offset_first_var();
  for (uint i = 0; i < nblocks; ++i) {
    actual_block_size = stop[i] - start[i] + 1;
    // G (nsamples, actual_block_size)
    read_block_initial(start[i], stop[i], standardize);
    VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
  }

  return;
}

void Data::calcu_vt_update(const Mat2D &T, const Mat2D &U, const Mat1D &svals, Mat2D &VT, bool standardize) {
  if (nblocks == 1) {
    cao.error("only one block exists. please use in-memory mode instead");
  }
  uint actual_block_size;
  check_file_offset_first_var();
  for (uint i = 0; i < nblocks; ++i) {
    actual_block_size = stop[i] - start[i] + 1;
    // G (nsamples, actual_block_size)
    read_block_update(start[i], stop[i], U, svals, VT, standardize);
    VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
  }

  return;
}

// S: signular values
// E: eigen values
void Data::write_eigs_files(const Mat1D &E, const Mat1D &S, const Mat2D &U, const Mat2D &V) {
  std::ofstream outs(params.fileout + ".sigvals");
  std::ofstream oute(params.fileout + ".eigvals");
  std::ofstream outu(params.fileout + ".eigvecs");
  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  if (outs.is_open()) {
    outs << '#' << U.rows() << ',' << V.rows() << '\n';
    outs << S.format(fmt) << '\n';
  }
  if (oute.is_open()) oute << (E * params.ploidy).format(fmt) << '\n';
  if (outu.is_open()) outu << U.format(fmt) << '\n';
  if (params.project == 0 && params.printv) {
    save_snps_in_bim();
    std::ofstream outv(params.fileout + ".loadings");
    if (outv.is_open()) outv << V.format(fmt) << '\n';
  }
  cao.print(tick.date(), "eigen vectors and values saved");
}

void Data::write_residuals(const Mat1D &S, const Mat2D &U, const Mat2D &VT) {
  // we always filter snps for in-core mode
  if (params.ld_stats == 1) {
    cao.print(tick.date(), "ld-stats=1: calculate standardized genotype matrix!");
  } else {
    cao.print(tick.date(), "ld-stats=0: calculate the ancestry adjusted LD matrix!");
  }
  std::setlocale(LC_ALL, "C");
  std::ios_base::sync_with_stdio(false);
  std::ofstream ofs(params.fileout + ".residuals", std::ios::binary);
  const uint64 ibyte = 4;
  const uint64 magic = ibyte * 2;
  uint64 bytes_per_snp = nsamples * ibyte;
  ofs.write((char *)&nsnps, ibyte);
  ofs.write((char *)&nsamples, ibyte);
  Eigen::VectorXf fg;
  uint64 idx;
  if (!params.out_of_core) {
    if (params.ld_stats == 0) G -= U * S.asDiagonal() * VT;  // get residuals matrix
    G.rowwise() -= G.colwise().mean();                       // Centering
    for (Eigen::Index ib = 0; ib < G.cols(); ib++) {
      fg = G.col(ib).cast<float>();
      if (params.perm) {
        idx = magic + (uint64)perm.indices()[ib] * bytes_per_snp;
        ofs.seekp(idx, std::ios_base::beg);
      }
      ofs.write((char *)fg.data(), bytes_per_snp);
    }
  } else {
    check_file_offset_first_var();
    int i = 0;
    for (uint b = 0; b < nblocks; ++b) {
      read_block_initial(start[b], stop[b], false);
      // G (nsamples, actual_block_size)
      if (params.ld_stats == 0) G -= U * S.asDiagonal() * VT.middleCols(start[b], G.cols());
      G.rowwise() -= G.colwise().mean();  // Centering
      for (Eigen::Index ib = 0; ib <= stop[b] - start[b]; ib++, i++) {
        fg = G.col(ib).cast<float>();
        if (params.perm) {
          idx = magic + (uint64)perm.indices()[i] * bytes_per_snp;
          ofs.seekp(idx, std::ios_base::beg);
        }
        ofs.write((char *)fg.data(), bytes_per_snp);
      }
    }
  }

  save_snps_in_bim();
  cao.print(tick.date(), "the LD matrix and SNPs info are saved");
}

void Data::update_batch_E(const Mat2D &U, const Mat1D &svals, const Mat2D &VT) {
  uint ks = svals.size();
  if (params.pcangsd) {
// for gp
#pragma omp parallel for
    for (uint j = 0; j < nsnps; ++j) {
      double p0, p1, p2;
      uint s = params.keepsnp ? keepSNPs[j] : j;
      for (uint i = 0; i < nsamples; ++i) {
        // Rescale individual allele frequencies
        double pt = 0.0;
        for (uint k = 0; k < ks; ++k) {
          pt += U(i, k) * svals(k) * VT(k, j);
        }
        pt = (pt + 2.0 * F(j)) / 2.0;
        pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
        // update E, which is G here
        p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
      }
    }
  } else {
// for gt
#pragma omp parallel for
    for (uint i = 0; i < nsnps; ++i) {
      for (uint j = 0; j < nsamples; ++j) {
        if (C[i * nsamples + j])  // no bool & 1
        {                         // sites need to be predicted
          G(j, i) = 0.0;
          for (uint k = 0; k < ks; ++k) {
            G(j, i) += U(j, k) * svals(k) * VT(k, i);
          }
          G(j, i) = fmin(fmax(G(j, i), -F(i)), 1 - F(i));
        }
      }
    }
  }
}

void Data::standardize_E() {
#pragma omp parallel for
  for (uint i = 0; i < nsnps; ++i) {
    for (uint j = 0; j < nsamples; ++j) {
      double sd = sqrt(F(i) * (1 - F(i)));
      // in case denominator is too small.
      if (sd > VAR_TOL) G(j, i) /= sd;
    }
  }
}

void Data::pcangsd_standardize_E(const Mat2D &U, const Mat1D &svals, const Mat2D &VT) {
  cao.print(tick.date(), "begin to standardize the matrix for pcangsd");
  uint nk = svals.size();
  Dc = Mat1D::Zero(nsamples);
#pragma omp parallel
  {
    Mat1D diag_private = Mat1D::Zero(nsamples);  // Thread private vector;
#pragma omp for
    for (uint j = 0; j < nsnps; j++) {
      double p0, p1, p2, pt, pSum, tmp;
      double norm = sqrt(2.0 * F(j) * (1.0 - F(j)));
      uint s = params.keepsnp ? keepSNPs[j] : j;
      for (uint i = 0; i < nsamples; i++) {
        // Rescale individual allele frequencies
        pt = 0.0;
        for (uint k = 0; k < nk; ++k) {
          pt += U(i, k) * svals(k) * VT(k, j);
        }
        pt = (pt + 2.0 * F(j)) / 2.0;
        pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
        // Update e
        p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        pSum = p0 + p1 + p2;
        G(i, j) = (p1 + 2 * p2) / pSum - 2.0 * F(j);
        if (norm > VAR_TOL) G(i, j) /= norm;

        // Update diag
        tmp = (0.0 - 2.0 * F(j)) * (0.0 - 2.0 * F(j)) * (p0 / pSum);
        tmp += (1.0 - 2.0 * F(j)) * (1.0 - 2.0 * F(j)) * (p1 / pSum);
        tmp += (2.0 - 2.0 * F(j)) * (2.0 - 2.0 * F(j)) * (p2 / pSum);
        diag_private[i] += tmp / (2.0 * F(j) * (1.0 - F(j)));
      }
    }
#pragma omp critical
    {
      for (uint i = 0; i < nsamples; i++) {
        Dc[i] += diag_private[i];  // Sum arrays for threads
      }
    }
  }
}

void Data::predict_missing_E(const Mat2D &U, uint64 start_idx, uint64 stop_idx) {
  // assume full size G is in memory
  // Eigen::ColPivHouseholderQR<Mat2D> qr(U);  // nsamples x nk
  const uint nk = U.cols();
// get beta by linear regresion on data without NAs
// Ux = g, x is the beta
#pragma omp parallel for
  for (uint j = start_idx; j <= stop_idx; j++) {
    // find non-missing samples
    Int1D idx, na_idx;
    for (uint i = 0; i < nsamples; i++) {
      if (!C(j * nsamples + i))
        idx.push_back(i);
      else
        na_idx.push_back(i);
    }
    if (idx.size() == 0) continue;
    Mat1D v = U(idx, Eigen::all).colPivHouseholderQr().solve(G(idx, j));
    // predict for NAs in E
    for (auto i : na_idx) {  // samples need to be predicted
      G(i, j) = 0.0;
      for (uint k = 0; k < nk; ++k) {
        G(i, j) += U(i, k) * v(k);
      }
      G(i, j) = fmin(fmax(G(i, j), -F(j)), 1 - F(j));
    }
  }
}
