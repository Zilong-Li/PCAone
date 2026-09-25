/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Data.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Data.hpp"

#include <omp.h>

#include "Cmd.hpp"
#include "Utils.hpp"

using namespace std;

void Data::prepare() {
  if (nsamples > nsnps) nsamples_ge_nsnps = true;

  if (!params.dopca) {  // for projection, read F from reference set
    cao.print(tick.date(), "read allele frequency from .mbim file: " + params.filebim);
    F = read_frq(params.filebim);
    if (!keepRefSNPs.empty()) {
      Mat1D Fnew(keepRefSNPs.size());
      for (int i = 0; i < (int)keepRefSNPs.size(); ++i) {
        Fnew(i) = F(keepRefSNPs[i]);
      }
      F = Fnew;
      cao.print(tick.date(), "keep AF of", F.size(), " matched sites in the reference");
    }
    if (!flipSNPs.empty()) {
      for (int i : flipSNPs) F(i) = 1.0 - F(i);
    }
  }

  if (!params.out_of_core) {
    tick.clock();
    read_all();
    readtime += tick.reltime();
    cao.print(tick.date(), "done reading all data");
    return;
  }

  // some common settings for out-of-core. LD keeps dopca on to estimate F, but
  // runs no PCA, so its blocks are sized like any other non-PCA run
  const bool pca_blocks = params.dopca && !params.ld;
  const bool exact = pca_blocks && params.svd_t == SvdType::FULL;
  if (exact) {
    // exact PCA (Exact.cpp): the N x N GRM and F, plus the loadings for -V,
    // are held besides one block. Without -m, blocks of ~64 MB (and at least
    // 2048 sites) keep the GRM update a large matrix product.
    const double n = std::max(1u, nsamples);
    const double fixed = n * n + 6.0 * nsnps + (params.printv ? (double)nsnps * params.k : 0.0);
    double b = std::max(2048.0, 8388608.0 / n);
    if (params.memory > 0) {
      if (params.memory * 134217728 > 1.1 * fixed + 64 * n)
        b = (params.memory * 134217728 - fixed) / n;
      else
        cao.warn("the exact PCA needs at least", fixed / 134217728, " GB for the", nsamples, " x", nsamples,
                 " GRM. -m is too small, using blocks of", (uint64)std::ceil(b), " sites");
    }
    blocksize = (uint)std::min<double>(std::max(1.0, std::ceil(b)), std::max(1u, nsnps));
  } else if (pca_blocks) {
    if (params.svd_t == SvdType::IRAM) {
      // ram of arnoldi = n * b * 8 / 1024 kb
      blocksize = (uint)ceil((double)params.memory * 134217728 / nsamples);
    } else {
      // ram of halko = (3*n*l + 2*m*l + 5*m + n*b)*8/1024 Kb
      // in doubles: the unsigned products wrapped once nsnps * l reached 2^31
      const double l = (double)params.k + params.oversamples;
      const double fixed = 3.0 * nsamples * l + 2.0 * nsnps * l + 5.0 * nsnps;  // doubles held besides the block
      double m = fixed / 134217728;
      if (params.memory > 1.1 * m)
        m = 0;
      else
        cao.warn("minimum RAM required is ", m, " GB. trying to allocate more RAM.");
      blocksize = (unsigned int)ceil(((m + params.memory) * 134217728 - fixed) / nsamples);
    }
  } else {
    // ram of non-pca run
    blocksize = (uint)ceil((double)params.memory * 134217728 / nsamples);
  }

  nblocks = (unsigned int)ceil((double)nsnps / blocksize);
  if (exact) {
    cao.print(tick.date(), "blocks for the exact PCA: blocksize =", blocksize, ", nblocks =", nblocks);
  } else {
    cao.print(tick.date(), "initial setting by -m/--memory: blocksize =", blocksize, ", nblocks =", nblocks,
              ", factor =", bandFactor);
    if (nblocks == 1) cao.error("only one block exists. please remove -m option");
  }
  if (pca_blocks && params.svd_t == SvdType::PCAoneAlg2) {
    // decrease blocksize for the winSVD
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

// filter snps, update keepSNPs, reassign nsnps;
void Data::filter_snps_resize_F() {
  if (!params.filterSNP) return;
  if (!(params.maf > 0 && params.maf <= 0.5)) cao.error("--maf has to be between (0, 0.5)");

  nsnps_all = F.size();
  Mat1D Fnew(F.size());  // make a temp F
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

// initially only works with plink inputs
// but can work with beagle file as long as there is beagle.gz.bim file
// TODO: always output mbim even though there is no beagle.gz.bim file
void Data::save_snps_in_mbim() {
  const bool has_metadata_header = params.file_t == FileType::PGEN;
  const std::string bim_path = params.filein + (has_metadata_header ? ".pvar" : ".bim");
  std::ifstream ifs_bim(bim_path);
  if (!ifs_bim.is_open()) {
    cao.warn(params.filein + ".bim/.pvar not found; skipping mbim output");
    return;
  }
  std::ofstream ofs_bim(params.fileout + ".mbim");
  std::string line;
  auto read_metadata_line = [&]() {
    while (getline(ifs_bim, line)) {
      if (has_metadata_header && !line.empty() && line[0] == '#') continue;
      if (has_metadata_header) line = pvar_line_to_bim_line(line, bim_path);
      return true;
    }
    return false;
  };
  const bool metadata_is_permuted =
      params.file_t == FileType::PLINK && params.perm && params.out_of_core && perm.indices().size() == nsnps;
  const bool frequency_is_permuted = params.perm && params.out_of_core && perm.indices().size() == nsnps;

  if (!params.filterSNP && !params.perm) {
    for (Eigen::Index j = 0; j < F.size() && read_metadata_line(); ++j) {
      ofs_bim << line << "\t" << F(j) << "\n";
    }
    ofs_bim.close();
    return;
  }

  std::vector<std::string> metadata(nsnps);
  if (metadata_is_permuted) {
    for (Eigen::Index permuted_idx = 0; permuted_idx < nsnps && read_metadata_line(); ++permuted_idx) {
      Eigen::Index original_idx = perm.indices()[permuted_idx];
      metadata[original_idx] = line;
    }
  } else if (params.filterSNP) {
    int kept = 0;
    for (int source_idx = 0; kept < (int)nsnps && read_metadata_line(); ++source_idx) {
      if (kept < (int)keepSNPs.size() && keepSNPs[kept] == source_idx) {
        metadata[kept++] = line;
      }
    }
  } else {
    for (Eigen::Index j = 0; j < nsnps && read_metadata_line(); ++j) {
      metadata[j] = line;
    }
  }

  std::vector<Eigen::Index> original_to_logical(nsnps);
  for (Eigen::Index logical_idx = 0; logical_idx < nsnps; ++logical_idx) {
    Eigen::Index original_idx = frequency_is_permuted ? perm.indices()[logical_idx] : logical_idx;
    original_to_logical[original_idx] = logical_idx;
  }

  for (Eigen::Index original_idx = 0; original_idx < nsnps; ++original_idx) {
    ofs_bim << metadata[original_idx] << "\t" << F(original_to_logical[original_idx]) << "\n";
  }
  ofs_bim.close();
  cao.print(tick.date(), "save matched sites in .mbim file and permutation mode is", params.perm);
}

/**
  T = U'/s
  VT = (U'/s) * G = T * G
  V = G' * (U/s) // calculate V is not a good idea
 **/
void Data::calcu_vt_initial(const Mat2D& T, Mat2D& VT, bool standardize) {
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

void Data::calcu_vt_update(const Mat2D& T, const Mat2D& U, const Mat1D& svals, Mat2D& VT, bool standardize) {
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
void Data::set_svd_transform(bool standardized) {
  // standardize_E() is a no-op unless params.scale is the genetic default, so a
  // decomposition that "standardized" under any other --scale did not actually
  // get the sqrt(ploidy)/sd transform that -P knows how to undo.
  svd_scale = (standardized && params.scale == SCALE_STANDARDIZE_GENETIC) ? SCALE_STANDARDIZE_GENETIC : params.scale;
  if (!standardized && params.scale == SCALE_STANDARDIZE_GENETIC) svd_scale = 0;  // centred only
  // fit_with_pi() builds the pcangsd matrix as dosage - 2f, on the 0..2 scale,
  // rather than PCAone's usual 0..1 coding.
  svd_gscale = params.pcangsd ? 2 : 1;
}

void Data::write_eigs_files(const Mat1D& E, const Mat1D& S, const Mat2D& U, const Mat2D& V) {
  std::ofstream outs(params.fileout + ".sigvals");
  std::ofstream oute(params.fileout + ".eigvals");
  std::ofstream outu(params.fileout + ".eigvecs");
  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  if (outs.is_open()) {
    // key=value fields are appended after nsamples,nsnps; older PCAone parses
    // the first two with stoi and ignores the rest. See UsvTransform.
    outs << '#' << U.rows() << ',' << V.rows() << ",scale=" << svd_scale << ",ploidy=" << params.ploidy
         << ",gscale=" << svd_gscale << '\n';
    outs << S.format(fmt) << '\n';
  }
  if (oute.is_open()) oute << E.format(fmt) << '\n';
  if (outu.is_open()) outu << U.format(fmt) << '\n';
  if (params.printv) {
    save_snps_in_mbim();
    std::ofstream outv(params.fileout + ".loadings");
    if (!outv.is_open()) cao.error("can not open " + params.fileout + ".loadings");
    if (params.perm && V.rows() == nsnps && perm.indices().size() == nsnps) {
      std::vector<Eigen::Index> original_to_logical(nsnps);
      for (Eigen::Index logical_idx = 0; logical_idx < nsnps; ++logical_idx) {
        original_to_logical[perm.indices()[logical_idx]] = logical_idx;
      }
      for (Eigen::Index original_idx = 0; original_idx < V.rows(); ++original_idx) {
        outv << V.row(original_to_logical[original_idx]).format(fmt) << '\n';
      }
    } else {
      outv << V.format(fmt) << '\n';
    }
  }

  cao.print(tick.date(), "eigen vectors and values saved");
}

void Data::fit_with_pi(const Mat2D& U, const Mat1D& svals, const Mat2D& VT) {
  if (params.verbose >= 3) cao.print(tick.date(), "call fit_with_pi");
  uint ks = svals.size();
  if (params.pcangsd) {  // for pcangsd with beagle input
    const Mat2D US = U * svals.asDiagonal();
    for_each_product_column(US, VT, 0, nsnps, [&](Eigen::Index jj, const auto& recon) {
      const uint j = (uint)jj;
      double p0, p1, p2;
      const uint original = unpermuted_snp_index(j);
      const double f = F(original);
      uint s = params.filterSNP ? keepSNPs[original] : original;
      for (uint i = 0; i < nsamples; ++i) {
        // Rescale individual allele frequencies
        double pt = recon(i);
        pt = (pt + 2.0 * f) / 2.0;
        pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
        // update E, which is G here
        p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        G(i, j) = (p1 + 2.0 * p2) / (p0 + p1 + p2) - 2.0 * f;
      }
    });
  }

  if (params.emu) {
#pragma omp parallel for
    for (uint i = 0; i < nsnps; ++i) {
      const uint original = unpermuted_snp_index(i);
      const double f = F(original);
      for (uint j = 0; j < nsamples; ++j) {
        if (C[(uint64)original * nsamples + j]) {  // sites need to be predicted
          G(j, i) = 0.0;
          for (uint k = 0; k < ks; ++k) {
            G(j, i) += U(j, k) * svals(k) * VT(k, i);
          }
          G(j, i) = fmin(fmax(G(j, i), -f), 1 - f);
        }
      }
    }
  }
}

void Data::standardize_E() {
  if (params.verbose >= 3) cao.print(tick.date(), "standardize the matrix");
  if (params.scale != -9) return;
#pragma omp parallel for
  for (uint i = 0; i < nsnps; ++i) {
    const double f = F(unpermuted_snp_index(i));
    const double sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) {  // in case denominator is too small.
      G.col(i) *= sqrt((double)params.ploidy) / sd;
    }
  }
}

// Decide how to scale a target genotype matrix so it matches the matrix a
// reference PCA decomposed.
//
// Before the transform was recorded in .sigvals, projection and selection
// scaled their G from this run's --scale, which silently mismatched whenever
// the reference had used anything else. That included the case where the user
// passed no --scale at all: a reference built with -D/--ld has standardisation
// turned off, and -D is a documented way to get the .mbim that -P needs.
bool Data::resolve_ref_scaling(const UsvTransform& t, const std::string& src, bool allow_dosage) const {
  if (!t.known) {
    cao.warn(src,
             " predates the recording of the PCA scaling, so the reference is assumed to have used this run's "
             "-C/--scale. that is wrong if the reference was built with -D/--ld, or with a non-default --scale; "
             "rerun it to record the transform.");
    return params.scale == SCALE_STANDARDIZE_GENETIC;
  }
  if (t.scale > 0)
    cao.error("the reference PCA used -C/--scale ", t.scale,
              ", which cannot be replayed on the target genotypes here");
  if (t.gscale != 1 && !(allow_dosage && t.gscale == 2))
    cao.error("the reference PCA decomposed dosages on the 0..2 scale (gscale=", t.gscale,
              "), which is not the 0..1 coding used here. this combination was silently mis-scaled before; "
              "rerun the reference on the same genotype coding");
  // mirrors FileUSV::check_transform(): standardize_E() and
  // pcangsd_standardize_E() disagree on what this matrix is, so it has no inverse
  if (t.gscale == 2 && t.scale == SCALE_STANDARDIZE_GENETIC)
    cao.error("the reference PCA standardised a dosage-scale matrix (gscale=2, scale=-9), which cannot be "
              "replayed on the target here. rerun the reference PCA so it decomposes centred dosages");
  if (t.ploidy != params.ploidy)
    cao.error("the reference PCA used ploidy ", t.ploidy, " but this run uses ", params.ploidy,
              "; pass --haploid consistently");
  const bool standardize = (t.scale == SCALE_STANDARDIZE_GENETIC);
  cao.print(tick.date(), "scaling the target genotypes as the reference did:",
            standardize ? "standardized" : (t.gscale == 2 ? "centred dosages (0..2)" : "centred only"),
            "(scale =", t.scale, ", ploidy =", t.ploidy, ", gscale =", t.gscale, ")");
  if (params.scale != t.scale)
    cao.warn("-C/--scale ", params.scale, " is ignored here; the scaling is taken from the reference PCA");
  return standardize;
}

// standardize_E() but driven by the reference's transform rather than this
// run's --scale, and using the reference's ploidy.
void Data::standardize_E_ref(const UsvTransform& t) {
  const double rploidy = sqrt((double)t.ploidy);
#pragma omp parallel for
  for (uint i = 0; i < nsnps; ++i) {
    const double f = F(i);
    const double sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) G.col(i) *= rploidy / sd;
  }
}

// the same, for one out-of-core block already read unstandardized
void Data::standardize_block_ref(const UsvTransform& t, uint64 start_idx, uint block_cols) {
  const double rploidy = sqrt((double)t.ploidy);
#pragma omp parallel for
  for (uint i = 0; i < block_cols; ++i) {
    const double f = F(start_idx + i);
    const double sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) G.col(i) *= rploidy / sd;
  }
}

void Data::pcangsd_standardize_E(const Mat2D& U, const Mat1D& svals, const Mat2D& VT) {
  if (params.scale != -9) return;
  cao.print(tick.date(), "begin to standardize the matrix for pcangsd procedure");
  const Mat2D US = U * svals.asDiagonal();
  Mat2D diag_threads = Mat2D::Zero(nsamples, omp_get_max_threads());  // a column per thread
  for_each_product_column(US, VT, 0, nsnps, [&](Eigen::Index jj, const auto& recon) {
    auto diag_private = diag_threads.col(omp_get_thread_num());
    {
      const uint j = (uint)jj;
      double p0, p1, p2, pt, pSum, tmp;
      const uint original = unpermuted_snp_index(j);
      const double f = F(original);
      double norm = sqrt(2.0 * f * (1.0 - f));
      uint s = params.filterSNP ? keepSNPs[original] : original;
      for (uint i = 0; i < nsamples; i++) {
        // Rescale individual allele frequencies
        pt = recon(i);
        pt = (pt + 2.0 * f) / 2.0;
        pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
        // Update e
        p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        pSum = p0 + p1 + p2;
        G(i, j) = (p1 + 2 * p2) / pSum - 2.0 * f;
        if (norm > VAR_TOL) G(i, j) /= norm;

        // Update diag
        tmp = (0.0 - 2.0 * f) * (0.0 - 2.0 * f) * (p0 / pSum);
        tmp += (1.0 - 2.0 * f) * (1.0 - 2.0 * f) * (p1 / pSum);
        tmp += (2.0 - 2.0 * f) * (2.0 - 2.0 * f) * (p2 / pSum);
        if (norm > VAR_TOL) diag_private(i) += tmp / (2.0 * f * (1.0 - f));  // f = 0 gave inf
      }
    }
  });
  Dc = diag_threads.rowwise().sum();  // in thread order, not in the order a critical section was entered
}

void Data::predict_missing_E(const Mat2D& U, uint64 start_idx, uint64 stop_idx) {
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
      if (!C((uint64)j * nsamples + i))
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
