/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/EvalAdmix.cpp
 * @author      Anders Albrechtsen
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 *
 * Correlation of residuals for evaluating the fit of a PCA/admixture model
 * (evalAdmix, Garcia-Erill & Albrechtsen 2020; projection estimator of
 * van Waaij et al. 2023, Genetics 225:iyad157).
 *
 * The statistic is  corres = bhat - chat  where, with P the projection onto
 * [PC_1..PC_k, 1] and R = G(I-P) the residuals,
 *
 *   bhat = cov2cor( Rtilde' Rtilde ),  Rtilde = column-centred R
 *   chat = cov2cor( (I-P) Dhat (I-P) ),  Dhat = diag(mean heterozygosity)
 *
 * Assuming no missing genotypes, Rtilde'Rtilde can be written without ever
 * forming R:
 *
 *   Rtilde'Rtilde = (I-P) [ G'G - M gbar gbar' ] (I-P)
 *
 * so one streaming pass accumulating the N x N Gram matrix G'G, the per-sample
 * mean genotype gbar and the per-sample mean heterozygosity is sufficient.
 * Memory is O(N^2), independent of the number of sites.
 *
 * Note PCAone codes genotypes on the 0..1 scale (BED2GENO = {1, NA, 0.5, 0}),
 * i.e. x = g/2. Both bhat and chat are correlation matrices, so the constant
 * factor between the x and g scales cancels and no rescaling is needed.
 ******************************************************************************/
#include "EvalAdmix.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "Utils.hpp"

// Read a whitespace/tab separated numeric matrix.
// NOTE: not using Utils' read_usv() here -- it maps a row-major buffer with
// Eigen::Map<Mat2D> (column-major), which scrambles any file with >1 column.
static Mat2D read_matrix_rowmajor(const std::string& path) {
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("can not open file: " + path);
  std::vector<std::vector<double>> rows;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    std::istringstream iss(line);
    std::vector<double> r;
    double v;
    while (iss >> v) r.push_back(v);
    if (!r.empty()) rows.push_back(std::move(r));
  }
  if (rows.empty()) cao.error("empty matrix file: " + path);
  const Eigen::Index n = rows.size(), k = rows[0].size();
  Mat2D M(n, k);
  for (Eigen::Index i = 0; i < n; ++i) {
    if ((Eigen::Index)rows[i].size() != k) cao.error("ragged matrix in file: " + path);
    for (Eigen::Index j = 0; j < k; ++j) M(i, j) = rows[i][j];
  }
  return M;
}

// turn a symmetric PSD matrix into a correlation matrix
static Mat2D cov2cor(const Mat2D& C) {
  const Eigen::Index n = C.rows();
  Mat1D s(n);
  for (Eigen::Index i = 0; i < n; ++i) s(i) = std::sqrt(std::max(C(i, i), 1e-300));
  Mat2D R(n, n);
  for (Eigen::Index i = 0; i < n; ++i)
    for (Eigen::Index j = 0; j < n; ++j) R(i, j) = C(i, j) / (s(i) * s(j));
  return R;
}

static void write_matrix(const std::string& fn, const Mat2D& M, const std::vector<std::string>& ids) {
  std::ofstream ofs(fn);
  if (!ofs.is_open()) cao.error("can not open file for writing: " + fn);
  ofs << std::fixed << std::setprecision(6);
  if (!ids.empty()) {
    for (size_t i = 0; i < ids.size(); ++i) ofs << (i ? "\t" : "") << ids[i];
    ofs << "\n";
  }
  for (Eigen::Index i = 0; i < M.rows(); ++i) {
    for (Eigen::Index j = 0; j < M.cols(); ++j) ofs << (j ? "\t" : "") << M(i, j);
    ofs << "\n";
  }
}

void run_evaladmix(Data* data, const Param& params) {
  const Eigen::Index N = data->nsamples;
  const uint M = data->nsnps;

  // ---- 1. principal component scores -------------------------------------
  std::string fpcs = params.fileU.empty() ? params.fileout + ".eigvecs" : params.fileU;
  Mat2D U = read_matrix_rowmajor(fpcs);  // N x kmax
  cao.print(tick.date(), "evalAdmix: read", U.rows(), "x", U.cols(), "PC scores from", fpcs);
  if (U.rows() != N) cao.error("number of samples in .eigvecs does not match the genotype file");
  Eigen::Index k = params.evaladmix_k > 0 ? params.evaladmix_k : U.cols();
  if (k > U.cols()) cao.error("--evaladmix-k is larger than the number of PCs available");
  cao.print(tick.date(), "evalAdmix: using", k, "PC(s) + intercept =", k + 1, "dimensions");

  // ---- 2. one streaming pass: Gram matrix, mean genotype, heterozygosity --
  Mat2D A = Mat2D::Zero(N, N);  // G'G
  Mat1D b = Mat1D::Zero(N);     // sum_s g_s
  Mat1D d = Mat1D::Zero(N);     // sum_s g_s .* (1 - g_s)   (0..1 scale)
  tick.clock();
  if (!params.out_of_core) {
    // G is nsamples x nsnps, raw (params.center forced false, not standardized)
    const Mat2D& G = data->G;
    A.noalias() = G * G.transpose();
    b = G.rowwise().sum();
    d = (G.array() * (1.0 - G.array())).rowwise().sum().matrix();
  } else {
    // Out-of-core. read_block_initial() estimates F on the fly and always
    // returns *centred* genotypes (it uses centered_geno_lookup regardless of
    // params.center), so both must be allocated here, and the centring has to
    // be undone to recover heterozygosity.
    //
    // A and b need no such correction: per-site centring subtracts c_s * 1 from
    // column s, and every resulting term carries a factor (I-P)1 = 0 because
    // the projection contains the intercept. So (I-P)[A - M gbar gbar'](I-P) is
    // identical for centred and raw genotypes. Only d is affected.
    data->F = Mat1D::Zero(data->nsnps);
    data->centered_geno_lookup = Arr2D::Zero(4, data->nsnps);
    data->check_file_offset_first_var();
    for (uint bi = 0; bi < data->nblocks; ++bi) {
      data->read_block_initial(data->start[bi], data->stop[bi], false);
      const Mat2D& G = data->G;
      A.noalias() += G * G.transpose();
      b += G.rowwise().sum();
      for (Eigen::Index i = 0; i < G.cols(); ++i) {
        const double f = data->F(data->start[bi] + i);
        d.array() += (G.col(i).array() + f) * (1.0 - (G.col(i).array() + f));
      }
    }
  }
  b /= (double)M;
  d /= (double)M;
  cao.print(tick.date(), "evalAdmix: accumulated summary statistics over", M, "sites in", tick.reltime(), "seconds");

  // ---- 3. projection onto [PCs, intercept] -------------------------------
  Mat2D V(N, k + 1);
  V.leftCols(k) = U.leftCols(k);
  V.col(k).setOnes();
  Mat2D VtV = V.transpose() * V;
  Mat2D P = V * VtV.completeOrthogonalDecomposition().pseudoInverse() * V.transpose();
  Mat2D IP = Mat2D::Identity(N, N) - P;

  // ---- 4. bhat and chat --------------------------------------------------
  Mat2D Cw = IP * (A - (double)M * b * b.transpose()) * IP;
  Mat2D bhat = cov2cor(Cw);
  Mat2D chat = cov2cor(IP * d.asDiagonal() * IP);
  Mat2D corres = bhat - chat;
  // a correlation cannot leave [-1, 1]
  corres = corres.cwiseMax(-1.0).cwiseMin(1.0);
  corres.diagonal().setZero();

  // ---- 5. output ---------------------------------------------------------
  std::vector<std::string> ids;
  if (params.file_t == FileType::PLINK) {
    std::ifstream ffam(params.filein + ".fam");
    std::string line, fid, iid;
    while (std::getline(ffam, line)) {
      std::istringstream iss(line);
      if (iss >> fid >> iid) ids.push_back(iid);
    }
  }
  if ((Eigen::Index)ids.size() != N) ids.clear();
  write_matrix(params.fileout + ".corres", corres, ids);
  Mat2D kin = corres * 0.5;  // correlation of residuals estimates 2*phi
  kin.diagonal().setZero();
  write_matrix(params.fileout + ".kinship", kin, ids);
  cao.print(tick.date(), "evalAdmix: correlation of residuals saved to", params.fileout + ".corres");
  cao.print(tick.date(), "evalAdmix: kinship (corres/2) saved to", params.fileout + ".kinship");
}
