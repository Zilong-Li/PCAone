/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Projection.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Projection.hpp"

#include "Cmd.hpp"
#include "Utils.hpp"

namespace {

inline double expected_centered_genotype(const Mat2D& gls, const Mat1D& freqs, uint sample_idx, uint raw_snp_idx,
                                         uint matched_snp_idx, double allele_freq) {
  allele_freq = std::fmin(std::fmax(allele_freq, 1e-4), 1.0 - 1e-4);
  const double gl0 = gls(2 * sample_idx + 0, raw_snp_idx);
  const double gl1 = gls(2 * sample_idx + 1, raw_snp_idx);
  const double p0 = gl0 * (1.0 - allele_freq) * (1.0 - allele_freq);
  const double p1 = gl1 * 2.0 * allele_freq * (1.0 - allele_freq);
  const double p2 = (1.0 - gl0 - gl1) * allele_freq * allele_freq;
  const double posterior_sum = p0 + p1 + p2;
  if (posterior_sum <= 0.0) return -2.0 * freqs(matched_snp_idx);
  return (p1 + 2.0 * p2) / posterior_sum - 2.0 * freqs(matched_snp_idx);
}

}  // namespace

/**
 * options:
 * 1: simple, assume no missingness
 * 2: like smartPCA, solving g=Vx, can take missing genotypes
 * 3: iterative GL-aware projection (EM): alternates between updating individual allele frequencies
 *    from current PC scores and re-solving for PC scores with updated expected genotypes (BEAGLE only)
   // NOTE: we don't support out-of-core for projection.
 */
void run_projection(Data* data, const Param& params) {
  BimMatch match;
  if (params.file_t == FileType::BEAGLE) {
    match = match_beagle_to_mbim(params.filein, params.filebim);
    if (match.bim_indices.empty())
      cao.error("no overlapped SNPs found between " + params.filein + " and " + params.filebim);
  } else {
    match = match_bim_to_mbim(params.filein + ".bim", params.filebim);
    if (match.bim_indices.empty())
      cao.error("no overlapped SNPs found between " + params.filein + ".bim and " + params.filebim);
  }
  if (!match.identical) {
    data->keepSNPs = match.bim_indices;
    data->keepRefSNPs = match.mbim_indices;
    for (int k = 0; k < (int)match.flip.size(); ++k) {
      if (match.flip[k]) data->flipSNPs.push_back(k);
    }
    cao.warn("SNP info is not fully identical between input and reference .mbim;\n projection will use ",
             match.bim_indices.size(), " overlapped sites only");
    if (!data->flipSNPs.empty())
      cao.warn(data->flipSNPs.size(), " SNPs have flipped ref/alt alleles and will be corrected");
  }
  cao.print(tick.date(), "run projection");
  data->prepare();
  if (params.project != 3) data->standardize_E();
  cao.print(tick.date(), "start parsing V:", params.fileV, ", S:", params.fileS);
  uint nsamples, nsnps;
  Mat1D S;
  read_sigvals(params.fileS, nsamples, nsnps, S);
  int K = fmin(S.size(), params.k);
  Mat2D V = read_eigvecs(params.fileV, nsnps, K);
  if (!match.identical) {
    Mat2D V_overlap(match.mbim_indices.size(), K);
    for (int i = 0; i < (int)match.mbim_indices.size(); ++i) {
      V_overlap.row(i) = V.row(match.mbim_indices[i]);
      if (match.flip[i]) V_overlap.row(i) = -V_overlap.row(i);
    }
    V = V_overlap;
  }
  Mat2D U(data->nsamples, K);

  if (params.project == 1) {
    if (data->p_miss > 0) cao.warn("there are missing genotypes. recommend using --project 2 or 3.");
    // get 1 / Singular = sqrt(Eigen * M)
    V = V * (S.array().inverse().matrix().asDiagonal());
    // G V = U D
    U = data->G * V;
  } else if (params.project == 2) {
    V = V * S.asDiagonal();
    if (data->p_miss == 0.0) {
      Eigen::ColPivHouseholderQR<Mat2D> qr(V);
#pragma omp parallel for
      for (uint i = 0; i < data->nsamples; i++) {
        // Vx = g
        U.row(i) = qr.solve(data->G.row(i).transpose());
      }
    } else {
#pragma omp parallel for
      for (uint i = 0; i < data->nsamples; i++) {
        // find non-missing snps for sample i
        Int1D idx;
        for (int j = 0; j < V.rows(); j++) {
          if (!data->C(j * data->nsamples + i)) idx.push_back(j);
        }
        // Vx = g
        U.row(i) = V(idx, Eigen::all).colPivHouseholderQr().solve(data->G(i, idx).transpose());
      }
    }
  } else if (params.project == 3) {
    // project == 3: iterative GL-aware projection (EM)
    // E-step: update expected genotype G using per-sample allele frequencies from U
    // M-step: solve V*S*u_i = G_std_i (standardized with reference F)
    if (params.file_t != FileType::BEAGLE)
      cao.error("--project 3 requires BEAGLE genotype likelihood input");
    // data->G is centered (expected_dosage - 2F), data->P holds raw GLs
    const int nsnps = (int)data->nsnps;
    const int nsamples = (int)data->nsamples;
    const bool filter = !data->keepSNPs.empty();

    // sd(j) = sqrt(F*(1-F)/ploidy): convert G_centered <-> G_std
    Mat1D sd(nsnps);
    for (int j = 0; j < nsnps; ++j)
      sd(j) = std::sqrt(data->F(j) * (1.0 - data->F(j)) / params.ploidy);

    const auto sdiag = S.asDiagonal();
    // V*S is fixed across iterations
    Mat2D VS = V * sdiag;
    Eigen::ColPivHouseholderQR<Mat2D> qr(VS);
    Mat2D G_std(data->nsamples, nsnps);
    Mat2D pred(data->nsamples, nsnps);

    U.setZero();
    for (uint iter = 0; iter < params.maxiter; ++iter) {
      Mat2D U_prev = U;

      // M-step: standardize G and solve for U
      G_std = data->G;
      for (int j = 0; j < nsnps; ++j)
        if (sd(j) > VAR_TOL) G_std.col(j) /= sd(j);
#pragma omp parallel for
      for (int i = 0; i < nsamples; ++i)
        U.row(i) = qr.solve(G_std.row(i).transpose());

      // E-step: update G using individual allele frequencies
      // pred(i,j) = predicted G_std; pt = (pred*sd + 2F)/2 = individual allele freq
      pred.noalias() = U * sdiag * V.transpose();  // nsamples x nsnps
#pragma omp parallel for
      for (int j = 0; j < nsnps; ++j) {
        uint s = filter ? data->keepSNPs[j] : (uint)j;
        for (int i = 0; i < nsamples; ++i) {
          const double pt = (pred(i, j) * sd(j) + 2.0 * data->F(j)) / 2.0;
          data->G(i, j) = expected_centered_genotype(data->P, data->F, i, s, j, pt);
        }
      }

      double delta = (U - U_prev).norm() / (U_prev.norm() + 1e-10);
      cao.print(tick.date(), "projection iter", iter + 1, ", delta =", delta);
      if (iter > 0 && delta < params.tol) break;
    }
  } else {
    cao.error("unsupported --project mode: " + std::to_string(params.project));
  }

  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  std::ofstream outu(params.fileout + ".eigvecs");
  if (outu.is_open()) outu << U.format(fmt) << '\n';
}
