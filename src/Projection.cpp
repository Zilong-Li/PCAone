/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Projection.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Projection.hpp"

#include <limits>

#include "Cmd.hpp"
#include "Common.hpp"
#include "Data.hpp"
#include "Utils.hpp"

Mat2D solve_bootstrap_projection_no_missing(const Mat2D& design, const Mat2D& G, const std::vector<uint>& counts) {
  const int M = design.rows();
  const int K = design.cols();
  const int N = G.rows();
  Mat2D ata = Mat2D::Zero(K, K);
  Mat2D atg = Mat2D::Zero(K, N);
  for (int j = 0; j < M; ++j) {
    const uint w = counts[j];
    if (w == 0) continue;
    ata.noalias() += (double)w * design.row(j).transpose() * design.row(j);
    atg.noalias() += (double)w * design.row(j).transpose() * G.col(j).transpose();
  }
  return ata.ldlt().solve(atg).transpose();
}

Mat2D solve_bootstrap_projection_missing(const Mat2D& design,
                                         const ArrBool& C,
                                         const Mat2D& G,
                                         const std::vector<uint>& counts) {
  const int M = design.rows();
  const int K = design.cols();
  const int N = G.rows();
  Mat2D U(N, K);
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    Mat2D ata = Mat2D::Zero(K, K);
    Mat1D atg = Mat1D::Zero(K);
    int observed = 0;
    for (int j = 0; j < M; ++j) {
      const uint w = counts[j];
      if (w == 0 || C(j * N + i)) continue;
      ata.noalias() += (double)w * design.row(j).transpose() * design.row(j);
      atg.noalias() += (double)w * design.row(j).transpose() * G(i, j);
      observed += w;
    }
    if (observed < K) {
      U.row(i).setConstant(std::numeric_limits<double>::quiet_NaN());
    } else {
      U.row(i) = ata.ldlt().solve(atg).transpose();
    }
  }
  return U;
}

void write_projection_bootstrap_stats(
    const Mat2D& design, const ArrBool& C, const Mat2D& G, const Mat2D& U0, uint nreps, const Param& params) {
  const int M = design.rows();
  const int N = G.rows();
  const int K = design.cols();
  if (nreps < 2) cao.error("--project-bootstrap must be at least 2");
  if (M < K) cao.error("projection bootstrap requires at least as many matched SNPs as PCs");

  cao.print(tick.date(), "run projection SNP bootstrap with", nreps, "replicates");
  const bool has_missing = C.size() && C.count() > 0;
  std::mt19937_64 rng(params.seed);
  std::uniform_int_distribution<int> snp_dist(0, M - 1);

  Mat2D sum = Mat2D::Zero(N, K);
  Mat2D sumsq = Mat2D::Zero(N, K);
  Mat2D diff_sumsq = Mat2D::Zero(N, K);
  Mat2D minv = Mat2D::Constant(N, K, std::numeric_limits<double>::infinity());
  Mat2D maxv = Mat2D::Constant(N, K, -std::numeric_limits<double>::infinity());

  std::vector<std::pair<int, int>> pairs;
  pairs.reserve(static_cast<size_t>(K) * (K - 1) / 2);
  for (int a = 0; a < K; ++a) {
    for (int b = a + 1; b < K; ++b) pairs.emplace_back(a, b);
  }
  Mat2D cross_sum = Mat2D::Zero(N, pairs.size());

  std::ofstream raw;
  if (params.project_bootstrap_save) {
    raw.open(params.fileout + ".proj.bootstrap.eigvecs");
    if (!raw.is_open()) cao.error("can not open " + params.fileout + ".proj.bootstrap.eigvecs");
    raw << "replicate\tsample";
    for (int k = 0; k < K; ++k) raw << "\tPC" << k + 1;
    raw << '\n';
  }

  std::vector<uint> counts(M);
  for (uint r = 0; r < nreps; ++r) {
    std::fill(counts.begin(), counts.end(), 0);
    for (int draw = 0; draw < M; ++draw) counts[snp_dist(rng)]++;

    Mat2D Ub = has_missing ? solve_bootstrap_projection_missing(design, C, G, counts)
                           : solve_bootstrap_projection_no_missing(design, G, counts);

    for (int k = 0; k < K; ++k) {
      if (Ub.col(k).allFinite() && Ub.col(k).dot(U0.col(k)) < 0.0) Ub.col(k) *= -1.0;
    }

    sum += Ub;
    sumsq += Ub.array().square().matrix();
    diff_sumsq += (Ub - U0).array().square().matrix();
    minv = minv.cwiseMin(Ub);
    maxv = maxv.cwiseMax(Ub);
    for (int p = 0; p < (int)pairs.size(); ++p) {
      cross_sum.col(p).array() += Ub.col(pairs[p].first).array() * Ub.col(pairs[p].second).array();
    }

    if (raw.is_open()) {
      for (int i = 0; i < N; ++i) {
        raw << r + 1 << '\t' << i + 1;
        for (int k = 0; k < K; ++k) raw << '\t' << Ub(i, k);
        raw << '\n';
      }
    }
  }

  const double R = (double)nreps;
  Mat2D mean = sum / R;
  Mat2D variance = (sumsq.array() - (sum.array().square() / R)).max(0.0) / std::max(1.0, R - 1.0);
  Mat2D sd = variance.array().sqrt().matrix();
  Mat2D rmsd = (diff_sumsq.array() / R).sqrt().matrix();

  std::ofstream out(params.fileout + ".proj.bootstrap.tsv");
  if (!out.is_open()) cao.error("can not open " + params.fileout + ".proj.bootstrap.tsv");
  out << "sample\tpc\tbaseline\tmean\tbootstrap_se\tmin\tmax\trmsd\n";
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      out << i + 1 << "\tPC" << k + 1 << '\t' << U0(i, k) << '\t' << mean(i, k) << '\t' << sd(i, k) << '\t'
          << minv(i, k) << '\t' << maxv(i, k) << '\t' << rmsd(i, k) << '\n';
    }
  }

  std::ofstream covout(params.fileout + ".proj.bootstrap.cov.tsv");
  if (!covout.is_open()) cao.error("can not open " + params.fileout + ".proj.bootstrap.cov.tsv");
  covout << "sample\tpc_x\tpc_y\tbaseline_x\tbaseline_y\tmean_x\tmean_y\tvar_x\tvar_y\tcov\tcorr\n";
  for (int i = 0; i < N; ++i) {
    for (int p = 0; p < (int)pairs.size(); ++p) {
      const int a = pairs[p].first;
      const int b = pairs[p].second;
      const double cov = (cross_sum(i, p) - R * mean(i, a) * mean(i, b)) / std::max(1.0, R - 1.0);
      const double denom = std::sqrt(variance(i, a) * variance(i, b));
      const double corr = denom > 0.0 ? cov / denom : std::numeric_limits<double>::quiet_NaN();
      covout << i + 1 << "\tPC" << a + 1 << "\tPC" << b + 1 << '\t' << U0(i, a) << '\t' << U0(i, b) << '\t'
             << mean(i, a) << '\t' << mean(i, b) << '\t' << variance(i, a) << '\t' << variance(i, b) << '\t' << cov
             << '\t' << corr << '\n';
    }
  }
  cao.print(tick.date(), "projection bootstrap diagnostics saved to", params.fileout + ".proj.bootstrap.tsv");
}

void solve_projection_scores(const Mat2D& V, const ArrBool& C, const Mat2D& G, Mat2D& U) {
  if (U.rows() == 0 || U.cols() == 0) return;
  double p_miss = C.size() ? (double)C.count() / (double)C.size() : 0.0;
  if (p_miss == 0.0) {
    Eigen::ColPivHouseholderQR<Mat2D> qr(V);
#pragma omp parallel for
    for (uint i = 0; i < (uint)U.rows(); i++) {
      U.row(i) = qr.solve(G.row(i).transpose());
    }
  } else {
#pragma omp parallel for
    for (uint i = 0; i < (uint)U.rows(); i++) {
      Int1D idx;
      for (int j = 0; j < V.rows(); j++) {
        if (!C(j * U.rows() + i)) idx.push_back(j);
      }
      U.row(i) = V(idx, Eigen::all).colPivHouseholderQr().solve(G(i, idx).transpose());
    }
  }
}

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
  } else if (params.file_t == FileType::PGEN) {
    match = match_pvar_to_mbim(params.filein + ".pvar", params.filebim);
    if (match.bim_indices.empty())
      cao.error("no overlapped SNPs found between " + params.filein + ".pvar and " + params.filebim);
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
    cao.warn("SNP info is not fully identical between input and reference .mbim");
    cao.print(tick.date(), "projection will use", match.bim_indices.size(), " overlapped sites. there are",
              data->flipSNPs.size(), " flipped alleles");
    if (!data->flipSNPs.empty())
      cao.warn(data->flipSNPs.size(), " SNPs have flipped ref/alt alleles and will be corrected");
  }
  cao.print(tick.date(), "run projection");
  data->prepare();  // read AF and resize F to matched size
  if (params.project != 3) data->standardize_E();
  cao.print(tick.date(), "start parsing V:", params.fileV, ", S:", params.fileS);
  uint nsamples, nsnps;
  Mat1D S;
  read_sigvals(params.fileS, nsamples, nsnps, S);
  // target number of PCs for getting individual allele frequency
  const int K = fmin(S.size(), params.k);
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
  double p_miss = data->C.size() ? (double)data->C.count() / (double)data->C.size() : data->p_miss;

  if (params.project == 1) {
    if (p_miss > 0) cao.warn("there are missing genotypes. recommend using --project 2 or 3.");
    // get 1 / Singular = sqrt(Eigen * M)
    V = V * (S.array().inverse().matrix().asDiagonal());
    // G V = U D
    U = data->G * V;
  } else if (params.project == 2) {
    V = V * S.asDiagonal();
    if (p_miss == 0.0) cao.warn("there is no missing genotypes");
    solve_projection_scores(V, data->C, data->G, U);
    if (params.project_bootstrap > 0)
      write_projection_bootstrap_stats(V, data->C, data->G, U, params.project_bootstrap, params);
  } else if (params.project == 3) {
    // project == 3: iterative GL-aware projection (EM)
    // E-step: update expected G (0, 1) using individual allele frequencies PI
    // M-step: solve U with new G
    if (params.file_t != FileType::BEAGLE) cao.error("--project 3 requires BEAGLE genotype likelihood input");
    const bool filter = !data->keepSNPs.empty();

    V = V * S.asDiagonal();  // VS
    solve_projection_scores(V, data->C, data->G, U);

    cao.print(tick.date(), "run EM to update expected G and solve U iteratively");
    // NOTE: First, we map G to domain [0,1]; Second, we can't standarize/scale G.
    for (uint iter = 0; iter < params.maxiter; ++iter) {
      Mat2D Uprev = U;

      // E-step: update G using individual allele frequencies
#pragma omp parallel for
      for (uint j = 0; j < data->nsnps; ++j) {
        // const double norm = sqrt(2.0 * data->F(j) * (1.0 - data->F(j)));
        uint s = filter ? data->keepSNPs[j] : (uint)j;
        for (uint i = 0; i < data->nsamples; ++i) {
          double pt = 0.0;
          for (int k = 0; k < K; ++k) {
            pt += U(i, k) * V(j, k);
          }
          // if (params.scale == SCALE_STANDARDIZE_GENETIC && norm > VAR_TOL) pt *= norm;
          pt = fmin(fmax(pt * 0.5 + data->F(j), 1e-4), 1.0 - 1e-4);  //
          const double p0 = data->P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
          const double p1 = data->P(2 * i + 1, s) * 2.0 * pt * (1.0 - pt);
          const double p2 = (1.0 - data->P(2 * i + 0, s) - data->P(2 * i + 1, s)) * pt * pt;
          const double psum = p0 + p1 + p2;
          if (!std::isfinite(psum) || psum <= 0.0) {
            data->C[j * data->nsamples + i] = 1;
            data->G(i, j) = 0.0;
            continue;
          }
          data->C[j * data->nsamples + i] = 0;
          data->G(i, j) = (p1 + 2.0 * p2) / (2.0 * psum) - data->F(j);  // domain (0,1)
          // if (params.scale == SCALE_STANDARDIZE_GENETIC && norm > VAR_TOL) data->G(i, j) /= norm;
        }
      }

      // M-step: solve for U using expected G
      solve_projection_scores(V, data->C, data->G, U);
      // // flip signs
      // for (int k = 0; k < K; ++k) {
      //   if (U.col(k).dot(Uprev.col(k)) < 0.0) U.col(k) *= -1.0;
      // }
      double denom = Uprev.norm();
      if (denom < 1e-12) denom = 1.0;
      double diff = (U - Uprev).norm() / denom;
      cao.print(tick.date(), "GL projection iter", iter + 1, ", diff =", diff);
      if (diff < params.tolem) break;
    }
  } else {
    cao.error("unsupported --project mode: " + std::to_string(params.project));
  }

  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  std::ofstream outu(params.fileout + ".eigvecs");
  if (outu.is_open()) outu << U.format(fmt) << '\n';
}
