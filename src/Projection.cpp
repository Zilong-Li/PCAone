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
      if (w == 0 || C((Eigen::Index)j * N + i)) continue;
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
  PortableRng rng(params.seed);  // uniform_int_distribution differs between standard libraries

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
    for (int draw = 0; draw < M; ++draw) counts[rng.below(M)]++;

    Mat2D Ub = has_missing ? solve_bootstrap_projection_missing(design, C, G, counts)
                           : solve_bootstrap_projection_no_missing(design, G, counts);

    // no sign alignment: the design V*S is fixed, so Ub is an ordinary
    // least-squares estimate with no arbitrary sign. Flipping it towards U0
    // would fold the replicates onto the baseline's side and bias the mean,
    // variance and covariance.

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
    // with fewer called sites than PCs the least-squares system has no unique
    // solution; the pivoted QR returned 0 for every PC (a sample with no calls
    // at all was projected to the origin), indistinguishable from a real score
    uint64 nunder = 0;
#pragma omp parallel for reduction(+ : nunder)
    for (uint i = 0; i < (uint)U.rows(); i++) {
      Int1D idx;
      for (int j = 0; j < V.rows(); j++) {
        if (!C((Eigen::Index)j * U.rows() + i)) idx.push_back(j);
      }
      if ((Eigen::Index)idx.size() < V.cols()) {
        U.row(i).setConstant(std::numeric_limits<double>::quiet_NaN());
        ++nunder;
        continue;
      }
      U.row(i) = V(idx, Eigen::all).colPivHouseholderQr().solve(G(i, idx).transpose());
    }
    if (nunder > 0)
      cao.warn(std::to_string(nunder) + " sample(s) have fewer called sites than PCs and are written as NA");
  }
}

// Per-site factor a_j for --project 3, where X_ref = a_j * (g/2 - f) is what the
// reference PCA decomposed in terms of PCAone's 0..1 centred coding. The EM
// regresses a_j * (E[g]/2 - f) on V*S, so U lands on the reference's scale, and
// maps the reconstruction back with pi = f + (U*S*V')_j / a_j. PCAngsd's
// Algorithm 1 is the a_j = 2 case: pi = f + recon/2 on centred dosages.
//
//   standardized (scale=-9)      sqrt(ploidy)/sd   as standardize_E_ref()
//   centred dosages (gscale=2)   2                 pcangsd: dosage - 2f
//   centred 0..1                 1                 -D/--ld, --scale 0
static Mat1D gl_site_scale(const Mat1D& F, bool standardize, const UsvTransform& t) {
  Mat1D a = Mat1D::Constant(F.size(), (double)t.gscale);
  if (!standardize) return a;
  const double rploidy = sqrt((double)t.ploidy);
  for (Eigen::Index j = 0; j < F.size(); ++j) {
    const double sd = sqrt(F(j) * (1.0 - F(j)));
    if (sd > VAR_TOL) a(j) = rploidy / sd;
  }
  return a;
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
  cao.print(tick.date(), "start parsing V:", params.fileV, ", S:", params.fileS);
  uint nsamples, nsnps;
  Mat1D S;
  // read the reference's transform before touching G: the target has to be put
  // on the same scale as the matrix that produced V and S, which is not
  // necessarily what this run's -C/--scale says. --project 3 also accepts a
  // dosage-scale (pcangsd) reference, and applies the transform per site
  // inside its EM (gl_site_scale) rather than to G up front.
  UsvTransform usv;
  read_sigvals(params.fileS, nsamples, nsnps, S, &usv);
  data->prepare();  // read AF and resize F to matched size
  const bool standardize = data->resolve_ref_scaling(usv, params.fileS, params.project == 3);
  if (params.project != 3 && standardize) data->standardize_E_ref(usv);
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
    // E-step: update expected G using individual allele frequencies PI
    // M-step: solve U with new G
    if (params.file_t != FileType::BEAGLE) cao.error("--project 3 requires BEAGLE genotype likelihood input");
    const bool filter = !data->keepSNPs.empty();

    // FileBeagle::read_all() leaves G = E[g]/2 - f; put it on the scale of the
    // matrix the reference decomposed before the first solve, as modes 1 and 2 do
    const Mat1D a = gl_site_scale(data->F, standardize, usv);
    data->G.array().rowwise() *= a.transpose().array();

    V = V * S.asDiagonal();  // VS
    solve_projection_scores(V, data->C, data->G, U);

    cao.print(tick.date(), "run EM to update expected G and solve U iteratively");
    for (uint iter = 0; iter < params.maxiter; ++iter) {
      Mat2D Uprev = U;

      // E-step: update G using individual allele frequencies, from U * V' formed
      // a panel of sites at a time instead of a scalar loop per genotype
      for_each_product_column(U, V.transpose(), 0, data->nsnps, [&](Eigen::Index jj, const auto& recon) {
        const uint j = (uint)jj;
        uint s = filter ? data->keepSNPs[j] : (uint)j;
        for (uint i = 0; i < data->nsamples; ++i) {
          double pt = recon(i);
          // pt is on the reference's scale, the same as G; undo a_j to get pi
          pt = fmin(fmax(pt / a(j) + data->F(j), 1e-4), 1.0 - 1e-4);
          const double p0 = data->P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
          const double p1 = data->P(2 * i + 1, s) * 2.0 * pt * (1.0 - pt);
          const double p2 = (1.0 - data->P(2 * i + 0, s) - data->P(2 * i + 1, s)) * pt * pt;
          const double psum = p0 + p1 + p2;
          if (!std::isfinite(psum) || psum <= 0.0) {
            data->C[(uint64)j * data->nsamples + i] = 1;
            data->G(i, j) = 0.0;
            continue;
          }
          data->C[(uint64)j * data->nsamples + i] = 0;
          data->G(i, j) = a(j) * ((p1 + 2.0 * p2) / (2.0 * psum) - data->F(j));
        }
      });

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

  std::ofstream outu(params.fileout + ".eigvecs");
  if (outu.is_open()) write_rows(outu, U);  // NA for a sample without a score
}
