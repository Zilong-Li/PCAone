#include "Selection.hpp"

#include "LinAlg.hpp"
#include "Stats.hpp"
#include "TextIO.hpp"
#include "Utils.hpp"

namespace {
// A site whose genotype column is constant carries no information: its loading
// is 0 and its residual variance is 0, so neither statistic is defined for it.
// Reporting 0 (p = 1) instead put it into the robust covariance, into the
// genomic-inflation median and into the output as a well-behaved site. pcadapt
// drops these (`zscores[pass, ]`) and reports NA; so do we now.
void warn_dropped_sites(uint64 ndrop, uint64 ntotal, const std::string& what) {
  if (ndrop > 0)
    cao.warn(ndrop, " of ", ntotal,
             " sites have no " + what +
                 " (monomorphic, or fully explained by the PCs). they are reported as NA and left out of the "
                 "statistic, the robust covariance and the inflation factor. consider --maf to remove them.");
}

// --maf removes sites before the scan, but every output is read by position
// against the .bim/.pvar and carries no IDs, so the rows after the first removed
// site belonged to later sites. Put the removed sites back as NA rows: line r
// of an output is row rows[r] of the statistic, all NA where it is negative.
std::vector<Eigen::Index> all_sites(const Data* data) {
  std::vector<Eigen::Index> rows;
  if (data->keepSNPs.empty()) return rows;
  rows.assign(data->nsnps_all, -1);
  for (size_t j = 0; j < data->keepSNPs.size(); ++j) rows[data->keepSNPs[j]] = (Eigen::Index)j;
  return rows;
}

void write_sites(const std::string& path,
                 const std::string& header,
                 const Eigen::Ref<const Mat2D>& M,
                 const std::vector<Eigen::Index>& rows) {
  TextFormat fmt;
  fmt.na = "NA";
  write_matrix(path, M, fmt, header, rows.empty() ? nullptr : &rows);
}

// FastPCA/Galinsky statistic: M * v_{jk}^2 ~ chi-squared(df=1) under the null,
// where v_{jk} is the SNP loading for SNP j along PC k.
void galinsky_selection_stat(Mat2D& V) {
#pragma omp parallel for
  for (int j = 0; j < V.rows(); j++) {
    for (int i = 0; i < V.cols(); i++) {
      V(j, i) = V(j, i) * V(j, i) * V.rows();
    }
  }
}

// ---------------------------------------------------------------------------
// Orthogonalized Gnanadesikan-Kettenring (Maronna & Zamar 2002), matching
// bigutilsr::dist_ogk() -- which is what pcadapt calls -- with its defaults
// niter = 2, beta = 0.9 and the tau scale of Yohai & Zamar (1998).
//
// The plain pairwise GK estimator this replaces was not even scale
// equivariant: rescaling one column of z-scores changed the Mahalanobis
// distances, so the statistic depended on an arbitrary normalization. The
// orthogonalization step fixes that, and the final reweighted estimate is a
// classical covariance of the retained rows, so it is positive definite by
// construction rather than by flooring eigenvalues.
// ---------------------------------------------------------------------------

constexpr double TAU_C1 = 4.5;
constexpr double TAU_C2 = 3.0;
constexpr double QNORM_3_4 = 0.6744897501960817;  // qnorm(3/4)
constexpr double INV_SQRT2 = 0.7071067811865475244;
constexpr double SQRT_2PI = 2.5066282746310005024;

// E[rho] for the consistency factor of the tau scale at the normal, exactly as
// robustbase::scaleTau2 computes it: Es2(c2) = Erho(c2 * qnorm(3/4)).
double tau_Erho(double b) {
  const double Phi = 0.5 * std::erfc(-b * INV_SQRT2);
  const double phi = std::exp(-0.5 * b * b) / SQRT_2PI;
  return 2.0 * ((1.0 - b * b) * Phi - b * phi + b * b) - 1.0;
}
const double TAU_ES2 = tau_Erho(TAU_C2 * QNORM_3_4);

// scaleTau2 with robustbase's defaults (c1 = 4.5, c2 = 3, consistency at the
// normal, iter = 1). Returns the scale, and the location too when asked.
// `buf` is scratch the caller owns; we reorder and overwrite it.
double scale_tau2(const double* x, Eigen::Index n, std::vector<double>& buf, double* mu_out = nullptr) {
  buf.assign(x, x + n);
  const double mu0 = median_inplace(buf);
  for (Eigen::Index i = 0; i < n; ++i) buf[i] = std::abs(x[i] - mu0);
  const double sigma0 = median_inplace(buf);  // MAD without the consistency factor
  if (!(sigma0 > 0.0)) {                      // more than half the column is a single value
    if (mu_out) *mu_out = mu0;
    return 0.0;
  }

  // location: a one-step weighted mean, weights falling to zero at c1 * sigma0
  const double inv = 1.0 / (sigma0 * TAU_C1);
  double sw = 0.0, sxw = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    const double t = std::abs(x[i] - mu0) * inv;
    const double u = 1.0 - t * t;
    if (u > 0.0) {
      const double w = u * u;
      sw += w;
      sxw += x[i] * w;
    }
  }
  const double mu = (sw > 0.0) ? sxw / sw : mu0;

  // scale: the winsorized second moment about that location
  const double c2sq = TAU_C2 * TAU_C2;
  double rho = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    const double z = (x[i] - mu) / sigma0;
    rho += std::min(z * z, c2sq);
  }
  if (mu_out) *mu_out = mu;
  return sigma0 * std::sqrt(rho / ((double)n * TAU_ES2));
}

// One Maronna-Zamar step, in place: scale every column to unit tau scale, form
// the pairwise GK matrix of the scaled columns, and rotate onto its
// eigenvectors. `tmp` is scratch of the same shape as Z.
void ogk_step(Mat2D& Z, Mat2D& tmp) {
  const int p = (int)Z.cols();
  const Eigen::Index n = Z.rows();

#pragma omp parallel
  {
    std::vector<double> buf;
#pragma omp for schedule(static)
    for (int j = 0; j < p; ++j) {
      const double s = scale_tau2(Z.col(j).data(), n, buf);
      if (s > 0.0) Z.col(j) /= s;
    }
  }

  // the diagonal is 1 by construction, the columns now having unit tau scale
  Mat2D S = Mat2D::Identity(p, p);
  Int1D pi, pj;
  pi.reserve(p * (p - 1) / 2);
  pj.reserve(p * (p - 1) / 2);
  for (int i = 0; i < p; ++i)
    for (int j = i + 1; j < p; ++j) pi.push_back(i), pj.push_back(j);
  const int npairs = (int)pi.size();

#pragma omp parallel
  {
    std::vector<double> buf;
    Mat1D plus(n), minus(n);
#pragma omp for schedule(dynamic)
    for (int t = 0; t < npairs; ++t) {
      const int i = pi[t], j = pj[t];
      plus = Z.col(i) + Z.col(j);
      minus = Z.col(i) - Z.col(j);
      const double sp = scale_tau2(plus.data(), n, buf);
      const double sm = scale_tau2(minus.data(), n, buf);
      const double cij = 0.25 * (sp * sp - sm * sm);
      S(i, j) = cij;
      S(j, i) = cij;
    }
  }

  Eigen::SelfAdjointEigenSolver<Mat2D> eig(S);
  if (eig.info() != Eigen::Success) cao.error("failed eigendecomposition of the pcadapt OGK matrix.");
  tmp.noalias() = Z * eig.eigenvectors();
  Z.swap(tmp);
}

// Robust location and scatter by OGK, followed by the hard-rejection
// reweighting that bigutilsr::covrob_ogk() applies.
void robust_cov_ogk(const Mat2D& U, Mat1D& wcenter, Mat2D& wcov, int niter, double beta) {
  const int p = (int)U.cols();
  const Eigen::Index n = U.rows();

  Mat2D Z = U, tmp(n, p);
  for (int it = 0; it < niter; ++it) ogk_step(Z, tmp);
  tmp.resize(0, 0);

  // location and scale of the orthogonalized data; the distance is then just
  // the standardized sum of squares in that basis
  Mat1D ctr(p), sg(p);
#pragma omp parallel
  {
    std::vector<double> buf;
#pragma omp for schedule(static)
    for (int j = 0; j < p; ++j) {
      double m = 0.0;
      sg(j) = scale_tau2(Z.col(j).data(), n, buf, &m);
      ctr(j) = m;
    }
  }

  Mat1D d(n);
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j) {
      if (!(sg(j) > 0.0)) continue;  // degenerate direction carries no information
      const double t = (Z(i, j) - ctr(j)) / sg(j);
      acc += t * t;
    }
    d(i) = acc;
  }
  Z.resize(0, 0);

  std::vector<double> dv(d.data(), d.data() + n);
  const double cdelta = median_inplace(dv) / qchisq(0.5, p);
  const double cutoff = qchisq(beta, p) * cdelta;

  Eigen::Index nkeep = 0;
  for (Eigen::Index i = 0; i < n; ++i)
    if (d(i) < cutoff) ++nkeep;
  const bool use_all = (nkeep <= (Eigen::Index)p);
  if (use_all)
    cao.warn("pcadapt: the OGK reweighting kept only ", nkeep,
             " sites, too few to estimate a covariance; using all sites instead");

  auto kept = [&](Eigen::Index i) { return use_all || d(i) < cutoff; };
  const double nk = (double)(use_all ? n : nkeep);

  wcenter = Mat1D::Zero(p);
  for (Eigen::Index i = 0; i < n; ++i)
    if (kept(i)) wcenter += U.row(i).transpose();
  wcenter /= nk;

  wcov = Mat2D::Zero(p, p);
#pragma omp parallel
  {
    Mat2D acc = Mat2D::Zero(p, p);
    Mat1D v(p);
#pragma omp for schedule(static)
    for (Eigen::Index i = 0; i < n; ++i) {
      if (!kept(i)) continue;
      v = U.row(i).transpose() - wcenter;
      acc.selfadjointView<Eigen::Lower>().rankUpdate(v);
    }
#pragma omp critical
    wcov += acc;
  }
  for (int a = 0; a < p; ++a)  // rankUpdate only filled the lower triangle
    for (int b = a + 1; b < p; ++b) wcov(a, b) = wcov(b, a);
  wcov /= nk;
}

// `keep[i] == 0` marks a site with no residual variance: it is left out of the
// robust fit and the inflation factor, and its outputs come back as NaN.
void pcadapt_selection_stats(
    const Mat2D& Z, const std::vector<char>& keep, Mat1D& stat, Mat1D& chi2_stat, Mat1D& pval, double& gif) {
  const int k = (int)Z.cols();
  const Eigen::Index n = Z.rows();

  // Sites with no residual variance have an undefined z-score. They must not
  // reach the robust covariance, the inflation median or the output: pcadapt
  // fits on `zscores[pass, ]` and leaves the rest NA, and so do we. The common
  // case drops nothing, and then no copy is made.
  std::vector<Eigen::Index> idx;
  for (Eigen::Index i = 0; i < n; ++i)
    if (keep[i]) idx.push_back(i);
  const Eigen::Index nk = (Eigen::Index)idx.size();
  if (nk <= k)
    cao.error("pcadapt: only ", nk, " of ", n, " sites have any residual variance, too few to fit a ", k,
              "-dimensional covariance");

  const bool subset = (nk != n);
  Mat2D Zs;
  if (subset) {
    Zs.resize(nk, k);
    for (Eigen::Index i = 0; i < nk; ++i) Zs.row(i) = Z.row(idx[i]);
  }
  const Mat2D& W = subset ? Zs : Z;

  Mat1D center;
  Mat2D cov;
  robust_cov_ogk(W, center, cov, 2, 0.9);  // bigutilsr::covrob_ogk defaults
  Mat2D inv_cov = cov.inverse();
  if (!inv_cov.allFinite()) cao.error("pcadapt: the robust covariance of the z-scores is singular; try a smaller -k");

  Mat1D s(nk);
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < nk; ++i) {
    Mat1D dz = W.row(i).transpose() - center;
    s(i) = dz.transpose() * inv_cov * dz;
  }

  std::vector<double> stat_vec(s.data(), s.data() + nk);
  gif = median_inplace(stat_vec) / qchisq(0.5, k);  // stat_vec is ours to reorder
  if (!(gif > 0.0) || !std::isfinite(gif)) gif = 1.0;

  const double NA = std::numeric_limits<double>::quiet_NaN();
  stat = Mat1D::Constant(n, NA);
  chi2_stat = Mat1D::Constant(n, NA);
  pval = Mat1D::Constant(n, NA);
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < nk; ++i) {
    const Eigen::Index r = idx[i];
    stat(r) = s(i);
    chi2_stat(r) = s(i) / gif;
    pval(r) = pchisq(chi2_stat(r), k, false);
  }
}
}  // namespace

void run_selection(Data* data, const Param& params) {
  cao.print(tick.date(), "run selection");
  data->prepare();
  cao.print(tick.date(), "parsing U:", params.fileU, ", E:", params.fileE);
  Mat1D E = read_eigvals(params.fileE);
  const int K = ref_pcs(params, E.size(), params.fileE);  // -k: the leading PCs, else all
  cao.print(tick.date(), "selection: using", K, "of the", E.size(), "PCs in the reference");
  Mat2D U = read_eigvecs(params.fileU, data->nsamples, K);
  Mat2D V(data->nsnps, K);
  Mat1D y_norm2(data->nsnps);
  uint j;

  // Selection reads only .eigvals and .eigvecs, so the reference's transform --
  // which lives in .sigvals -- has to be fetched separately. Without it G was
  // scaled from this run's -C/--scale and silently mismatched the reference.
  UsvTransform usv;
  if (!params.fileS.empty()) {
    uint rn, rm;
    Mat1D rs;
    read_sigvals(params.fileS, rn, rm, rs, &usv);
  }
  const bool standardize = data->resolve_ref_scaling(usv, params.fileS.empty() ? "the reference PCA" : params.fileS);

  // V = G'U as one product in column panels of G; it was M matrix-vector
  // products, each writing a strided row of V
  if (!params.out_of_core) {
    if (standardize) data->standardize_E_ref(usv);
    mul_Xt_Y(data->G, U, V);
#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < data->nsnps; j++) y_norm2(j) = data->G.col(j).squaredNorm();
  } else {
    data->check_file_offset_first_var();
    for (uint b = 0; b < data->nblocks; b++) {
      // read unstandardized and scale here, so the factor follows the
      // reference rather than this run's --scale, which is what the reader
      // would consult
      data->read_block_initial(data->start[b], data->stop[b], false);
      uint64 actual_block_size = data->stop[b] - data->start[b] + 1;
      if (standardize) data->standardize_block_ref(usv, data->start[b], actual_block_size);
      mul_Xt_Y(data->G.leftCols(actual_block_size), U, V.middleRows(data->start[b], actual_block_size));
#pragma omp parallel for private(j) schedule(static)
      for (j = 0; j < actual_block_size; j++) y_norm2(j + data->start[b]) = data->G.col(j).squaredNorm();
    }
  }

  // V holds the regression coefficients U' g_j of each site on the orthonormal
  // PCs. Galinsky turns them into normalized loadings; pcadapt keeps the
  // coefficients themselves, and needs ||U' g||^2 to form the residual sum of
  // squares, so the division below belongs inside the Galinsky branch only.
  // Dividing first made the pcadapt residual variance ~ ||g||^2 for every site,
  // which deflates z exactly where the PCs explain the most -- the sites the
  // scan is for -- by a per-site factor no genomic-inflation step can undo.
  const double NA = std::numeric_limits<double>::quiet_NaN();
  const std::vector<Eigen::Index> rows = all_sites(data);
  if (!data->keepSNPs.empty())
    cao.warn(std::to_string(data->nsnps_all - data->nsnps) +
             " sites removed by --maf are written as NA rows, so every output keeps one row per site of the input");
  if (params.selection == 1) {
    // a new vector: E = E.head(K) * m resized E before reading it, so with -k
    // below the reference's PC count every statistic came from freed memory
    const Mat1D s2 = E.head(K) * (double)V.rows();         // squared singular values
    V.array().rowwise() /= s2.transpose().array().sqrt();  // divide by singular values
    cao.print(tick.date(), "calculate galinksky statistics");
    galinsky_selection_stat(V);
    uint64 ndrop = 0;
    for (j = 0; j < data->nsnps; ++j)
      if (!(y_norm2(j) > VAR_TOL)) V.row(j).setConstant(NA), ++ndrop;
    warn_dropped_sites(ndrop, data->nsnps, "variance");
    write_sites(params.fileout + ".galinsky", "#FastPCA/Galinsky selection statistic for each site and PC\n", V, rows);
    const Mat2D P = V.unaryExpr([](double x) { return pchisq(x, 1, false); });  // NaN in, NaN out
    write_sites(params.fileout + ".galinsky.pval",
                "#P-value for the FastPCA/Galinsky selection statistic for each site and PC\n", P, rows);
  } else if (params.selection == 2) {
    cao.print(tick.date(), "calculate pcadapt statistics");
    const int dof = static_cast<int>(data->nsamples) - K;
    if (dof <= 0) cao.error("pcadapt selection requires nsamples > K.");

    Mat2D Z = V;
    std::vector<char> keep(data->nsnps, 1);
#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < data->nsnps; ++j) {
      double rss = y_norm2(j) - V.row(j).squaredNorm();
      rss = std::max(rss, 0.0);
      double sigma = std::sqrt(rss / dof);
      if (sigma > VAR_TOL) {
        Z.row(j) /= sigma;
      } else {
        Z.row(j).setConstant(NA);  // no residual variance: z is not defined here
        keep[j] = 0;
      }
    }
    uint64 ndrop = 0;
    for (j = 0; j < data->nsnps; ++j) ndrop += (keep[j] ? 0 : 1);
    warn_dropped_sites(ndrop, data->nsnps, "residual variance");

    Mat1D stat, chi2_stat, pval;
    double gif = 1.0;
    pcadapt_selection_stats(Z, keep, stat, chi2_stat, pval, gif);

    write_sites(params.fileout + ".zscore", "#pcadapt z-scores for each site and PC\n", Z, rows);
    write_sites(params.fileout + ".pcadapt", "#pcadapt raw squared Mahalanobis statistic for each site\n", stat, rows);
    write_sites(params.fileout + ".pcadapt.chi2",
                "#pcadapt chi-square statistic after genomic inflation correction for each site\n", chi2_stat, rows);
    write_sites(params.fileout + ".pcadapt.pval", "#pcadapt p-value for each site\n", pval, rows);
    std::ofstream outg(params.fileout + ".pcadapt.gif");
    if (outg.is_open()) outg << "#genomic inflation factor\n" << gif << '\n';
  }
}
