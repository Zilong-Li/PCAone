#include "Selection.hpp"

#include "Utils.hpp"

namespace {
// Eigen's IOFormat prints a NaN as "nan"; write NA instead. That is what
// pcadapt reports for the sites it drops, and what both R and pandas read back
// as missing. Stream precision 6 is what the IOFormat this replaces used, so
// every finite value comes out byte for byte unchanged.
template <typename Derived>
void write_rows(std::ostream& os, const Eigen::DenseBase<Derived>& M) {
  const std::streamsize old = os.precision(6);
  for (Eigen::Index i = 0; i < M.rows(); ++i) {
    for (Eigen::Index j = 0; j < M.cols(); ++j) {
      if (j) os << '\t';
      const double v = M(i, j);
      if (std::isfinite(v))
        os << v;
      else
        os << "NA";
    }
    os << '\n';
  }
  os.precision(old);
}

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
}  // namespace

void run_selection(Data* data, const Param& params) {
  cao.print(tick.date(), "run selection");
  data->prepare();
  cao.print(tick.date(), "parsing U:", params.fileU, ", E:", params.fileE);
  Mat1D E = read_eigvals(params.fileE);
  int K = fmin(E.size(), params.k);
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
  const bool standardize =
      data->resolve_ref_scaling(usv, params.fileS.empty() ? "the reference PCA" : params.fileS);

  if (!params.out_of_core) {
    if (standardize) data->standardize_E_ref(usv);
#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < data->nsnps; j++) {
      V.row(j) = U.transpose() * data->G.col(j);
      y_norm2(j) = data->G.col(j).squaredNorm();
    }
  } else {
    data->check_file_offset_first_var();
    for (uint b = 0; b < data->nblocks; b++) {
      // read unstandardized and scale here, so the factor follows the
      // reference rather than this run's --scale, which is what the reader
      // would consult
      data->read_block_initial(data->start[b], data->stop[b], false);
      uint64 actual_block_size = data->stop[b] - data->start[b] + 1;
      if (standardize) data->standardize_block_ref(usv, data->start[b], actual_block_size);
#pragma omp parallel for private(j) schedule(static)
      for (j = 0; j < actual_block_size; j++) {
        V.row(j + data->start[b]) = U.transpose() * data->G.col(j);
        y_norm2(j + data->start[b]) = data->G.col(j).squaredNorm();
      }
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
  if (params.selection == 1) {
    E = E.head(K) * V.rows();                             // downscale
    V.array().rowwise() /= E.transpose().array().sqrt();  // divide by singular values
    cao.print(tick.date(), "calculate galinksky statistics");
    std::ofstream out(params.fileout + ".galinsky");
    galinsky_selection_stat(V);
    uint64 ndrop = 0;
    for (j = 0; j < data->nsnps; ++j)
      if (!(y_norm2(j) > VAR_TOL)) V.row(j).setConstant(NA), ++ndrop;
    warn_dropped_sites(ndrop, data->nsnps, "variance");
    out << "#FastPCA/Galinsky selection statistic for each site and PC\n";
    write_rows(out, V);
    std::ofstream outp(params.fileout + ".galinsky.pval");
    if (outp.is_open()) {
      Mat2D P = V.unaryExpr([](double x) { return pchisq(x, 1, false); });  // NaN in, NaN out
      outp << "#P-value for the FastPCA/Galinsky selection statistic for each site and PC\n";
      write_rows(outp, P);
    }
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

    std::ofstream outz(params.fileout + ".zscore");
    std::ofstream out(params.fileout + ".pcadapt");
    std::ofstream outc(params.fileout + ".pcadapt.chi2");
    std::ofstream outp(params.fileout + ".pcadapt.pval");
    std::ofstream outg(params.fileout + ".pcadapt.gif");

    if (outz.is_open()) {
      outz << "#pcadapt z-scores for each site and PC\n";
      write_rows(outz, Z);
    }
    if (out.is_open()) {
      out << "#pcadapt raw squared Mahalanobis statistic for each site\n";
      write_rows(out, stat);
    }
    if (outc.is_open()) {
      outc << "#pcadapt chi-square statistic after genomic inflation correction for each site\n";
      write_rows(outc, chi2_stat);
    }
    if (outp.is_open()) {
      outp << "#pcadapt p-value for each site\n";
      write_rows(outp, pval);
    }
    if (outg.is_open()) {
      outg << "#genomic inflation factor\n";
      outg << gif << '\n';
    }
  }
}
