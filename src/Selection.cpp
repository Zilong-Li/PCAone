#include "Selection.hpp"

#include "Utils.hpp"

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

  if (!params.out_of_core) {
    data->standardize_E();
#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < data->nsnps; j++) {
      V.row(j) = U.transpose() * data->G.col(j);
      y_norm2(j) = data->G.col(j).squaredNorm();
    }
  } else {
    data->check_file_offset_first_var();
    for (uint b = 0; b < data->nblocks; b++) {
      data->read_block_initial(data->start[b], data->stop[b], true);
      uint64 actual_block_size = data->stop[b] - data->start[b] + 1;
#pragma omp parallel for private(j) schedule(static)
      for (j = 0; j < actual_block_size; j++) {
        V.row(j + data->start[b]) = U.transpose() * data->G.col(j);
        y_norm2(j + data->start[b]) = data->G.col(j).squaredNorm();
      }
    }
  }

  E = E.head(K) * V.rows();                             // downscale
  V.array().rowwise() /= E.transpose().array().sqrt();  // divid by singluar values
  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  if (params.selection == 1) {
    cao.print(tick.date(), "calculate galinksky statistics");
    std::ofstream out(params.fileout + ".galinsky");
    galinsky_selection_stat(V);
    out << "#FastPCA/Galinsky selection statistic for each site and PC\n";
    out << V.format(fmt) << '\n';
    std::ofstream outp(params.fileout + ".galinsky.pval");
    if (outp.is_open()) {
      Mat2D P = V.unaryExpr([](double x) { return pchisq(x, 1, false); });
      outp << "#P-value for the FastPCA/Galinsky selection statistic for each site and PC\n";
      outp << P.format(fmt) << '\n';
    }
  } else if (params.selection == 2) {
    cao.print(tick.date(), "calculate pcadapt statistics");
    const int dof = static_cast<int>(data->nsamples) - K;
    if (dof <= 0) cao.error("pcadapt selection requires nsamples > K.");

    Mat2D Z = V;
#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < data->nsnps; ++j) {
      double rss = y_norm2(j) - V.row(j).squaredNorm();
      rss = std::max(rss, 0.0);
      double sigma = std::sqrt(rss / dof);
      if (sigma > VAR_TOL) {
        Z.row(j) /= sigma;
      } else {
        Z.row(j).setZero();
      }
    }

    Mat1D stat, chi2_stat, pval;
    double gif = 1.0;
    pcadapt_selection_stats(Z, stat, chi2_stat, pval, gif);

    std::ofstream outz(params.fileout + ".zscore");
    std::ofstream out(params.fileout + ".pcadapt");
    std::ofstream outc(params.fileout + ".pcadapt.chi2");
    std::ofstream outp(params.fileout + ".pcadapt.pval");
    std::ofstream outg(params.fileout + ".pcadapt.gif");

    if (outz.is_open()) {
      outz << "#pcadapt z-scores for each site and PC\n";
      outz << Z.format(fmt) << '\n';
    }
    if (out.is_open()) {
      out << "#pcadapt raw squared Mahalanobis statistic for each site\n";
      out << stat.format(fmt) << '\n';
    }
    if (outc.is_open()) {
      outc << "#pcadapt chi-square statistic after genomic inflation correction for each site\n";
      outc << chi2_stat.format(fmt) << '\n';
    }
    if (outp.is_open()) {
      outp << "#pcadapt p-value for each site\n";
      outp << pval.format(fmt) << '\n';
    }
    if (outg.is_open()) {
      outg << "#genomic inflation factor\n";      
      outg << gif << '\n';
    }
  }
}
