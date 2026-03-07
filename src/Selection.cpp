#include "Selection.hpp"

#include "Utils.hpp"

void run_selection(Data* data, const Param& params) {
  cao.print(tick.date(), "run selection");
  data->prepare();
  cao.print(tick.date(), "start parsing U:", params.fileU, ", E:", params.fileE);
  Mat1D E = read_eigvals(params.fileE);
  int K = fmin(E.size(), params.k);
  Mat2D U = read_eigvecs(params.fileU, data->nsamples, K);
  Mat2D V(data->nsnps, K);

  if (!params.out_of_core) {
    data->standardize_E();
#pragma omp parallel for
    for (uint j = 0; j < data->nsnps; j++) {
      V.row(j) = U.transpose() * data->G.col(j);
    }
  } else {
    data->check_file_offset_first_var();
    for (uint b = 0; b < data->nblocks; b++) {
      data->read_block_initial(data->start[b], data->stop[b], true);
      uint64 actual_block_size = data->stop[b] - data->start[b] + 1;
#pragma omp parallel for
      for (uint j = 0; j < actual_block_size; j++) {
        V.row(j + data->start[b]) = U.transpose() * data->G.col(j);
      }
    }
  }

  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  std::ofstream outz(params.fileout + ".zscore");
  if (outz.is_open()) outz << V.format(fmt) << '\n';
  if (params.selection == 1) {
    std::ofstream out(params.fileout + ".galinsky");
    E = E.head(K) * V.rows();
    V.array().rowwise() /= E.transpose().array().sqrt();
    galinsky_selection_scan(V);
    out << "#Pvalue of each site for each PC with galinksky selection test\n";
    out << V.format(fmt) << '\n';
  }

}
