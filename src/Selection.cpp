#include "Selection.hpp"

#include "Utils.hpp"

void run_selection(Data* data, const Param& params) {
  cao.print(tick.date(), "run selection");
  data->prepare();
  cao.print(tick.date(), "start parsing U:", params.fileU, ", S:", params.fileS);
  Mat2D U = read_eigvecs(params.fileU, data->nsamples, params.k);
  Mat2D V(data->nsnps, params.k);

  Eigen::ColPivHouseholderQR<Mat2D> qr(U);

  if (!params.out_of_core) {
    data->standardize_E();
#pragma omp parallel for
    for (uint j = 0; j < data->nsnps; j++) {
      V.row(j) = qr.solve(data->G.col(j));
    }
  } else {
    data->check_file_offset_first_var();
    for (uint b = 0; b < data->nblocks; b++) {
      data->read_block_initial(data->start[b], data->stop[b], true);
      uint64 actual_block_size = data->stop[b] - data->start[b] + 1;
#pragma omp parallel for
      for (uint j = 0; j < actual_block_size; j++) {
        V.row(j + data->start[b]) = qr.solve(data->G.col(j));
      }
    }
  }

  Eigen::IOFormat fmt(6, Eigen::DontAlignCols, "\t", "\n");
  std::ofstream out(params.fileout + ".zscore");
  if (out.is_open()) out << V.format(fmt) << '\n';
}
