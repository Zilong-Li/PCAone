#include "Projection.hpp"

#include "Utils.hpp"

/**
 * options:
 * 1: simple, assume no missingness
 * 2: smartPCA, solving g=Vx, can take missing genotypes
 * 3: OADP, laser2, can take missing genotypes
 */
void run_projection(Data* data, const Param& params) {
  cao.print(tick.date(), "run projection");
  // get 1 /  Singular = sqrt(Eigen * M)
  Mat1D S =
      1.0 /
      (read_eigvals(params.fileS) * data->nsnps / params.ploidy).array().sqrt();

  int pcs = S.size();
  Mat2D V = read_eigvecs(params.fileV, data->nsnps, pcs);  // M x K

  // X V = U D
  data->standardize_E();
  Mat2D U = (data->G * V) * S.asDiagonal();
  data->write_eigs_files(1.0 / S.array(), U, V);
}
