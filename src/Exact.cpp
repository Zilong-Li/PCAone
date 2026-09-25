/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Exact.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Exact.hpp"

#include <omp.h>

#include <Spectra/SymEigsSolver.h>

#include "Utils.hpp"

namespace {

// y = K x for the dense symmetric K, one column panel of K per thread
// (Spectra's DenseSymMatProd runs Eigen's symv on one thread)
class ParallelSymProd {
 public:
  using Scalar = double;
  explicit ParallelSymProd(const Mat2D& K) : K_(K) {}
  Eigen::Index rows() const { return K_.rows(); }
  Eigen::Index cols() const { return K_.cols(); }
  void perform_op(const double* x_in, double* y_out) const {
    const Eigen::Index n = K_.rows();
    Eigen::Map<const Mat1D> x(x_in, n);
    Eigen::Map<Mat1D> y(y_out, n);
    const Eigen::Index bs = std::max<Eigen::Index>(64, n / (4 * omp_get_max_threads()));
    const Eigen::Index nb = (n + bs - 1) / bs;
#pragma omp parallel for schedule(static)
    for (Eigen::Index b = 0; b < nb; ++b) {
      const Eigen::Index c = b * bs, w = std::min(bs, n - c);
      y.segment(c, w).noalias() = K_.middleCols(c, w).transpose() * x;  // K is symmetric
    }
  }

 private:
  const Mat2D& K_;
};

// Largest k eigenpairs of the symmetric K, in decreasing order. Lanczos
// (Spectra) to 1e-12 agrees with the full eigendecomposition to ~1e-13 and is
// far cheaper once N is large: 1 s against 50 s at N = 4000 on 2 threads.
void top_eigenpairs(const Mat2D& K, Eigen::Index k, Eigen::Index ncv, Mat1D& evals, Mat2D& U) {
  const Eigen::Index n = K.rows();
  ncv = std::min(std::max(ncv, 2 * k + 1), n);
  if (n > 1000 && ncv < n) {
    ParallelSymProd op(K);
    Spectra::SymEigsSolver<ParallelSymProd> eigs(op, k, ncv);
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestAlge, 1000, 1e-12);
    if (eigs.info() == Spectra::CompInfo::Successful) {
      evals = eigs.eigenvalues();
      U = eigs.eigenvectors();
      return;
    }
    cao.warn("the Lanczos solver did not converge on the GRM. falling back to the full eigendecomposition");
  }
  Eigen::SelfAdjointEigenSolver<Mat2D> eig(K);
  if (eig.info() != Eigen::Success) cao.error("failed eigendecomposition of the sample covariance matrix.");
  evals = eig.eigenvalues().tail(k).reverse();
  U = eig.eigenvectors().rightCols(k).rowwise().reverse();
}

// the sign rule of flip_UV(): the entry of largest magnitude of each PC is
// positive. V is formed from the flipped U, so it follows.
void flip_U(Mat2D& U) {
  for (Eigen::Index i = 0; i < U.cols(); ++i) {
    Eigen::Index x;
    U.col(i).cwiseAbs().maxCoeff(&x);
    if (U(x, i) < 0) U.col(i) *= -1;
  }
}

// K holds G G' in its lower triangle: mirror it and scale by 1/M
void finish_grm(Mat2D& K, double M) {
  mirror_lower(K);
  K /= M;
}

bool standardizes(const Param& params) {
  return params.file_t == FileType::PLINK || params.file_t == FileType::BGEN || params.file_t == FileType::PGEN;
}

// G is never held whole: each block is read (centred and standardized as
// in-core), added to the lower triangle of K = G G', and dropped. A second pass
// forms the loadings V = G' U / s only when -V asks for them.
void run_exact_streaming(Data* data, const Param& params) {
  const bool standardized = standardizes(params);
  const uint N = data->nsamples;
  const double M = data->nsnps;
  const Eigen::Index k = params.k;
  const double gb = (double)N * N * 8 / 1073741824;
  cao.print(tick.date(), "running exact PCA: streaming the", N, " x", N, " GRM of", gb, " GB over", data->nblocks,
            " blocks of", data->blocksize, " sites");
  Mat2D K;
  try {
    K = Mat2D::Zero(N, N);
  } catch (const std::bad_alloc&) {
    cao.error("out of memory allocating the", N, " x", N, " GRM of", gb, " GB for the exact PCA. use --svd 2 instead");
  }
  data->check_file_offset_first_var();
  for (uint b = 0; b < data->nblocks; ++b) {
    tick.clock();
    data->read_block_initial(data->start[b], data->stop[b], standardized);
    data->readtime += tick.reltime();
    const Eigen::Index bs = data->stop[b] - data->start[b] + 1;
    syrk_lower_add(K, data->G.leftCols(bs));
  }
  finish_grm(K, M);

  Mat1D evals;
  Mat2D U;
  top_eigenpairs(K, k, params.ncv, evals, U);
  K.resize(0, 0);  // free the GRM before the loadings pass
  evals = evals.cwiseMax(0.0);
  const Mat1D svals = (evals.array() * M).sqrt();
  flip_U(U);

  Mat2D V(data->nsnps, params.printv ? k : 0);  // write_eigs_files() takes nsnps from V.rows()
  if (params.printv) {
    cao.print(tick.date(), "second pass over the data for the loadings (-V)");
    data->check_file_offset_first_var();
    for (uint b = 0; b < data->nblocks; ++b) {
      tick.clock();
      data->read_block_initial(data->start[b], data->stop[b], standardized);
      data->readtime += tick.reltime();
      const Eigen::Index bs = data->stop[b] - data->start[b] + 1;
      mul_Xt_Y(data->G.leftCols(bs), U, V.middleRows(data->start[b], bs));
    }
    for (Eigen::Index i = 0; i < k; ++i)
      if (svals(i) > 0) V.col(i) /= svals(i);
  }
  data->set_svd_transform(standardized);
  data->write_eigs_files(evals, svals, U, V);
}

void run_exact_incore(Data* data, const Param& params) {
  const bool standardized = standardizes(params);
  if (standardized) data->standardize_E();
  cao.print(tick.date(), "running exact PCA with in-core eigendecomposition (PLINK-like).");
  const Eigen::Index ncomp = std::min<Eigen::Index>(params.k, std::min<Eigen::Index>(data->G.rows(), data->G.cols()));
  Mat1D evals(ncomp), svals(ncomp);
  Mat2D U(data->nsamples, ncomp), V(data->nsnps, ncomp);
  if (data->nsamples <= data->nsnps) {
    // the same GRM and solver as the streaming path, from the G in memory
    Mat2D K = Mat2D::Zero(data->nsamples, data->nsamples);
    syrk_lower_add(K, data->G);
    finish_grm(K, data->nsnps);
    top_eigenpairs(K, ncomp, params.ncv, evals, U);
    evals = evals.cwiseMax(0.0);
    svals = (evals.array() * data->nsnps).sqrt();
    mul_Xt_Y(data->G, U, V);
    for (Eigen::Index i = 0; i < ncomp; ++i) {
      if (svals(i) > 0) V.col(i) /= svals(i);
    }
  } else {
    Mat2D K = (data->G.transpose() * data->G) / data->nsnps;
    Eigen::SelfAdjointEigenSolver<Mat2D> eig(K);
    if (eig.info() != Eigen::Success) cao.error("failed eigendecomposition of the feature covariance matrix.");
    for (Eigen::Index i = 0; i < ncomp; ++i) {
      Eigen::Index idx = eig.eigenvalues().size() - 1 - i;
      evals(i) = std::max(0.0, eig.eigenvalues()(idx));
      V.col(i) = eig.eigenvectors().col(idx);
    }
    svals = (evals.array() * data->nsnps).sqrt();
    U.noalias() = data->G * V;
    for (Eigen::Index i = 0; i < ncomp; ++i) {
      if (svals(i) > 0) U.col(i) /= svals(i);
    }
  }
  flip_UV(U, V);
  data->set_svd_transform(standardized);
  data->write_eigs_files(evals, svals, U, V);
}

}  // namespace

void run_pca_exact(Data* data, const Param& params) {
  if (params.out_of_core)
    run_exact_streaming(data, params);
  else
    run_exact_incore(data, params);
}
