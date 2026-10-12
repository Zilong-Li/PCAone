/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/LinAlg.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 *
 * Dense kernels shared by the PCA methods: products threaded where Eigen's
 * GEMM is not, the sign convention of the singular vectors, and the measures
 * of convergence between iterations.
 ******************************************************************************/
#ifndef PCAONE_LINALG_
#define PCAONE_LINALG_

#include "Common.hpp"

// out = X' * Y for a wide X (samples x sites) and a thin Y, so out is tall.
// Eigen's multithreaded GEMM packs a kc x X.cols() panel of X' for a tall
// result, up to min(n, 256) * m doubles: nearly a second copy of X when n is
// small. Column panels of X, one thread each, need a few MB and run 2-4x
// faster. Every row of out comes from one panel, over the whole of n, so the
// result does not depend on the number of threads.
void mul_Xt_Y(const Eigen::Ref<const Mat2D>& X, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out);

// out = X * Y, or out += X * Y with add, for a tall X (samples x sites) and a
// thin Y, so out is samples x (k + oversamples). Eigen's multithreaded GEMM
// splits such a product over the columns of out only, at most cols / 4
// threads: 5 for -k 10 whatever -n is. Row panels of X, one per thread, use all
// of them and read X once. Each panel is a single-threaded product whose depth
// blocking does not depend on its height, and the panels start at multiples of
// every GEBP kernel height, so with 192 rows or more (two panels) the result is
// that of one thread, whatever -n is. Fewer rows go to Eigen's GEMM as before.
// With an external BLAS (MKL, OpenBLAS, Accelerate) the BLAS threads it instead.
void mul_X_Y(const Eigen::Ref<const Mat2D>& X,
             const Eigen::Ref<const Mat2D>& Y,
             Eigen::Ref<Mat2D> out,
             bool add = false);
// the number of row panels mul_X_Y uses for X with this many rows, here
Eigen::Index mul_X_Y_panels(Eigen::Index rows);

// K += X * X' on the lower triangle of the square K (the strict upper triangle
// is left untouched), for X of K.rows() x b. Eigen's rankUpdate runs on one
// thread and its GEMM does twice the flops, so the lower triangle is cut into
// tiles, each its own GEMM on one thread. A BLAS build hands it to ?syrk.
void syrk_lower_add(Eigen::Ref<Mat2D> K, const Eigen::Ref<const Mat2D>& X);

// copy the strict lower triangle of the square K onto its upper triangle
void mirror_lower(Eigen::Ref<Mat2D> K);

// Call f(j, x) for j in [0, count), where x is column `first + j` of A * B
// (A is n x k, B is k x m). The product is formed a panel of ~1 MB at a time in
// a thread-local buffer: one GEMM per panel in place of the scalar triple loops
// the reconstructions U * S * V' were written as. Runs in parallel over panels.
template <typename DerivedB, typename Func>
void for_each_product_column(
    const Mat2D& A, const Eigen::MatrixBase<DerivedB>& B, Eigen::Index first, Eigen::Index count, Func f) {
  const Eigen::Index bs = std::max<Eigen::Index>(1, (Eigen::Index(1) << 17) / std::max<Eigen::Index>(1, A.rows()));
  const Eigen::Index nb = (count + bs - 1) / bs;
#pragma omp parallel
  {
    Mat2D panel;
#pragma omp for schedule(static)
    for (Eigen::Index b = 0; b < nb; ++b) {
      const Eigen::Index c = b * bs, w = std::min<Eigen::Index>(bs, count - c);
      panel.noalias() = A * B.middleCols(first + c, w);
      for (Eigen::Index t = 0; t < w; ++t) f(c + t, panel.col(t));
    }
  }
}

// centre every column of X, and scale it to unit sd unless the sd is below tol
void standardize(Mat2D& X, double tol = 1e-10);

// Sign correction to make the SVD deterministic: the entry of largest
// magnitude of each column of U (of V with ubase = false) is made positive,
// and the matching singular vector flipped with it. V may be k x m or m x k.
void flip_UV(Mat2D& U, Mat2D& V, bool ubase = true);

// --- convergence between iterations ------------------------------------------

// mean explained variance: the mean norm of the projections of the columns of
// Y on the span of X; 1 when the two spans agree
double mev(const Mat2D& X, const Mat2D& Y);

// for each column of X, the least sum of squared differences to any column of
// Y, either sign. Y is the truth, X the test
Mat1D minSSE(const Mat2D& X, const Mat2D& Y);

double rmse1d(const Mat1D& x, const Mat1D& y);

#endif  // PCAONE_LINALG_
