/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/LinAlg.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "LinAlg.hpp"

#include <omp.h>

#include <utility>

#include "Utils.hpp"

// The row panels of mul_X_Y: their height, a multiple of 48 (of every GEBP
// kernel height, so the kernels take the rows of the whole product in the same
// groups), and how many have rows. A last panel under 48 rows joins the one
// before: a one-row panel would be a GEMV, which sums in another order.
static std::pair<Eigen::Index, Eigen::Index> mul_X_Y_split(Eigen::Index n) {
#ifdef EIGEN_USE_BLAS
  return {n, 1};  // the BLAS threads the product itself
#else
  // panels of at least 96 rows, at most one per thread
  const Eigen::Index T = omp_in_parallel() ? 1 : omp_get_max_threads();
  const Eigen::Index panels = std::min<Eigen::Index>(T, n / 96);
  if (panels <= 1) return {n, 1};
  const Eigen::Index h = ((n + panels - 1) / panels + 47) / 48 * 48;
  Eigen::Index used = (n + h - 1) / h;
  if (used > 1 && n - (used - 1) * h < 48) --used;
  return {h, used};
#endif
}

Eigen::Index mul_X_Y_panels(Eigen::Index rows) { return mul_X_Y_split(rows).second; }

void mul_X_Y(const Eigen::Ref<const Mat2D>& X, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out, bool add) {
  const Eigen::Index n = X.rows();
  // plain locals: Apple Clang cannot capture structured bindings in an OpenMP loop
  const auto split = mul_X_Y_split(n);
  const Eigen::Index h = split.first, panels = split.second;
  if (panels <= 1) {
    if (add)
      out.noalias() += X * Y;
    else
      out.noalias() = X * Y;
    return;
  }
#pragma omp parallel for schedule(static)
  for (Eigen::Index p = 0; p < panels; ++p) {
    const Eigen::Index r = p * h, m = p + 1 < panels ? h : n - r;
    if (add)
      out.middleRows(r, m).noalias() += X.middleRows(r, m) * Y;
    else
      out.middleRows(r, m).noalias() = X.middleRows(r, m) * Y;
  }
}

void mul_Xt_Y(const Eigen::Ref<const Mat2D>& X, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out) {
  const Eigen::Index m = X.cols();
  // ~2 MB of X per panel, and at least 64 columns
  const Eigen::Index bs = std::max<Eigen::Index>(64, (Eigen::Index(1) << 18) / std::max<Eigen::Index>(1, X.rows()));
  const Eigen::Index nb = (m + bs - 1) / bs;
#pragma omp parallel for schedule(static)
  for (Eigen::Index b = 0; b < nb; ++b) {
    const Eigen::Index c = b * bs, w = std::min(bs, m - c);
    out.middleRows(c, w).noalias() = X.middleCols(c, w).transpose() * Y;
  }
}

void syrk_lower_add(Eigen::Ref<Mat2D> K, const Eigen::Ref<const Mat2D>& X) {
  const Eigen::Index n = K.rows(), m = X.cols();
  if (n == 0 || m == 0) return;
#ifdef EIGEN_USE_BLAS
  K.selfadjointView<Eigen::Lower>().rankUpdate(X);  // ?syrk, threaded by the BLAS
#else
  // Tiles of the lower triangle, 256 rows each (one tile up to n = 512). When
  // there are fewer than 2 tiles per thread, the sites are also split into D
  // chunks, each added into its own copy of K (at most 1 GB of copies) and
  // summed at the end: a few large products run faster than many small ones.
  const Eigen::Index T = omp_get_max_threads();
  const Eigen::Index ts = n <= 512 ? n : 256;
  const Eigen::Index nt = (n + ts - 1) / ts, ntiles = nt * (nt + 1) / 2;
  const Eigen::Index D = std::max<Eigen::Index>(
      1, std::min<Eigen::Index>({T, (2 * T + ntiles - 1) / ntiles, 1 + (Eigen::Index(1) << 27) / (n * n), m / 256}));
  // chunk 0 goes into K itself. not Kd(D - 1, Mat2D::Zero(n, n)): its fill
  // value is an n x n temporary even when D = 1, a second N x N matrix at the
  // peak of the dense --evaladmix and of the exact PCA
  std::vector<Mat2D> Kd(D - 1);
  for (auto& Kc : Kd) Kc.setZero(n, n);
  #pragma omp parallel for schedule(dynamic, 1)
  for (Eigen::Index task = 0; task < ntiles * D; ++task) {
    const Eigen::Index t = task % ntiles, d = task / ntiles;
    // t -> tile (i, j), j <= i, row by row of the triangle of tiles
    Eigen::Index i = (Eigen::Index)((std::sqrt(8.0 * t + 1) - 1) / 2);
    while (i * (i + 1) / 2 > t) --i;
    while ((i + 1) * (i + 2) / 2 <= t) ++i;
    const Eigen::Index j = t - i * (i + 1) / 2;
    const Eigen::Index r = i * ts, h = std::min(ts, n - r), c = j * ts, w = std::min(ts, n - c);
    const Eigen::Index c0 = m * d / D, c1 = m * (d + 1) / D;
    const auto Xd = X.middleCols(c0, c1 - c0);
    Eigen::Map<Mat2D, 0, Eigen::OuterStride<>> Kt(d == 0 ? K.data() : Kd[d - 1].data(), n, n,
                                                  Eigen::OuterStride<>(d == 0 ? K.outerStride() : n));
    if (i == j)
      Kt.block(r, r, h, h).selfadjointView<Eigen::Lower>().rankUpdate(Xd.middleRows(r, h));
    else
      Kt.block(r, c, h, w).noalias() += Xd.middleRows(r, h) * Xd.middleRows(c, w).transpose();
  }
  if (D > 1) {
  #pragma omp parallel for schedule(dynamic, 16)
    for (Eigen::Index j = 0; j < n; ++j)
      for (const auto& Q : Kd) K.col(j).tail(n - j) += Q.col(j).tail(n - j);
  }
#endif
}

void mirror_lower(Eigen::Ref<Mat2D> K) {
  const Eigen::Index n = K.rows();
#pragma omp parallel for schedule(dynamic, 64)
  for (Eigen::Index j = 1; j < n; ++j)
    for (Eigen::Index i = 0; i < j; ++i) K(i, j) = K(j, i);
}

void standardize(Mat2D& X, double tol) {
  double sqrt_rdf = sqrt(X.rows() - 1.0);
  // if X is centered, then we can convert norm to sd
  for (Eigen::Index j = 0; j < X.cols(); j++) {
    double mean = X.col(j).mean();
    X.col(j).array() -= mean;
    double sd = X.col(j).norm() / sqrt_rdf;  // sd
    if (sd > tol) X.col(j) /= sd;
  }
}

// Sign correction to ensure deterministic output from SVD.
// see https://www.kite.com/python/docs/sklearn.utils.extmath.svd_flip
void flip_UV(Mat2D& U, Mat2D& V, bool ubase) {
  if (ubase) {
    Eigen::Index x, i;
    for (i = 0; i < U.cols(); ++i) {
      U.col(i).cwiseAbs().maxCoeff(&x);
      if (U(x, i) < 0) {
        U.col(i) *= -1;
        if (V.cols() == U.cols()) {
          V.col(i) *= -1;
        } else if (V.rows() == U.cols()) {
          V.row(i) *= -1;
        } else {
          cao.error("the dimention of U and V have different k ranks.\n");
        }
      }
    }
  } else {
    Eigen::Index x, i;
    for (i = 0; i < V.cols(); ++i) {
      if (V.cols() == U.cols()) {
        V.col(i).cwiseAbs().maxCoeff(&x);
        if (V(x, i) < 0) {
          U.col(i) *= -1;
          V.col(i) *= -1;
        }
      } else if (V.rows() == U.cols()) {
        V.row(i).cwiseAbs().maxCoeff(&x);
        if (V(i, x) < 0) {
          U.col(i) *= -1;
          V.row(i) *= -1;
        }
      } else {
        cao.error("the dimention of U and V have different k ranks.\n");
      }
    }
  }
}

double mev(const Mat2D& X, const Mat2D& Y) {
  double res = 0;
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    res += (X.transpose() * Y.col(i)).norm();
  }
  return res / X.cols();
}

// Y is the truth matrix, X is the test matrix
Mat1D minSSE(const Mat2D& X, const Mat2D& Y) {
  Eigen::Index w1, w2;
  Mat1D res(X.cols());
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    // test against the original matrix to find the index with mincoeff
    ((-Y).colwise() + X.col(i)).array().square().colwise().sum().minCoeff(&w1);
    // test against the flipped matrix with the opposite sign
    (Y.colwise() + X.col(i)).array().square().colwise().sum().minCoeff(&w2);
    // get the minSSE value for X.col(i) against -Y.col(w1)
    auto val1 = (-Y.col(w1) + X.col(i)).array().square().sum();
    // get the minSSE value for X.col(i) against Y.col(w2)
    auto val2 = (Y.col(w2) + X.col(i)).array().square().sum();
    if (w1 != w2 && val1 > val2)
      res[i] = val2;
    else
      res[i] = val1;
  }
  return res;
}

double rmse1d(const Mat1D& x, const Mat1D& y) { return sqrt((x - y).array().square().sum() / x.size()); }
