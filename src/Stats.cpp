/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Stats.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Stats.hpp"

#include <limits>

#include "kfunc.h"
#include "Utils.hpp"

double chisq1d(const double x) {
  double p = kf_gammaq(1.0 / 2.0, x / 2.0);  // nan expected
  return std::isnan(p) ? 1.0 : p;            // if nan, then retrun 1.0
}

double pchisq(double x, int df, bool lower_tail) {
  // Validate degrees of freedom
  if (df <= 0) return std::numeric_limits<double>::quiet_NaN();

  // Handle edge cases
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (x <= 0.0) return lower_tail ? 0.0 : 1.0;
  if (std::isinf(x)) return lower_tail ? 1.0 : 0.0;

  double s = df / 2.0;
  double z = x / 2.0;

  // Use the complementary function for numerical stability:
  // - When p is close to 1 (large x), kf_gammap loses precision
  // - When p is close to 0 (small x), kf_gammaq loses precision
  // Switch based on which tail is smaller to get best precision,
  // then flip if needed.
  double result;
  bool use_lower = (x < df + 1.0);  // heuristic: switch near mode

  if (use_lower) {
    result = kf_gammap(s, z);  // lower incomplete gamma
    return lower_tail ? result : 1.0 - result;
  } else {
    result = kf_gammaq(s, z);  // upper incomplete gamma
    return lower_tail ? 1.0 - result : result;
  }
}

double qchisq(double p, int df) {
  if (df <= 0) return std::numeric_limits<double>::quiet_NaN();  // ADD THIS
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return std::numeric_limits<double>::infinity();

  double s = df / 2.0;

  // Initial guess using Wilson-Hilferty approximation
  // For chi-sq(df), approximate quantile:
  //   x ≈ df * (1 - 2/(9*df) + z_p * sqrt(2/(9*df)))^3
  // where z_p is the standard normal quantile of p.
  // Approximate z_p using rational approximation (Abramowitz & Stegun 26.2.23)
  double t;
  if (p < 0.5) {
    t = std::sqrt(-2.0 * std::log(p));
    t = t - (2.515517 + t * (0.802853 + t * 0.010328)) / (1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308)));
    t = -t;  // negative side
  } else {
    t = std::sqrt(-2.0 * std::log(1.0 - p));
    t = t - (2.515517 + t * (0.802853 + t * 0.010328)) / (1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308)));
  }

  double a = 2.0 / (9.0 * df);
  double x = df * std::pow(1.0 - a + t * std::sqrt(a), 3.0);
  if (x <= 0.0) x = 0.01;  // fallback

  // Log of chi-sq PDF normalization: log(2^(df/2) * Gamma(df/2))
  double log_norm = s * std::log(2.0) + std::lgamma(s);

  // Newton-Raphson iteration
  for (int iter = 0; iter < 100; ++iter) {
    double cdf = kf_gammap(s, x / 2.0);
    double err = cdf - p;

    // chi-sq PDF at x: f(x) = x^(s-1) * exp(-x/2) / (2^s * Gamma(s))
    double log_pdf = (s - 1.0) * std::log(x) - x / 2.0 - log_norm;
    double pdf = std::exp(log_pdf);

    // acceptable if converged, but should check |err| first
    if (pdf < 1e-300) break;  // avoid division by zero

    double delta = err / pdf;
    delta = std::max(delta, -x * 0.9);  // don't step to negative
    x -= delta;

    if (x <= 0.0) x = 1e-10;  // keep positive

    if (std::fabs(delta) < 1e-12 * (1.0 + x)) break;  // +1 for small x
  }

  return x;
}

void emMAF_with_GL(Mat1D& F, const Mat2D& P, int maxiter, double tolmaf) {
  uint nsnps = P.cols();
  uint nsamples = P.rows() / 2;
  Mat1D Ft = Mat1D::Zero(nsnps);
  double scale = 1.0 / (2.0 * nsamples);
  double diff;
  // run EM to estimate allele frequencies
  for (int it = 0; it < maxiter; it++) {
#pragma omp parallel for
    for (uint j = 0; j < nsnps; j++) {
      Ft(j) = F(j);
      double p0, p1, p2, pt = 0.0;
      for (uint i = 0; i < nsamples; i++) {
        p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
        p1 = P(2 * i + 1, j) * 2.0 * F(j) * (1.0 - F(j));
        p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
        pt += (p1 + 2.0 * p2) / (p0 + p1 + p2);
      }
      F(j) = pt * scale;
    }
    // calculate differences between iterations
    diff = sqrt((F - Ft).array().square().sum() / nsnps);
    // Check for convergence
    if (diff < tolmaf) {
      cao.print(tick.date(), "EM (MAF) converged at iteration:", it + 1);
      break;
    } else if (it == (maxiter - 1)) {
      cao.print(tick.date(), "EM (MAF) did not converge");
    }
  }
}
