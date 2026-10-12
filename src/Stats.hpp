/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Stats.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_STATS_
#define PCAONE_STATS_

#include <stdexcept>
#include <type_traits>

#include "Common.hpp"

/// the p-value of a 1-degree chi-squared statistic; 1 for nan
double chisq1d(const double x);

// Chi-square CDF: P(X <= x) where X ~ chi-sq(df)
// chi-sq CDF = regularized lower incomplete gamma P(df/2, x/2)
// pchisq(x, df, lower_tail=true)  = kf_gammap(df/2, x/2)
// pchisq(x, df, lower_tail=false) = kf_gammaq(df/2, x/2)
double pchisq(double x, int df, bool lower_tail = true);

// Chi-square quantile function via Newton's method
// Find x such that P(X <= x) = p, where X ~ chi-sq(df)
// Uses the chi-square PDF for the Newton update:
//   f(x) = x^(df/2-1) * exp(-x/2) / (2^(df/2) * Gamma(df/2))
double qchisq(double p, int df);

/// estimate population allele frequency with genotype likelihoods
void emMAF_with_GL(Mat1D& F, const Mat2D& P, int maxiter, double tolmaf);

// Median of a buffer the caller is done with: reorders it in place.
//
// std::nth_element is O(n) where std::sort is O(n log n), and the value is
// identical -- after nth_element every element below `mid` is <= v[mid], so the
// largest of them is the other central order statistic. scale_tau2() (Selection.cpp)
// calls this millions of times, so the difference is the runtime of --selection 2.
template <typename T>
auto median_inplace(std::vector<T>& v) {
  static_assert(!std::is_same_v<T, bool>, "Boolean type is not supported");

  if (v.empty()) {
    throw std::invalid_argument("Cannot calculate median of an empty vector");
  }

  const size_t n = v.size();
  const size_t mid = n / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  if (n % 2 != 0) return v[mid];
  const T lo = *std::max_element(v.begin(), v.begin() + mid);
  return static_cast<T>((lo + v[mid]) / static_cast<T>(2));
}

#endif  // PCAONE_STATS_
