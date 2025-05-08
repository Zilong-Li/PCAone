/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Common.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAone_Common_H
#define PCAone_Common_H

#include <Eigen/Dense>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#define ZSTD_STATIC_LINKING_ONLY
#include "zstd.h"

#define MAF(a) ((a) > 0.5 ? (1 - a) : (a))

using Mat2D = Eigen::MatrixXd;
using Mat1D = Eigen::VectorXd;
using Arr2D = Eigen::ArrayXXd;
using Arr1D = Eigen::ArrayXd;
using ArrBool = Eigen::Array<bool, Eigen::Dynamic, 1>;
using PermMat = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;

using uchar = unsigned char;
using uint = unsigned int;
using uint64 = unsigned long long;
using Int1D = std::vector<int>;
using Int2D = std::vector<Int1D>;
using Double1D = std::vector<double>;
using String1D = std::vector<std::string>;
using UMapIntInt = std::unordered_map<int, int>;
using UMapIntDouble = std::unordered_map<int, double>;
using UMapIntString = std::unordered_map<int, std::string>;
using Pds = std::pair<double, std::string>;
using UMapIntPds = std::unordered_map<int, Pds>;


/**
 * Recode genotype codes to allelic dosages of first allele in .bim file (A1),
 * similarly to .raw files generated with '--recode A' in PLINK. A coding for
 * the missing value needs to be provided in 'na_value'.
 * 00 ->  2 (copies of A1)
 * 10 ->  1 (copy of A1)
 * 11 ->  0 (copy of A1)
 * 01 ->  3 (missing)
 */
const double BED_MISSING_VALUE = -9;
const double BED2GENO[4] = {1, BED_MISSING_VALUE, 0.5, 0};

inline UMapIntInt vector2map(const Int1D& v) {
  UMapIntInt m;
  int i = 0;
  for (const auto& k : v) m.insert({k, i++});
  return m;
}

template <typename T>
inline std::vector<size_t> sortidx(const std::vector<T>& v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
}

struct Line {
  std::string data;
  operator std::string const&() const { return data; }
  friend std::istream& operator>>(std::istream& ifs, Line& line) { return std::getline(ifs, line.data); }
};

// C++ equivalent of Python's math.isclose function
// https://github.com/jamadagni/areclose/blob/master/areclose.hpp
struct AreClose {
  AreClose() : abs_tol{0}, rel_tol{1e-09} {}

  AreClose(double absTol, double relTol) : abs_tol{absTol}, rel_tol{relTol} {
    if (invalidTolerance(abs_tol) || invalidTolerance(rel_tol)) abs_tol = rel_tol = NAN;
    // this makes operator() always return false
  }

  bool operator()(double a, double b) const {
    // The following is needed to catch infinities of the same sign since their
    // difference is NAN. And rarely, finite inputs may also be equal.
    if (a == b) return true;

    // Since an infinity would have an infinite relative tolerance,
    // any finite number would be considered relatively close to an infinity.
    // Further we need to catch infinities of opposite signs whose
    // difference is infinite and would compare equal to their relative
    // tolerance. The following is needed for both those cases:
    if (std::isinf(a) || std::isinf(b)) return false;

    // The below is effectively the same as:
    // abs(a - b) <= max(abs_tol, rel_tol * max(abs(a), abs(b)))
    double diff = std::fabs(a - b);
    return diff <= abs_tol || diff <= std::fabs(rel_tol * a) || diff <= std::fabs(rel_tol * b);
  }

  static bool invalidTolerance(double tol) {
    // A tolerance may be zero to switch off that particular type of tolerance checking
    // but if it is non-zero then it should be finite and at least 4 × machine epsilon
    // https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/roots_noderiv/root_termination.html
    return !(tol == 0.0 || (tol >= 4 * DBL_EPSILON && std::isfinite(tol)));
  }

 private:
  double abs_tol, rel_tol;
};

#endif  // PCAone_Common_H
