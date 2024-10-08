/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Common.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAone_Common_H
#define PCAone_Common_H

#include <Eigen/Dense>
#include <algorithm>
#include <cstdio>
#include <iterator>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#define ZSTD_STATIC_LINKING_ONLY
#include "zstd.h"

#define MAF(a) ((a) > 0.5 ? (1 - a) : (a))

typedef Eigen::MatrixXd Mat2D;
typedef Eigen::VectorXd Mat1D;
typedef Eigen::ArrayXXd Arr2D;
typedef Eigen::ArrayXd Arr1D;
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrBool;
using PermMat = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
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
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
}

struct Line {
  std::string data;
  operator std::string const &() const { return data; }
  friend std::istream& operator>>(std::istream& ifs, Line& line) {
    return std::getline(ifs, line.data);
  }
};

// zstd deccompression buffer
struct ZstdDS {
  ZstdDS() {
    buffInTmp.reserve(buffInSize);
    buffOutTmp.reserve(buffOutSize);
  }
  ~ZstdDS() {
    ZSTD_freeDCtx(dctx);
    if (fclose(fin)) {
      perror("fclose error");
      exit(1);
    }
  }
  FILE* fin = nullptr;
  size_t const buffInSize = ZSTD_DStreamInSize();
  size_t const buffOutSize = ZSTD_DStreamOutSize();
  ZSTD_DCtx* const dctx = ZSTD_createDCtx();
  size_t lastRet = 1;
  std::string buffCur = "";
  std::string buffLine, buffInTmp, buffOutTmp;
};

// zstd deccompression buffer
struct ZstdCS {
  ZstdCS() {
    buffInTmp.reserve(buffInSize);
    buffOutTmp.reserve(buffOutSize);
  }
  ~ZstdCS() {
    ZSTD_freeCCtx(cctx);
    if (fclose(fout)) {
      perror("fclose error");
      exit(1);
    }
  }
  FILE* fout = nullptr;
  size_t const buffInSize = ZSTD_CStreamInSize();
  size_t const buffOutSize = ZSTD_CStreamOutSize();
  ZSTD_CCtx* const cctx = ZSTD_createCCtx();
  size_t lastRet = 1;
  std::string buffCur = "";
  std::string buffLine, buffInTmp, buffOutTmp;
};

#endif  // PCAone_Common_H
