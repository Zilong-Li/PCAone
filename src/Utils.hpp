#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include <sys/utsname.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <stdexcept>

#include "Common.hpp"
#include "Logger.hpp"
#include "Timer.hpp"
#include "zstd.h"

// MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
Logger cao;  // logger
Timer tick;  // Timer
#else
extern Timer tick;
extern Logger cao;
#endif

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

struct SNPld {
  std::vector<int> pos;          // pos of each SNP
  std::vector<int> end_pos;      // 0-based index for last snp pos
  std::vector<std::string> chr;  // chr sequences
  Double1D af;                   // allele frequency
  std::vector<int> ws;           //  the snp index, i.e the index for lead SNP
  std::vector<int> we;  // the number of SNPs (including lead SNP) in a window
};

std::string get_machine();

void fcloseOrDie(FILE* file);

FILE* fopenOrDie(const char* filename, const char* instruction);

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file);

size_t count_lines(const std::string& fpath);

std::string timestamp();

void flip_UV(MyMatrix& U, MyMatrix& V, bool ubase = true);

void flip_Omg(MyMatrix& Omg2, MyMatrix& Omg);

void flip_Y(const MyMatrix& X, MyMatrix& Y);

double rmse(const MyMatrix& X, const MyMatrix& Y);

Eigen::VectorXd minSSE(const MyMatrix& X, const MyMatrix& Y);

double mev(const MyMatrix& X, const MyMatrix& Y);

void mev_rmse_byk(const MyMatrix& X, const MyMatrix& Y, MyVector& Vm,
                  MyVector& Vr);

double get_median(std::vector<double> v);

std::vector<std::string> split_string(const std::string& s,
                                      const std::string& separators);

void make_plink2_eigenvec_file(int K, std::string fout, const std::string& fin,
                               const std::string& fam);

struct ZstdBuffer {
  ZstdBuffer() {
    buffInTmp.reserve(buffInSize);
    buffOutTmp.reserve(buffOutSize);
  }
  ~ZstdBuffer() {
    ZSTD_freeDCtx(dctx);
    fcloseOrDie(fin);
  }
  FILE* fin = nullptr;
  size_t const buffInSize = ZSTD_DStreamInSize();
  size_t const buffOutSize = ZSTD_DStreamOutSize();
  ZSTD_DCtx* const dctx = ZSTD_createDCtx();
  size_t lastRet = 1;
  std::string buffCur = "";
  std::string buffLine, buffInTmp, buffOutTmp;
};

#endif  // PCAONE_UTILES_
