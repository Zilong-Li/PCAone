#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include <sys/utsname.h>

#include <Eigen/Dense>
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
#include <string>
#include <unordered_map>
#include <vector>

#include "Logger.hpp"
#include "Timer.hpp"

// MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
Logger cao;  // logger
Timer tick;  // Timer
#else
extern Timer tick;
extern Logger cao;
#endif

typedef Eigen::MatrixXd MyMatrix;
typedef Eigen::VectorXd MyVector;
typedef Eigen::ArrayXXd MyArrayX;
typedef Eigen::ArrayXd MyArray;
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

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

#define MAF(a) ((a) > 0.5 ? (1 - a) : (a))

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
  std::vector<int> ws;           //  the snp index, i.e the index for lead SNP
  std::vector<int> we;  // the number of SNPs (including lead SNP) in a window
  Double1D af;          // allele frequency
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

inline UMapIntInt vector2map(const Int1D& v) {
  UMapIntInt m;
  int i = 0;
  for (const auto& k : v) m.insert({k, i++});
  return m;
}

// get the chr, bp, id
struct BIM {
  std::string data;
  operator std::string const &() const { return data; }
  friend std::istream& operator>>(std::istream& ifs, BIM& BIM) {
    std::getline(ifs, BIM.data);
    auto token = split_string(BIM.data, " \t");
    BIM.data = token[0] + "\t" + token[3] + "\t" + token[1];
    return ifs;
  }
};

#endif  // PCAONE_UTILES_
