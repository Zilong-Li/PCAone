#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include "Logger.hpp"
#include "Timer.hpp"
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
#include <sys/utsname.h>
#include <unordered_map>
#include <vector>

// MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
Logger cao; // logger
Timer tick; // Timer
#else
extern Timer tick;
extern Logger cao;
#endif

typedef Eigen::MatrixXd MyMatrix;
typedef Eigen::VectorXd MyVector;
typedef Eigen::ArrayXXd MyArrayX;
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
using Int1D = std::vector<int>;
using Int2D = std::vector<Int1D>;
using UMapIntInt = std::unordered_map<int, int>;
using UMapIntDouble = std::unordered_map<int, double>;

#define MAF(a) ((a) > 0.5 ? (1 - a) : (a))

template<typename T>
inline std::vector<size_t> sortidx(const std::vector<T> & v)
{
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
    return idx;
}

inline std::vector<size_t> sortidx(const UMapIntDouble & m)
{
    std::vector<double> ps;
    for(auto it = m.begin(); it != m.end(); it++) ps.push_back(it->second);
    return sortidx(ps);
}

struct Line
{
    std::string data;
    operator std::string const &() const
    {
        return data;
    }
    friend std::istream & operator>>(std::istream & is, Line & line)
    {
        return std::getline(is, line.data);
    }
};

std::string get_machine();

void fcloseOrDie(FILE * file);

FILE * fopenOrDie(const char * filename, const char * instruction);

size_t freadOrDie(void * buffer, size_t sizeToRead, FILE * file);

size_t count_lines(const std::string & fpath);

std::string timestamp();

void flip_UV(MyMatrix & U, MyMatrix & V, bool ubase = true);

void flip_Omg(MyMatrix & Omg2, MyMatrix & Omg);

void flip_Y(const MyMatrix & X, MyMatrix & Y);

double rmse(const MyMatrix & X, const MyMatrix & Y);

Eigen::VectorXd minSSE(const MyMatrix & X, const MyMatrix & Y);

double mev(const MyMatrix & X, const MyMatrix & Y);

void mev_rmse_byk(const MyMatrix & X, const MyMatrix & Y, MyVector & Vm, MyVector & Vr);

double get_median(std::vector<double> v);

std::vector<std::string> split_string(const std::string & s, const std::string & separators);

void get_snp_pos_bim(const std::string & filebim,
                     Int1D & pos,
                     Int1D & chr_pos_end,
                     std::vector<std::string> & chrs,
                     bool header = false,
                     Int1D idx = Int1D{0, 3});

#endif // PCAONE_UTILES_
