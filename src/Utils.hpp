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

MyVector cor_cross(const MyMatrix & X, const MyVector & Y);

void cor_by_window(const std::string & filein,
                   const std::string & fileout,
                   MyMatrix & X,
                   const std::vector<int> & ws,
                   const std::vector<int> & we,
                   double r2_tol);

void calc_ld_metrics(const std::string & filein,
                     const std::string & fileout,
                     MyMatrix & G,
                     const MyMatrix & U,
                     const MyVector & S,
                     const MyMatrix & V,
                     const std::vector<int> & snp_pos,
                     const std::vector<int> & chr_pos_end,
                     int ld_window_bp,
                     double r2_tol);

std::vector<std::string> split_string(const std::string & s, const std::string & separators);

void get_snp_pos_bim(const std::string & filebim, std::vector<int> & pos, std::vector<int> & chr_pos_end);

#endif // PCAONE_UTILES_
