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
#include <sys/utsname.h>

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

template<typename MatrixType>
inline void permute_matrix(MatrixType & G,
                           Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> & P,
                           bool bycol = true)
{
    if(bycol)
    {
        P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(G.cols());
        P.setIdentity();
        auto rng = std::default_random_engine{};
        std::shuffle(P.indices().data(), P.indices().data() + P.indices().size(), rng);
        G = G * P; // permute columns in-place
    }
    else
    {
        P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(G.rows());
        P.setIdentity();
        auto rng = std::default_random_engine{};
        std::shuffle(P.indices().data(), P.indices().data() + P.indices().size(), rng);
        G = P * G; // permute rows in-place
    }
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

#endif // PCAONE_UTILES_
