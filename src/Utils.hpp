#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <iostream>
#include <iterator>


const std::string colwarn = "\x1B[33mwarning: ";
const std::string colerror = "\x1B[31merror: ";
const std::string coldone = "\x1B[32m";
const std::string colend = "\033[0m";

typedef Eigen::MatrixXd MyMatrix;
typedef Eigen::VectorXd MyVector;
typedef Eigen::ArrayXXd MyArrayX;
typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

struct Line
{
    std::string data;
    operator std::string const&() const
    {
        return data;
    }
    friend std::istream& operator>>(std::istream& is, Line& line)
    {
        return std::getline(is, line.data);
    }
};

void fcloseOrDie(FILE* file);

FILE* fopenOrDie(const char* filename, const char* instruction);

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file);

size_t count_lines(const std::string& fpath);

std::string timestamp();

void permute_plink(std::string& fin, const std::string& fout, uint gb, uint nbands);

void flip_UV(MyMatrix& U, MyMatrix& V, bool ubase = true);

void flip_Omg(MyMatrix& Omg2, MyMatrix& Omg);

void flip_Y(const MyMatrix& X, MyMatrix& Y);

double rmse(const MyMatrix& X, const MyMatrix& Y);

Eigen::VectorXd minSSE(const MyMatrix& X, const MyMatrix& Y);

double mev(const MyMatrix& X, const MyMatrix& Y);

void mev_rmse_byk(const MyMatrix& X, const MyMatrix& Y, MyVector& Vm, MyVector& Vr);

double get_median(std::vector<double> v);

#endif // PCAONE_UTILES_
