#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include <sys/utsname.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

#include "Common.hpp"
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

std::string get_machine();

void fcloseOrDie(FILE* file);

FILE* fopenOrDie(const char* filename, const char* instruction);

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file);

size_t count_lines(const std::string& fpath);

std::string timestamp();

void flip_UV(Mat2D& U, Mat2D& V, bool ubase = true);

void flip_Omg(Mat2D& Omg2, Mat2D& Omg);

void flip_Y(const Mat2D& X, Mat2D& Y);

double rmse(const Mat2D& X, const Mat2D& Y);

Eigen::VectorXd minSSE(const Mat2D& X, const Mat2D& Y);

double mev(const Mat2D& X, const Mat2D& Y);

void mev_rmse_byk(const Mat2D& X, const Mat2D& Y, Mat1D& Vm, Mat1D& Vr);

double get_median(std::vector<double> v);

std::vector<std::string> split_string(const std::string& s,
                                      const std::string& separators);

void make_plink2_eigenvec_file(int K, std::string fout, const std::string& fin,
                               const std::string& fam);

bool isZstdCompressed(const char* filename);

Mat2D read_usv(const std::string& path);

Mat1D read_eigvals(const std::string& path);

Mat2D read_eigvecs(const std::string& path, int n, int k);

Mat1D read_frq(const std::string& path);

void check_bim_vs_mbim(const std::string& bim_file,
                       const std::string& mbim_file);
#endif  // PCAONE_UTILES_
