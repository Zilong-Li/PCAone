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

bool isZstdCompressed(const char* filename);

MyMatrix read_usv(const std::string& path);

#endif  // PCAONE_UTILES_
