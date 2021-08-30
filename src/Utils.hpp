#ifndef __EMU_UTILES__
#define __EMU_UTILES__

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <cstring>
#include <string>
#include <chrono>
#include <random>
#include <cmath>
#include <algorithm>
#include <iterator>

using namespace std;
using namespace Eigen;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

template<typename M>
M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<float> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ' ')) {
            values.push_back(std::stof(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}

struct Param {
    string intype = ""; // bfile, bgen, beagle, csv
    string bed_prefix;
    string pgen_prefix;
    string bgen;
    string beagle;
    string csvfile;
    string outfile;
    uint k = 10;
    uint p = 20;  // maximum number of power iterations
    uint threads = 1;
    uint blocksize = 0;
    uint bands = 128;
    // for emu iteration
    uint maxiter = 100;
    double alpha = 0.001;
    double tol_emu = 1e-4;
    // for pcangsd
    double tol_pcangsd = 1e-4;
    double tolmaf = 1e-5;
    // can be tol_emu or tol_pcangsd
    double tol;
    // for arnoldi
    uint ncv = 20;   // max(20, 2*k + 1)
    uint imaxiter = 500;
    double itol = 1e-6;
    // for halko
    uint oversamples = 10;
    double tol_halko = 1e-4;

    double memory = 2; // 2 G
    bool runem = false;
    bool batch = true; // if load all matrix into RAM.
    bool fast = false;
    bool emu = false;
    bool pcangsd  = false; // read GP field for PCAngsd instead of GT.
    bool halko = true;
    bool arnoldi = false;
    bool verbose = false;
};

struct Line
{
    std::string data;
    operator std::string const&() const {return data;}
    friend std::istream& operator>>(std::istream& is, Line& line) {
        return std::getline(is, line.data);
    }
};

size_t count_lines(const string& fpath);
string timestamp();
void permute_plink(string& fin);
void flip_UV(MatrixXf& U, MatrixXf& V, bool ubase = true);
void flip_Y(const MatrixXf& X, MatrixXf& Y);
double rmse(const MatrixXf& X, const MatrixXf& Y);
double mev(const MatrixXf& X, const MatrixXf& Y);
void mev_rmse_byk(const MatrixXf& X, const MatrixXf& Y, VectorXd& Vm, VectorXd& Vr);

#endif
