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
#include <clocale>

const double VAR_TOL = 1e-9;

using namespace std;
using namespace Eigen;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

template<typename M>
M load_csv (const std::string & path) {
    std::ifstream in;
    in.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(in, line)) {
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
    uint maxp = 20;  // maximum number of power iterations
    uint threads = 1;
    uint blocksize = 0;
    uint bands = 64;
    // for emu iteration
    uint maxiter = 100;
    double alpha = 0.001;
    // can be tol_emu or tol_pcangsd
    double tol = 1e-5;
    double tolmaf = 1e-4;
    // for arnoldi
    uint ncv = 20;   // max(20, 2*k + 1)
    uint imaxiter = 1000;
    double itol = 1e-6;
    // for halko
    uint oversamples = 10;
    double tol_halko = 1e-4;
    uint buffer = 2;

    double memory = 2; // 2 G
    bool printv = false;
    bool runem = false;
    bool batch = true; // if load all matrix into RAM.
    bool noshuffle = false;
    bool fast = false;
    bool emu = false;
    bool pcangsd  = false; // read GP field for PCAngsd instead of GT.
    bool halko = false;
    bool arnoldi = true;
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

void permute_plink2(string& fin, uint gb = 2);

void permute_plink(string& fin, uint blocksize=1);

void flip_UV(MatrixXd& U, MatrixXd& V, bool ubase = true);

void flip_Omg(MatrixXd& Omg2, MatrixXd& Omg);

void flip_Y(const MatrixXd& X, MatrixXd& Y);

double rmse(const MatrixXd& X, const MatrixXd& Y);

double mev(const MatrixXd& X, const MatrixXd& Y);

void mev_rmse_byk(const MatrixXd& X, const MatrixXd& Y, VectorXd& Vm, VectorXd& Vr);

#endif
