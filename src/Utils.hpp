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
    string intype = ""; // beagle, bfile, pfile, bgen
    string bed_prefix;
    string pgen_prefix;
    string bgen;
    string beagle;
    string outfile;
    string truefile;
    uint k;
    uint p = 2;
    uint threads = 1;
    uint blocksize = 0;
    // for emu iteration
    uint maxiter = 100;
    double alpha = 0.001;
    double tol = 5e-7;
    // for pcangsd
    double tol2 = 1e-4;
    double tolmaf = 1e-5;
    // for arnoldi
    uint ncv = 20;   // max(20, 2*k + 1)
    uint imaxiter = 500;
    double itol = 1e-6;
    // for halko
    uint oversamples = 10;

    bool batch = true; // if load all matrix into RAM.
    bool emu = false;
    bool pcangsd  = false; // read GP field for PCAngsd instead of GT.
    bool halko = true;
    bool arnoldi = false;
    bool verbose = false;
};

class MeasureTime {

  public:
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
    time_t start_time_info, end_time_info;

    void init() {
      auto start = std::chrono::system_clock::now(); // wall clock
      start_time_info = std::chrono::system_clock::to_time_t( start );
      begin = std::chrono::steady_clock::now(); // to measure elapsed time
    }

    void stop(){
      auto endtime = std::chrono::system_clock::now();
      end_time_info = std::chrono::system_clock::to_time_t( endtime );
      end = std::chrono::steady_clock::now();
    }

    MeasureTime(void);
    ~MeasureTime(void);
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

#endif
