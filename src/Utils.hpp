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
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

const double VAR_TOL = 1e-9;


typedef Eigen::MatrixXd MyMatrix;
typedef Eigen::VectorXd MyVector;
typedef Eigen::ArrayXXd MyArrayX;
typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

enum class FileType
{
    None,
    PLINK,
    CSV,
    BEAGLE,
    BGEN
};

struct Param
{
    FileType intype = FileType::None; // PLINK, CSV, BEAGLE, BGEN
    std::string bed_prefix = "";
    std::string pgen_prefix = "";
    std::string bgen = "";
    std::string beagle = "";
    std::string csvfile = "";
    std::string outfile = "pcaone";
    std::string tmpfile = "";
    //
    uint64 nsamples = 0;
    uint64 nsnps = 0;
    uint k = 10;
    uint maxp = 20; // maximum number of power iterations
    uint threads = 10;
    uint blocksize = 0;
    uint bands = 64;
    // for emu iteration
    uint maxiter = 100;
    double alpha = 0.001;
    // can be tol_emu or tol_pcangsd
    double tolem = 1e-4;
    double tolmaf = 1e-4;
    double maf = 0.0;
    // for arnoldi
    uint ncv = 20; // max(20, 2*k + 1)
    uint imaxiter = 1000;
    double itol = 1e-6;
    // for halko
    uint oversamples = 10;
    double tol = 1e-4;
    uint buffer = 2;

    double memory = 0; // 0 for disable
    bool cpmed = false;
    bool printv = false;
    bool runem = false;
    bool batch = true; // if load all matrix into RAM.
    bool noshuffle = false;
    bool fast = true;
    bool emu = false;
    bool pcangsd = false; // read GP field for PCAngsd instead of GT.
    bool halko = false;
    bool arnoldi = false;
    bool verbose = false;
    bool printu = false;
};

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

void permute_plink(std::string& fin, const std::string& fout, uint gb = 2);

void flip_UV(MyMatrix& U, MyMatrix& V, bool ubase = true);

void flip_Omg(MyMatrix& Omg2, MyMatrix& Omg);

void flip_Y(const MyMatrix& X, MyMatrix& Y);

double rmse(const MyMatrix& X, const MyMatrix& Y);

double mev(const MyMatrix& X, const MyMatrix& Y);

void mev_rmse_byk(const MyMatrix& X, const MyMatrix& Y, MyVector& Vm, MyVector& Vr);

double get_median(std::vector<double> v);

#endif // PCAONE_UTILES_
