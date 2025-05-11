/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Utils.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Utils.hpp"

#include <sys/utsname.h>  // POSIX

#include <cstring>  // strtok_r
#include <fstream>

#include "kfunc.h"

using namespace std;

std::string get_machine() {
  struct utsname unameData;
  if (uname(&unameData) != 0) {
    perror("uname");
    exit(EXIT_FAILURE);
  }
  std::string machine{unameData.machine};
  std::string node{unameData.nodename};
  std::string release{unameData.release};
  std::string version{unameData.version};
  std::string sysname{unameData.sysname};
  return "Machine name: " + machine + "\nNode name: " + node + "\nOperating system release: " + release +
         "\nOperating system version: " + version + "\nOperating system name: " + sysname + "\n";
}

int tgets(gzFile gz, char** buf, uint64* size) {
  int rlen = 0;
  char* tok = gzgets(gz, *buf + rlen, *size - rlen);  // return buf or NULL
  if (!tok) return rlen;
  int tmp = tok ? strlen(tok) : 0;
  if (tok[tmp - 1] != '\n') {
    // expand buf size if no end-of-line found
    rlen += tmp;
    *size *= 2;
    *buf = (char*)realloc(*buf, *size);
  }
  rlen += tmp;
  return rlen;
}

size_t count_lines(const std::string& fpath) {
  std::ifstream fin(fpath);
  if (!fin.is_open()) throw invalid_argument("can not open " + fpath);
  size_t count = 0;
  std::string line;
  while (getline(fin, line)) {
    count++;
  }
  return count;
}

std::string timestamp() {
  auto t1 = std::chrono::system_clock::now();
  std::time_t tc = std::chrono::system_clock::to_time_t(t1);
  std::string str(std::ctime(&tc));
  str.pop_back();  // str[str.size() - 1] = '.';
  str = std::string("[") + str + std::string("] ");
  return str;
}

// Sign correction to ensure deterministic output from SVD.
// see https://www.kite.com/python/docs/sklearn.utils.extmath.svd_flip
void flip_UV(Mat2D& U, Mat2D& V, bool ubase) {
  if (ubase) {
    Eigen::Index x, i;
    for (i = 0; i < U.cols(); ++i) {
      U.col(i).cwiseAbs().maxCoeff(&x);
      if (U(x, i) < 0) {
        U.col(i) *= -1;
        if (V.cols() == U.cols()) {
          V.col(i) *= -1;
        } else if (V.rows() == U.cols()) {
          V.row(i) *= -1;
        } else {
          cao.error("the dimention of U and V have different k ranks.\n");
        }
      }
    }
  } else {
    Eigen::Index x, i;
    for (i = 0; i < V.cols(); ++i) {
      if (V.cols() == U.cols()) {
        V.col(i).cwiseAbs().maxCoeff(&x);
        if (V(x, i) < 0) {
          U.col(i) *= -1;
          V.col(i) *= -1;
        }
      } else if (V.rows() == U.cols()) {
        V.row(i).cwiseAbs().maxCoeff(&x);
        if (V(i, x) < 0) {
          U.col(i) *= -1;
          V.row(i) *= -1;
        }
      } else {
        cao.error("the dimention of U and V have different k ranks.\n");
      }
    }
  }
}

void flip_Y(const Mat2D& X, Mat2D& Y) {
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    // if signs of half of values are flipped then correct signs.
    if ((X.col(i) - Y.col(i)).array().abs().sum() > 2 * (X.col(i) + Y.col(i)).array().abs().sum()) {
      Y.col(i) *= -1;
    }
  }
}

double rmse(const Mat2D& X, const Mat2D& Y) {
  Mat2D Z = Y;
  flip_Y(X, Z);
  return sqrt((X - Z).array().square().sum() / (X.cols() * X.rows()));
}

double rmse1d(const Mat1D& x, const Mat1D& y) { return sqrt((x - y).array().square().sum() / x.size()); }

// Y is the truth matrix, X is the test matrix
Mat1D minSSE(const Mat2D& X, const Mat2D& Y) {
  Eigen::Index w1, w2;
  Mat1D res(X.cols());
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    // test against the original matrix to find the index with mincoeff
    ((-Y).colwise() + X.col(i)).array().square().colwise().sum().minCoeff(&w1);
    // test against the flipped matrix with the opposite sign
    (Y.colwise() + X.col(i)).array().square().colwise().sum().minCoeff(&w2);
    // get the minSSE value for X.col(i) against -Y.col(w1)
    auto val1 = (-Y.col(w1) + X.col(i)).array().square().sum();
    // get the minSSE value for X.col(i) against Y.col(w2)
    auto val2 = (Y.col(w2) + X.col(i)).array().square().sum();
    if (w1 != w2 && val1 > val2)
      res[i] = val2;
    else
      res[i] = val1;
  }
  return res;
}

double mev(const Mat2D& X, const Mat2D& Y) {
  double res = 0;
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    res += (X.transpose() * Y.col(i)).norm();
  }
  return res / X.cols();
}

void mev_rmse_byk(const Mat2D& X, const Mat2D& Y, Mat1D& Vm, Mat1D& Vr) {
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    Vm(i) = 1 - mev(X.leftCols(i + 1), Y.leftCols(i + 1));
    Vr(i) = rmse(X.leftCols(i + 1), Y.leftCols(i + 1));
  }
}

double get_median(std::vector<double> v) {
  size_t n = v.size();
  if (n == 0) {
    return 0;
  }
  std::sort(v.begin(), v.end());
  if (n % 2 == 0) {
    return (v[n / 2 - 1] + v[n / 2]) / 2.0;
  } else {
    return v[n / 2];
  }
}

String1D split_string(const std::string& s, const std::string& separators) {
  String1D ret;
  bool is_seperator[256] = {false};
  for (auto& ch : separators) {
    is_seperator[(unsigned int)ch] = true;
  }
  int begin = 0;
  for (int i = 0; i <= (int)s.size(); i++) {
    if (is_seperator[(uint8_t)s[i]] || i == (int)s.size()) {
      ret.push_back(std::string(s.begin() + begin, s.begin() + i));
      begin = i + 1;
    }
  }
  return ret;
}

void make_plink2_eigenvec_file(int K, std::string fout, const std::string& fin, const std::string& fam) {
  std::ifstream ifam(fam);
  std::ifstream ifin(fin);
  std::ofstream ofs(fout);
  ofs << "#FID\tIID";
  for (int i = 0; i < K; i++) ofs << "\tPC" << i + 1;
  ofs << "\n";
  std::string line1, line2, sep{" \t"};
  while (getline(ifam, line1)) {
    auto tokens = split_string(line1, sep);
    getline(ifin, line2);
    ofs << tokens[0] + "\t" + tokens[1] + "\t" << line2 << std::endl;
  }
}

bool isZstdCompressed(const char* filename) {
  FILE* file = fopen(filename, "rb");
  if (!file) return false;

  char magicNumber[4];
  if (fread(magicNumber, 1, 4, file) != 4) {
    fclose(file);
    return false;
  }

  bool isCompressed = (ZSTD_isFrame(magicNumber, 4) != 0);

  fclose(file);
  return isCompressed;
}

// parse table file generally
Mat2D read_usv(const std::string& path) {
  const char sep = '\t';
  bool is_seperator[256] = {false};
  is_seperator[(unsigned int)sep] = true;
  double val;
  Double1D V;
  int j{0}, i{0}, k1{0}, k2{0}, begin{0}, k{0};

  std::ifstream fin(path);
  std::string line;
  while (std::getline(fin, line)) {
    for (begin = 0, i = 0; i <= (int)line.size(); i++) {
      if (is_seperator[(uint8_t)line[i]] || i == (int)line.size()) {
        val = std::stod(std::string(line.begin() + begin, line.begin() + i));
        V.push_back(val);
        begin = i + 1;
        if (j % 2 == 1) {
          k1++;
        } else {
          k2++;
        }
      }
    }
    if (j > 0 && k1 > 0 && k2 > 0) {
      if (k1 != k2) cao.error("the columns are not aligned!\n =>" + path);
      k = k1;
      k1 = 0, k2 = 0;
    }
    j++;
  }
  return Eigen::Map<Mat2D>(V.data(), j, k);
}

// parse .eigvals file
Mat1D read_eigvals(const std::string& path) {
  double val;
  Double1D V;
  std::ifstream fin(path);
  std::string line;
  while (getline(fin, line)) {
    val = std::stod(line);
    V.push_back(val);
  }
  return Eigen::Map<Mat1D>(V.data(), V.size());
}

// parse .eigvecs or .loadings file assume rows and cols are known
Mat2D read_eigvecs(const std::string& path, int n, int k) {
  Mat2D M(n, k);
  const char sep = '\t';
  bool is_seperator[256] = {false};
  is_seperator[(unsigned int)sep] = true;
  int begin{0}, j{0}, i{0}, k1{0};
  std::ifstream fin(path);
  std::string line;
  while (std::getline(fin, line)) {
    for (begin = 0, k1 = 0, i = 0; i <= (int)line.size(); i++) {
      if (is_seperator[(uint8_t)line[i]] || i == (int)line.size()) {
        M(j, k1) = std::stod(std::string(line.begin() + begin, line.begin() + i));
        begin = i + 1;
        k1++;
      }
    }
    if (k1 != k) cao.error("the columns are not aligned!\n =>" + path);
    j++;
  }

  if (j != n) cao.error("the number of rows differs from the input!\n =>" + path);

  return M;
}

// parse AF
Mat1D read_frq(const std::string& path) {
  const std::string sep{"\t"};
  double val;
  Double1D V;
  std::ifstream fin(path);
  std::string line;
  while (getline(fin, line)) {
    auto tokens = split_string(line, sep);
    if ((int)tokens.size() != 7) cao.error("the input file is not valid!\n => " + path);
    val = std::stod(tokens[6]);
    V.push_back(val);
  }
  return Eigen::Map<Mat1D>(V.data(), V.size());
}

void check_bim_vs_mbim(const std::string& bim_file, const std::string& mbim_file) {
  const std::string sep{"\t"};
  std::ifstream fb(bim_file), fm(mbim_file);
  std::string lb, lm;
  while (getline(fb, lb) && getline(fm, lm)) {
    auto id1 = split_string(lb, sep)[1];
    auto id2 = split_string(lm, sep)[1];
    if (id1 != id2)
      cao.error(
          "the bim files are not matched!\n"
          "id1 vs id2 => " +
          id1 + " vs " + id2);
  }
}

void parse_beagle_file(Mat2D& P, gzFile fp, const int nsamples, const int nsnps) {
  // Mat2D P(nsamples * 2, nsnps);  // genotype likelihood
  const char* delims = "\t \n";
  char *original, *buffer, *tok;
  uint64 bufsize = (uint64)128 * 1024 * 1024;
  original = buffer = (char*)calloc(bufsize, sizeof(char));
  tgets(fp, &buffer, &bufsize);
  int i = 0, j = 0;
  // read all GL data into P
  while (tgets(fp, &buffer, &bufsize)) {
    if (buffer != original) original = buffer;
    tok = strtok_r(buffer, delims, &buffer);
    tok = strtok_r(NULL, delims, &buffer);
    tok = strtok_r(NULL, delims, &buffer);
    for (i = 0; i < nsamples; i++) {
      tok = strtok_r(NULL, delims, &buffer);
      P(2 * i + 0, j) = strtod(tok, NULL);
      tok = strtok_r(NULL, delims, &buffer);
      P(2 * i + 1, j) = strtod(tok, NULL);
      tok = strtok_r(NULL, delims, &buffer);
    }
    buffer = original;
    j++;
  }
  free(buffer);
  if (nsnps != j) {
    cao.error("something wrong parsing beagle");
  }
}

String1D parse_beagle_samples(const std::string& fin) {
  const char* delims = "\t \n";
  char *buffer, *tok;
  uint64 bufsize = (uint64)128 * 1024 * 1024;
  buffer = (char*)calloc(bufsize, sizeof(char));
  gzFile fp = gzopen(fin.c_str(), "r");
  tgets(fp, &buffer, &bufsize);
  strtok_r(buffer, delims, &buffer);
  int nCol = 1;
  String1D res;
  while ((tok = strtok_r(NULL, delims, &buffer))) {
    nCol++;
    if ((nCol - 1) % 3 == 0) res.push_back(std::string(tok));
  }
  gzclose(fp);
  if (nCol % 3) cao.error("Number of columns should be a multiple of 3.");
  return res;
}

void write_eigvecs2_beagle(const Mat2D& U, const std::string& fin, const std::string& fout) {
  std::ofstream feig2(fout);
  if (!feig2.is_open()) cao.error("can not open " + fout);
  feig2 << "#FID\tIID";
  int i = 0;
  for (i = 0; i < U.rows(); i++) feig2 << "\tPC" << i + 1;
  feig2 << "\n";
  i = 0;
  for (const auto& s : parse_beagle_samples(fin)) {
    feig2 << s << "\t" << s << "\t" << U.row(i) << "\n";
    i++;
  }
}

double chisq1d(const double x) {
  double p = kf_gammaq(1.0 / 2.0, x / 2.0);  // nan expected
  return std::isnan(p) ? 1.0 : p;            // if nan, then retrun 1.0
}

void moveFile(const std::filesystem::path& source, const std::filesystem::path& destination) {
  try {
    if (std::filesystem::exists(destination)) {
      std::filesystem::remove(destination);  // Remove the destination file if it exists
    }
    std::filesystem::rename(source, destination);  // Move the file to the new location
  } catch (const std::filesystem::filesystem_error& e) {
    cao.error("Filesystem error: ", e.what());
  } catch (const std::exception& e) {
    cao.error("Error: ", e.what());
  }
}
