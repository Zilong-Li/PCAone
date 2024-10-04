/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Utils.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Utils.hpp"

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
  return "Machine name: " + machine + "\nNode name: " + node +
         "\nOperating system release: " + release +
         "\nOperating system version: " + version +
         "\nOperating system name: " + sysname + "\n";
}

void fcloseOrDie(FILE* file) {
  if (!fclose(file)) {
    return;
  };
  /* error */
  perror("fclose error");
  exit(1);
}

FILE* fopenOrDie(const char* filename, const char* instruction) {
  FILE* const inFile = fopen(filename, instruction);
  if (inFile) return inFile;
  /* error */
  perror(filename);
  exit(1);
}

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file) {
  size_t const readSize = fread(buffer, 1, sizeToRead, file);
  if (readSize == sizeToRead) return readSize; /* good */
  if (feof(file)) return readSize;             /* good, reached end of file */
  /* error */
  perror("fread");
  exit(1);
}

size_t count_lines(const std::string& fpath) {
  std::ifstream in(fpath);
  if (!in.is_open()) throw invalid_argument("can not open " + fpath);
  size_t count = 0;
  std::string line;
  while (getline(in, line)) {
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
void flip_UV(MyMatrix& U, MyMatrix& V, bool ubase) {
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
          throw std::runtime_error(
              "the dimention of U and V have different k ranks.\n");
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
        throw std::runtime_error(
            "the dimention of U and V have different k ranks.\n");
      }
    }
  }
}

void flip_Omg(MyMatrix& Omg2, MyMatrix& Omg) {
  for (Eigen::Index i = 0; i < Omg.cols(); ++i) {
    // if signs of half of values are flipped then correct signs.
    if ((Omg2.col(i) - Omg.col(i)).array().abs().sum() >
        2 * (Omg2.col(i) + Omg.col(i)).array().abs().sum()) {
      Omg.col(i) *= -1;
    }
  }
  Omg2 = Omg;
}

void flip_Y(const MyMatrix& X, MyMatrix& Y) {
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    // if signs of half of values are flipped then correct signs.
    if ((X.col(i) - Y.col(i)).array().abs().sum() >
        2 * (X.col(i) + Y.col(i)).array().abs().sum()) {
      Y.col(i) *= -1;
    }
  }
}

double rmse(const MyMatrix& X, const MyMatrix& Y) {
  MyMatrix Z = Y;
  flip_Y(X, Z);
  return sqrt((X - Z).array().square().sum() / (X.cols() * X.rows()));
}

// Y is the truth matrix, X is the test matrix
Eigen::VectorXd minSSE(const MyMatrix& X, const MyMatrix& Y) {
  Eigen::Index w1, w2;
  Eigen::VectorXd res(X.cols());
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

double mev(const MyMatrix& X, const MyMatrix& Y) {
  double res = 0;
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    res += (X.transpose() * Y.col(i)).norm();
  }
  return res / X.cols();
}

void mev_rmse_byk(const MyMatrix& X, const MyMatrix& Y, MyVector& Vm,
                  MyVector& Vr) {
  for (Eigen::Index i = 0; i < X.cols(); ++i) {
    Vm(i) = 1 - mev(X.leftCols(i + 1), Y.leftCols(i + 1));
    Vr(i) = rmse(X.leftCols(i + 1), Y.leftCols(i + 1));
  }
}

double get_median(std::vector<double> v) {
  size_t n = v.size();
  if (n == 0) {
    return 0;
  } else {
    std::sort(v.begin(), v.end());
    if (n % 2 == 0) {
      return (v[n / 2 - 1] + v[n / 2]) / 2.0;
    } else {
      return v[n / 2];
    }
  }
}

std::vector<std::string> split_string(const std::string& s,
                                      const std::string& separators) {
  std::vector<std::string> ret;
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

void make_plink2_eigenvec_file(int K, std::string fout, const std::string& fin,
                               const std::string& fam) {
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

// parse eigvecs,eigvals and loadings
MyMatrix read_usv(const std::string& path) {
  const char sep = '\t';
  bool is_seperator[256] = {false};
  is_seperator[(unsigned int)sep] = true;
  double val;
  Double1D V;
  int j{0}, i{0}, k1{0}, k2{0}, begin{0};

  std::ifstream fin(path);
  std::string line;
  while (std::getline(fin, line)) {
    k1 = 0, k2 = 0;
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
    if (k1 > 0 && k2 > 0) {
      if (k1 != k2) cao.error("the columns are not aligned! =>" + path);
    }
    j++;
  }
  return Eigen::Map<MyMatrix>(V.data(), j, k1);
}
