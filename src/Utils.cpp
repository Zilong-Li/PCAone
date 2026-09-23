/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Utils.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Utils.hpp"

#include <sys/utsname.h>

#include <cstddef>
#include <cstdlib>  // strtod
#include <cstring>  // strtok_r
#include <fstream>

#include "Common.hpp"
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

void standardize(Mat2D& X, double tol) {
  double sqrt_rdf = sqrt(X.rows() - 1.0);
  // if X is centered, then we can convert norm to sd
  for (Eigen::Index j = 0; j < X.cols(); j++) {
    double mean = X.col(j).mean();
    X.col(j).array() -= mean;
    double sd = X.col(j).norm() / sqrt_rdf;  // sd
    if (sd > tol) X.col(j) /= sd;
  }
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
  exit(4);  // error fread
}

size_t fwriteOrDie(const void* buffer, size_t sizeToWrite, FILE* file) {
  size_t const writtenSize = fwrite(buffer, 1, sizeToWrite, file);
  if (writtenSize == sizeToWrite) return sizeToWrite; /* good */
  /* error */
  perror("fwrite");
  exit(5);  // error fwrite
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

void make_plink2_eigenvec_from_psam(int K, const std::string& fout, const std::string& fin, const std::string& fpsam) {
  std::ifstream ipsam(fpsam);
  std::ifstream ifin(fin);
  std::ofstream ofs(fout);
  ofs << "#FID\tIID";
  for (int i = 0; i < K; i++) ofs << "\tPC" << i + 1;
  ofs << "\n";
  std::string header, line1, line2, sep{" \t"};
  // skip ## comment lines, stop at column header (starts with single '#')
  while (std::getline(ipsam, header))
    if (header.size() >= 1 && header[0] == '#' && (header.size() < 2 || header[1] != '#')) break;
  bool has_fid = (header.size() >= 4 && header.substr(0, 4) == "#FID");
  while (std::getline(ipsam, line1)) {
    if (line1.empty() || line1[0] == '#') continue;
    auto tokens = split_string(line1, sep);
    std::getline(ifin, line2);
    std::string fid = "0", iid;
    if (has_fid && tokens.size() >= 2) {
      fid = tokens[0];
      iid = tokens[1];
    } else if (!tokens.empty()) {
      iid = tokens[0];
    }
    ofs << fid << "\t" << iid << "\t" << line2 << "\n";
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
      const char* p = line.c_str();
      if (is_seperator[(uint8_t)line[i]] || i == (int)line.size()) {
        char* end;
        val = std::strtod(p + begin, &end);
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
  // a single-row file never enters the branch above, so k is still 0
  if (j == 1) k = (int)V.size();
  if (j < 1 || k < 1) cao.error("can not parse a matrix from file\n =>" + path);
  if ((int)V.size() != j * k) cao.error("the columns are not aligned!\n =>" + path);
  // V is filled row by row, so it must be mapped as row-major. Mapping it with
  // Eigen::Map<Mat2D> (column-major) transposes the contents of any file with
  // more than one column.
  using MatRowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  return Eigen::Map<MatRowMajor>(V.data(), j, k);
}

// parse .sigvals file
//
// header is  #nsamples,nsnps[,key=value]...
// The key=value fields were added later; std::stoi stops at the next comma, so
// an older PCAone reads a newer file correctly and simply ignores them, and a
// newer PCAone reading an older file just finds none.
void read_sigvals(const std::string& path, uint& N, uint& M, Mat1D& S, UsvTransform* transform) {
  double val;
  Double1D V;
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("can not open the singular values file\n => " + path);
  std::string line;
  getline(fin, line);
  // parse line #nsamples,nsnps
  size_t comma = line.find(',');
  N = std::stoi(line.substr(1, comma - 1));
  M = std::stoi(line.substr(comma + 1));
  if (transform != nullptr) {
    for (const auto& tok : split_string(line, ",")) {
      const size_t eq = tok.find('=');
      if (eq == std::string::npos) continue;
      const std::string key = tok.substr(0, eq), val_s = tok.substr(eq + 1);
      if (key == "scale") transform->scale = std::stoi(val_s);
      else if (key == "ploidy") transform->ploidy = std::stoi(val_s);
      else if (key == "gscale") transform->gscale = std::stoi(val_s);
      else continue;
      transform->known = true;
    }
  }
  // parse the rest
  while (getline(fin, line)) {
    val = std::stod(line);
    V.push_back(val);
  }
  S = Eigen::Map<Mat1D>(V.data(), V.size());
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
  cao.print(tick.date(), "read matrix from " + path);
  if (n <= 0 || k <= 0) cao.error("read_eigvecs: asked for a ", n, " x ", k, " matrix from ", path);

  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("can not open ", path);

  // Every bound is checked BEFORE the write it guards. The previous version
  // filled M(j, k1) for every line in the file and only compared j with n after
  // the loop, so a .eigvecs from a bigger cohort wrote past the end of the
  // matrix -- silently, since Eigen has no bounds check under -DNDEBUG. Its
  // inner loop had the mirror problem: it stopped only on k1 == k, so a short
  // line ran i past the end of the string, and the "columns not aligned" check
  // below it could never fire.
  Mat2D M(n, k);
  std::string line;
  int j = 0;
  while (std::getline(fin, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
    if (line.find_first_not_of(" \t") == std::string::npos) continue;  // blank line
    if (j >= n)
      cao.error(path, " has more than the ", n, " rows this dataset needs. is it from a different cohort?");
    const char* p = line.c_str();
    for (int c = 0; c < k; ++c) {
      char* end = nullptr;
      const double v = std::strtod(p, &end);  // strtod skips leading blanks itself
      if (end == p)
        cao.error(path, ": row ", j + 1, " has only ", c, " values but ", k,
                  " are needed. rerun the reference with a larger -k, or use a smaller one here");
      M(j, c) = v;
      p = end;
    }
    ++j;
  }
  if (j != n) cao.error(path, " has ", j, " rows but this dataset needs ", n);

  return M;
}

// parse AF
Mat1D read_frq(const std::string& path) {
  const std::string sep{" \t"};
  double val;
  Double1D V;
  std::ifstream fin(path);
  // without this the file simply reads as empty, F stays size 0, and the first
  // F(snp_idx) downstream segfaults with no message. -P/--USV points filebim at
  // <prefix>.mbim, which only a run with -D/--ld writes, so a prefix from a
  // plain PCA run lands here.
  if (!fin.is_open()) cao.error("can not open the allele frequency file\n => " + path);
  std::string line;
  while (getline(fin, line)) {
    auto tokens = split_string(line, sep);
    if ((int)tokens.size() != 7) cao.error("the input file is not valid!\n => " + path);
    val = std::stod(tokens[6]);
    V.push_back(val);
  }
  if (V.empty()) cao.error("no allele frequencies found in\n => " + path);
  return Eigen::Map<Mat1D>(V.data(), V.size());
}

std::string decode_beagle_allele(const std::string& allele) {
  if (allele == "0") return "A";
  if (allele == "1") return "C";
  if (allele == "2") return "G";
  if (allele == "3") return "T";
  return allele;
}

namespace {
String1D pvar_tokens_to_bim_tokens(const String1D& tokens, const std::string& path, const std::string& line) {
  if ((int)tokens.size() < 5) cao.error("the input pvar file is not valid!\n => " + path + "\n" + line);
  std::string alt = tokens[4];
  size_t comma = alt.find(',');
  if (comma != std::string::npos) alt = alt.substr(0, comma);
  // PGEN reader uses ALT allele dosage (allele_idx=1), so expose ALT as BIM A1.
  return String1D{tokens[0], tokens[2], "0", tokens[1], alt, tokens[3]};
}

std::string bim_match_key(const String1D& tokens, const std::string& path, const std::string& line) {
  if ((int)tokens.size() < 6) cao.error("the input bim file is not valid!\n => " + path + "\n" + line);
  return tokens[0] + "_" + tokens[3] + "_" + tokens[4] + "_" + tokens[5];
}
std::string bim_flip_key(const String1D& tokens, const std::string& path, const std::string& line) {
  if ((int)tokens.size() < 6) cao.error("the input bim file is not valid!\n => " + path + "\n" + line);
  return tokens[0] + "_" + tokens[3] + "_" + tokens[5] + "_" + tokens[4];
}
}  // namespace

std::string pvar_line_to_bim_line(const std::string& line, const std::string& path) {
  const std::string sep{" \t"};
  auto tokens = pvar_tokens_to_bim_tokens(split_string(line, sep), path, line);
  return tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3] + "\t" + tokens[4] + "\t" + tokens[5];
}

BimMatch match_bim_to_mbim(const std::string& bim_file, const std::string& mbim_file) {
  const std::string sep{" \t"};
  std::ifstream fb(bim_file), fm(mbim_file);
  if (!fb.is_open()) cao.error("can not open " + bim_file);
  if (!fm.is_open()) cao.error("can not open " + mbim_file);

  String1D bim_keys, mbim_keys;
  std::string line;
  while (getline(fb, line)) {
    auto tokens = split_string(line, sep);
    bim_keys.push_back(bim_match_key(tokens, bim_file, line));
  }
  std::unordered_map<std::string, int> mbim_lookup, mbim_flip_lookup;
  while (getline(fm, line)) {
    auto tokens = split_string(line, sep);
    if ((int)tokens.size() != 7) cao.error("the input file is not valid!\n => " + mbim_file);
    int idx = (int)mbim_keys.size();
    std::string key = bim_match_key(tokens, mbim_file, line);
    if (!mbim_lookup.insert({key, idx}).second)
      cao.error("duplicate SNP records found in " + mbim_file + " for key " + key);
    mbim_flip_lookup.insert({bim_flip_key(tokens, mbim_file, line), idx});
    mbim_keys.push_back(std::move(key));
  }

  BimMatch match;
  match.identical = bim_keys.size() == mbim_keys.size();
  if (match.identical) {
    for (int i = 0; i < (int)bim_keys.size(); ++i) {
      if (bim_keys[i] != mbim_keys[i]) {
        match.identical = false;
        break;
      }
    }
  }

  match.bim_indices.reserve(std::min(bim_keys.size(), mbim_keys.size()));
  match.mbim_indices.reserve(std::min(bim_keys.size(), mbim_keys.size()));
  match.flip.reserve(std::min(bim_keys.size(), mbim_keys.size()));
  for (int i = 0; i < (int)bim_keys.size(); ++i) {
    auto it = mbim_lookup.find(bim_keys[i]);
    if (it != mbim_lookup.end()) {
      match.bim_indices.push_back(i);
      match.mbim_indices.push_back(it->second);
      match.flip.push_back(false);
    } else {
      auto it2 = mbim_flip_lookup.find(bim_keys[i]);
      if (it2 != mbim_flip_lookup.end()) {
        match.bim_indices.push_back(i);
        match.mbim_indices.push_back(it2->second);
        match.flip.push_back(true);
      }
    }
  }

  return match;
}

BimMatch match_pvar_to_mbim(const std::string& pvar_file, const std::string& mbim_file) {
  const std::string sep{" \t"};
  std::ifstream fp(pvar_file), fm(mbim_file);
  if (!fp.is_open()) cao.error("can not open " + pvar_file);
  if (!fm.is_open()) cao.error("can not open " + mbim_file);

  String1D pvar_keys, mbim_keys;
  std::string line;
  while (getline(fp, line)) {
    if (line.empty() || line[0] == '#') continue;
    auto tokens = pvar_tokens_to_bim_tokens(split_string(line, sep), pvar_file, line);
    pvar_keys.push_back(bim_match_key(tokens, pvar_file, line));
  }
  std::unordered_map<std::string, int> mbim_lookup, mbim_flip_lookup;
  while (getline(fm, line)) {
    auto tokens = split_string(line, sep);
    if ((int)tokens.size() != 7) cao.error("the input file is not valid!\n => " + mbim_file);
    int idx = (int)mbim_keys.size();
    std::string key = bim_match_key(tokens, mbim_file, line);
    if (!mbim_lookup.insert({key, idx}).second)
      cao.error("duplicate SNP records found in " + mbim_file + " for key " + key);
    mbim_flip_lookup.insert({bim_flip_key(tokens, mbim_file, line), idx});
    mbim_keys.push_back(std::move(key));
  }

  BimMatch match;
  match.identical = pvar_keys.size() == mbim_keys.size();
  if (match.identical) {
    for (int i = 0; i < (int)pvar_keys.size(); ++i) {
      if (pvar_keys[i] != mbim_keys[i]) {
        match.identical = false;
        break;
      }
    }
  }

  match.bim_indices.reserve(std::min(pvar_keys.size(), mbim_keys.size()));
  match.mbim_indices.reserve(std::min(pvar_keys.size(), mbim_keys.size()));
  match.flip.reserve(std::min(pvar_keys.size(), mbim_keys.size()));
  for (int i = 0; i < (int)pvar_keys.size(); ++i) {
    auto it = mbim_lookup.find(pvar_keys[i]);
    if (it != mbim_lookup.end()) {
      match.bim_indices.push_back(i);
      match.mbim_indices.push_back(it->second);
      match.flip.push_back(false);
    } else {
      auto it2 = mbim_flip_lookup.find(pvar_keys[i]);
      if (it2 != mbim_flip_lookup.end()) {
        match.bim_indices.push_back(i);
        match.mbim_indices.push_back(it2->second);
        match.flip.push_back(true);
      }
    }
  }

  return match;
}

BimMatch match_beagle_to_mbim(const std::string& beagle_file, const std::string& mbim_file) {
  // Parse BEAGLE rows as marker_allele1_allele2 so matching respects allele order.
  gzFile fp = gzopen(beagle_file.c_str(), "r");
  if (!fp) cao.error("can not open " + beagle_file);
  uint64 bufsize = (uint64)128 * 1024 * 1024;
  char* original = (char*)calloc(bufsize, sizeof(char));
  char* buffer = original;
  const char* delims = "\t \n";
  tgets(fp, &buffer, &bufsize);  // skip header line
  if (buffer != original) original = buffer;
  buffer = original;
  String1D beagle_markers;
  while (tgets(fp, &buffer, &bufsize)) {
    if (buffer != original) original = buffer;
    char* marker = strtok_r(buffer, delims, &buffer);
    char* allele1 = strtok_r(NULL, delims, &buffer);
    char* allele2 = strtok_r(NULL, delims, &buffer);
    if (!marker || !allele1 || !allele2) cao.error("invalid BEAGLE record while matching markers:\n => " + beagle_file);
    beagle_markers.push_back(std::string(marker) + "_" + decode_beagle_allele(allele1) + "_" +
                             decode_beagle_allele(allele2));
    buffer = original;
  }
  gzclose(fp);
  free(original);

  // Parse mbim: chr id cM bp a1 a2 maf
  // Build two lookups: normal key chr_pos_a1_a2 and flipped key chr_pos_a2_a1
  const std::string sep{" \t"};
  std::ifstream fm(mbim_file);
  if (!fm.is_open()) cao.error("can not open " + mbim_file);
  std::unordered_map<std::string, int> mbim_lookup, mbim_flip_lookup;
  String1D mbim_keys;
  std::string line;
  while (getline(fm, line)) {
    auto tokens = split_string(line, sep);
    if ((int)tokens.size() != 7) cao.error("the input file is not valid!\n => " + mbim_file);
    int idx = (int)mbim_keys.size();
    std::string key = bim_match_key(tokens, mbim_file, line);
    if (!mbim_lookup.insert({key, idx}).second)
      cao.error("duplicate SNP records found in " + mbim_file + " for key " + key);
    mbim_flip_lookup.insert({bim_flip_key(tokens, mbim_file, line), idx});
    mbim_keys.push_back(std::move(key));
  }

  BimMatch match;
  match.identical = beagle_markers.size() == mbim_keys.size();
  if (match.identical) {
    for (int i = 0; i < (int)beagle_markers.size(); ++i) {
      if (beagle_markers[i] != mbim_keys[i]) {
        match.identical = false;
        break;
      }
    }
  }

  match.bim_indices.reserve(std::min(beagle_markers.size(), mbim_keys.size()));
  match.mbim_indices.reserve(std::min(beagle_markers.size(), mbim_keys.size()));
  match.flip.reserve(std::min(beagle_markers.size(), mbim_keys.size()));
  for (int i = 0; i < (int)beagle_markers.size(); ++i) {
    auto it = mbim_lookup.find(beagle_markers[i]);
    if (it != mbim_lookup.end()) {
      match.bim_indices.push_back(i);
      match.mbim_indices.push_back(it->second);
      match.flip.push_back(false);
    } else {
      auto it2 = mbim_flip_lookup.find(beagle_markers[i]);
      if (it2 != mbim_flip_lookup.end()) {
        match.bim_indices.push_back(i);
        match.mbim_indices.push_back(it2->second);
        match.flip.push_back(true);
      }
    }
  }

  return match;
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

// modified from https://github.com/facebook/zstd/blob/dev/examples/streaming_compression.c
void zstd_compress_file(const std::string& fname, std::string outname, int level = 3) {
  ZstdCS zbuf;  // zstd compression buffer
  zbuf.fout = fopenOrDie(outname.c_str(), "wb");

  // compression parameters
  ZSTD_CCtx_setParameter(zbuf.cctx, ZSTD_c_compressionLevel, level);
  ZSTD_CCtx_setParameter(zbuf.cctx, ZSTD_c_checksumFlag, 1);  // Add content checksum for integrity
  ZSTD_CCtx_setParameter(zbuf.cctx, ZSTD_c_nbWorkers, 0);     // single-threaded mode

  size_t const toRead = zbuf.buffInSize;
  auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
  auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));

  FILE* const fin = fopenOrDie(fname.c_str(), "rb");

  // this loop read one buffer chunk, compress it and write out
  for (;;) {
    size_t read = freadOrDie(buffIn, toRead, fin);
    /* Select the flush mode.
     * If the read may not be finished (read == toRead) we use
     * ZSTD_e_continue. If this is the last chunk, we use ZSTD_e_end.
     * Zstd optimizes the case where the first flush mode is ZSTD_e_end,
     * since it knows it is compressing the entire source in one pass.
     */
    int const lastChunk = (read < toRead);
    ZSTD_EndDirective const mode = lastChunk ? ZSTD_e_end : ZSTD_e_continue;
    /* Set the input buffer to what we just read.
     * We compress until the input buffer is empty, each time flushing the
     * output.
     */
    ZSTD_inBuffer input = {buffIn, read, 0};
    int finished;
    do {
      /* Compress into the output buffer and write all of the output to
       * the file so we can reuse the buffer next iteration.
       */
      ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
      zbuf.lastRet = ZSTD_compressStream2(zbuf.cctx, &output, &input, mode);
      if (ZSTD_isError(zbuf.lastRet)) cao.error("Error: ZSTD compression failed");
      fwriteOrDie(buffOut, output.pos, zbuf.fout);
      /* If we're on the last chunk we're finished when zstd returns 0,
       * which means its consumed all the input AND finished the frame.
       * Otherwise, we're finished when we've consumed all the input.
       */
      finished = lastChunk ? (zbuf.lastRet == 0) : (input.pos == input.size);
    } while (!finished);

    if (input.pos != input.size) cao.error("Impossible: zstd only returns 0 when the input is completely consumed!");
    if (lastChunk) break;
  }

  fcloseOrDie(fin);
}

void emMAF_with_GL(Mat1D& F, const Mat2D& P, int maxiter, double tolmaf) {
  uint nsnps = P.cols();
  uint nsamples = P.rows() / 2;
  Mat1D Ft = Mat1D::Zero(nsnps);
  double scale = 1.0 / (2.0 * nsamples);
  double diff;
  // run EM to estimate allele frequencies
  for (int it = 0; it < maxiter; it++) {
#pragma omp parallel for
    for (uint j = 0; j < nsnps; j++) {
      Ft(j) = F(j);
      double p0, p1, p2, pt = 0.0;
      for (uint i = 0; i < nsamples; i++) {
        p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
        p1 = P(2 * i + 1, j) * 2.0 * F(j) * (1.0 - F(j));
        p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
        pt += (p1 + 2.0 * p2) / (p0 + p1 + p2);
      }
      F(j) = pt * scale;
    }
    // calculate differences between iterations
    diff = sqrt((F - Ft).array().square().sum() / nsnps);
    // Check for convergence
    if (diff < tolmaf) {
      cao.print(tick.date(), "EM (MAF) converged at iteration:", it + 1);
      break;
    } else if (it == (maxiter - 1)) {
      cao.print(tick.date(), "EM (MAF) did not converge");
    }
  }
}

double qchisq(double p, int df) {
  if (df <= 0) return std::numeric_limits<double>::quiet_NaN();  // ADD THIS
  if (p <= 0.0) return 0.0;
  if (p >= 1.0) return std::numeric_limits<double>::infinity();

  double s = df / 2.0;

  // Initial guess using Wilson-Hilferty approximation
  // For chi-sq(df), approximate quantile:
  //   x ≈ df * (1 - 2/(9*df) + z_p * sqrt(2/(9*df)))^3
  // where z_p is the standard normal quantile of p.
  // Approximate z_p using rational approximation (Abramowitz & Stegun 26.2.23)
  double t;
  if (p < 0.5) {
    t = std::sqrt(-2.0 * std::log(p));
    t = t - (2.515517 + t * (0.802853 + t * 0.010328)) / (1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308)));
    t = -t;  // negative side
  } else {
    t = std::sqrt(-2.0 * std::log(1.0 - p));
    t = t - (2.515517 + t * (0.802853 + t * 0.010328)) / (1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308)));
  }

  double a = 2.0 / (9.0 * df);
  double x = df * std::pow(1.0 - a + t * std::sqrt(a), 3.0);
  if (x <= 0.0) x = 0.01;  // fallback

  // Log of chi-sq PDF normalization: log(2^(df/2) * Gamma(df/2))
  double log_norm = s * std::log(2.0) + std::lgamma(s);

  // Newton-Raphson iteration
  for (int iter = 0; iter < 100; ++iter) {
    double cdf = kf_gammap(s, x / 2.0);
    double err = cdf - p;

    // chi-sq PDF at x: f(x) = x^(s-1) * exp(-x/2) / (2^s * Gamma(s))
    double log_pdf = (s - 1.0) * std::log(x) - x / 2.0 - log_norm;
    double pdf = std::exp(log_pdf);

    // acceptable if converged, but should check |err| first
    if (pdf < 1e-300) break;  // avoid division by zero

    double delta = err / pdf;
    delta = std::max(delta, -x * 0.9);  // don't step to negative
    x -= delta;

    if (x <= 0.0) x = 1e-10;  // keep positive

    if (std::fabs(delta) < 1e-12 * (1.0 + x)) break;  // +1 for small x
  }

  return x;
}

double pchisq(double x, int df, bool lower_tail = true) {
  // Validate degrees of freedom
  if (df <= 0) return std::numeric_limits<double>::quiet_NaN();

  // Handle edge cases
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (x <= 0.0) return lower_tail ? 0.0 : 1.0;
  if (std::isinf(x)) return lower_tail ? 1.0 : 0.0;

  double s = df / 2.0;
  double z = x / 2.0;

  // Use the complementary function for numerical stability:
  // - When p is close to 1 (large x), kf_gammap loses precision
  // - When p is close to 0 (small x), kf_gammaq loses precision
  // Switch based on which tail is smaller to get best precision,
  // then flip if needed.
  double result;
  bool use_lower = (x < df + 1.0);  // heuristic: switch near mode

  if (use_lower) {
    result = kf_gammap(s, z);  // lower incomplete gamma
    return lower_tail ? result : 1.0 - result;
  } else {
    result = kf_gammaq(s, z);  // upper incomplete gamma
    return lower_tail ? 1.0 - result : result;
  }
}

void galinsky_selection_stat(Mat2D& V) {
// FastPCA/Galinsky statistic: M * v_{jk}^2 ~ chi-squared(df=1) under the null,
// where v_{jk} is the SNP loading for SNP j along PC k.
#pragma omp parallel for
  for (int j = 0; j < V.rows(); j++) {
    for (int i = 0; i < V.cols(); i++) {
      V(j, i) = V(j, i) * V(j, i) * V.rows();
    }
  }
}

namespace {

// ---------------------------------------------------------------------------
// Orthogonalized Gnanadesikan-Kettenring (Maronna & Zamar 2002), matching
// bigutilsr::dist_ogk() -- which is what pcadapt calls -- with its defaults
// niter = 2, beta = 0.9 and the tau scale of Yohai & Zamar (1998).
//
// The plain pairwise GK estimator this replaces was not even scale
// equivariant: rescaling one column of z-scores changed the Mahalanobis
// distances, so the statistic depended on an arbitrary normalization. The
// orthogonalization step fixes that, and the final reweighted estimate is a
// classical covariance of the retained rows, so it is positive definite by
// construction rather than by flooring eigenvalues.
// ---------------------------------------------------------------------------

constexpr double TAU_C1 = 4.5;
constexpr double TAU_C2 = 3.0;
constexpr double QNORM_3_4 = 0.6744897501960817;   // qnorm(3/4)
constexpr double INV_SQRT2 = 0.7071067811865475244;
constexpr double SQRT_2PI = 2.5066282746310005024;

// E[rho] for the consistency factor of the tau scale at the normal, exactly as
// robustbase::scaleTau2 computes it: Es2(c2) = Erho(c2 * qnorm(3/4)).
double tau_Erho(double b) {
  const double Phi = 0.5 * std::erfc(-b * INV_SQRT2);
  const double phi = std::exp(-0.5 * b * b) / SQRT_2PI;
  return 2.0 * ((1.0 - b * b) * Phi - b * phi + b * b) - 1.0;
}
const double TAU_ES2 = tau_Erho(TAU_C2 * QNORM_3_4);

// scaleTau2 with robustbase's defaults (c1 = 4.5, c2 = 3, consistency at the
// normal, iter = 1). Returns the scale, and the location too when asked.
// `buf` is scratch the caller owns; we reorder and overwrite it.
double scale_tau2(const double* x, Eigen::Index n, std::vector<double>& buf, double* mu_out = nullptr) {
  buf.assign(x, x + n);
  const double mu0 = median_inplace(buf);
  for (Eigen::Index i = 0; i < n; ++i) buf[i] = std::abs(x[i] - mu0);
  const double sigma0 = median_inplace(buf);  // MAD without the consistency factor
  if (!(sigma0 > 0.0)) {  // more than half the column is a single value
    if (mu_out) *mu_out = mu0;
    return 0.0;
  }

  // location: a one-step weighted mean, weights falling to zero at c1 * sigma0
  const double inv = 1.0 / (sigma0 * TAU_C1);
  double sw = 0.0, sxw = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    const double t = std::abs(x[i] - mu0) * inv;
    const double u = 1.0 - t * t;
    if (u > 0.0) {
      const double w = u * u;
      sw += w;
      sxw += x[i] * w;
    }
  }
  const double mu = (sw > 0.0) ? sxw / sw : mu0;

  // scale: the winsorized second moment about that location
  const double c2sq = TAU_C2 * TAU_C2;
  double rho = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    const double z = (x[i] - mu) / sigma0;
    rho += std::min(z * z, c2sq);
  }
  if (mu_out) *mu_out = mu;
  return sigma0 * std::sqrt(rho / ((double)n * TAU_ES2));
}

// One Maronna-Zamar step, in place: scale every column to unit tau scale, form
// the pairwise GK matrix of the scaled columns, and rotate onto its
// eigenvectors. `tmp` is scratch of the same shape as Z.
void ogk_step(Mat2D& Z, Mat2D& tmp) {
  const int p = (int)Z.cols();
  const Eigen::Index n = Z.rows();

#pragma omp parallel
  {
    std::vector<double> buf;
#pragma omp for schedule(static)
    for (int j = 0; j < p; ++j) {
      const double s = scale_tau2(Z.col(j).data(), n, buf);
      if (s > 0.0) Z.col(j) /= s;
    }
  }

  // the diagonal is 1 by construction, the columns now having unit tau scale
  Mat2D S = Mat2D::Identity(p, p);
  Int1D pi, pj;
  pi.reserve(p * (p - 1) / 2);
  pj.reserve(p * (p - 1) / 2);
  for (int i = 0; i < p; ++i)
    for (int j = i + 1; j < p; ++j) pi.push_back(i), pj.push_back(j);
  const int npairs = (int)pi.size();

#pragma omp parallel
  {
    std::vector<double> buf;
    Mat1D plus(n), minus(n);
#pragma omp for schedule(dynamic)
    for (int t = 0; t < npairs; ++t) {
      const int i = pi[t], j = pj[t];
      plus = Z.col(i) + Z.col(j);
      minus = Z.col(i) - Z.col(j);
      const double sp = scale_tau2(plus.data(), n, buf);
      const double sm = scale_tau2(minus.data(), n, buf);
      const double cij = 0.25 * (sp * sp - sm * sm);
      S(i, j) = cij;
      S(j, i) = cij;
    }
  }

  Eigen::SelfAdjointEigenSolver<Mat2D> eig(S);
  if (eig.info() != Eigen::Success) cao.error("failed eigendecomposition of the pcadapt OGK matrix.");
  tmp.noalias() = Z * eig.eigenvectors();
  Z.swap(tmp);
}

// Robust location and scatter by OGK, followed by the hard-rejection
// reweighting that bigutilsr::covrob_ogk() applies.
void robust_cov_ogk(const Mat2D& U, Mat1D& wcenter, Mat2D& wcov, int niter, double beta) {
  const int p = (int)U.cols();
  const Eigen::Index n = U.rows();

  Mat2D Z = U, tmp(n, p);
  for (int it = 0; it < niter; ++it) ogk_step(Z, tmp);
  tmp.resize(0, 0);

  // location and scale of the orthogonalized data; the distance is then just
  // the standardized sum of squares in that basis
  Mat1D ctr(p), sg(p);
#pragma omp parallel
  {
    std::vector<double> buf;
#pragma omp for schedule(static)
    for (int j = 0; j < p; ++j) {
      double m = 0.0;
      sg(j) = scale_tau2(Z.col(j).data(), n, buf, &m);
      ctr(j) = m;
    }
  }

  Mat1D d(n);
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j) {
      if (!(sg(j) > 0.0)) continue;  // degenerate direction carries no information
      const double t = (Z(i, j) - ctr(j)) / sg(j);
      acc += t * t;
    }
    d(i) = acc;
  }
  Z.resize(0, 0);

  std::vector<double> dv(d.data(), d.data() + n);
  const double cdelta = median_inplace(dv) / qchisq(0.5, p);
  const double cutoff = qchisq(beta, p) * cdelta;

  Eigen::Index nkeep = 0;
  for (Eigen::Index i = 0; i < n; ++i)
    if (d(i) < cutoff) ++nkeep;
  const bool use_all = (nkeep <= (Eigen::Index)p);
  if (use_all)
    cao.warn("pcadapt: the OGK reweighting kept only ", nkeep,
             " sites, too few to estimate a covariance; using all sites instead");

  auto kept = [&](Eigen::Index i) { return use_all || d(i) < cutoff; };
  const double nk = (double)(use_all ? n : nkeep);

  wcenter = Mat1D::Zero(p);
  for (Eigen::Index i = 0; i < n; ++i)
    if (kept(i)) wcenter += U.row(i).transpose();
  wcenter /= nk;

  wcov = Mat2D::Zero(p, p);
#pragma omp parallel
  {
    Mat2D acc = Mat2D::Zero(p, p);
    Mat1D v(p);
#pragma omp for schedule(static)
    for (Eigen::Index i = 0; i < n; ++i) {
      if (!kept(i)) continue;
      v = U.row(i).transpose() - wcenter;
      acc.selfadjointView<Eigen::Lower>().rankUpdate(v);
    }
#pragma omp critical
    wcov += acc;
  }
  for (int a = 0; a < p; ++a)  // rankUpdate only filled the lower triangle
    for (int b = a + 1; b < p; ++b) wcov(a, b) = wcov(b, a);
  wcov /= nk;
}
}  // namespace

void pcadapt_selection_stats(const Mat2D& Z, const std::vector<char>& keep, Mat1D& stat, Mat1D& chi2_stat,
                             Mat1D& pval, double& gif) {
  const int k = (int)Z.cols();
  const Eigen::Index n = Z.rows();

  // Sites with no residual variance have an undefined z-score. They must not
  // reach the robust covariance, the inflation median or the output: pcadapt
  // fits on `zscores[pass, ]` and leaves the rest NA, and so do we. The common
  // case drops nothing, and then no copy is made.
  std::vector<Eigen::Index> idx;
  for (Eigen::Index i = 0; i < n; ++i)
    if (keep[i]) idx.push_back(i);
  const Eigen::Index nk = (Eigen::Index)idx.size();
  if (nk <= k)
    cao.error("pcadapt: only ", nk, " of ", n,
              " sites have any residual variance, too few to fit a ", k, "-dimensional covariance");

  const bool subset = (nk != n);
  Mat2D Zs;
  if (subset) {
    Zs.resize(nk, k);
    for (Eigen::Index i = 0; i < nk; ++i) Zs.row(i) = Z.row(idx[i]);
  }
  const Mat2D& W = subset ? Zs : Z;

  Mat1D center;
  Mat2D cov;
  robust_cov_ogk(W, center, cov, 2, 0.9);  // bigutilsr::covrob_ogk defaults
  Mat2D inv_cov = cov.inverse();
  if (!inv_cov.allFinite())
    cao.error("pcadapt: the robust covariance of the z-scores is singular; try a smaller -k");

  Mat1D s(nk);
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < nk; ++i) {
    Mat1D dz = W.row(i).transpose() - center;
    s(i) = dz.transpose() * inv_cov * dz;
  }

  std::vector<double> stat_vec(s.data(), s.data() + nk);
  gif = median_inplace(stat_vec) / qchisq(0.5, k);  // stat_vec is ours to reorder
  if (!(gif > 0.0) || !std::isfinite(gif)) gif = 1.0;

  const double NA = std::numeric_limits<double>::quiet_NaN();
  stat = Mat1D::Constant(n, NA);
  chi2_stat = Mat1D::Constant(n, NA);
  pval = Mat1D::Constant(n, NA);
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < nk; ++i) {
    const Eigen::Index r = idx[i];
    stat(r) = s(i);
    chi2_stat(r) = s(i) / gif;
    pval(r) = pchisq(chi2_stat(r), k, false);
  }
}
