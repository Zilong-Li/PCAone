/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Utils.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Utils.hpp"

#include <sys/utsname.h>

#include <cstddef>
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

// parse .sigvals file
void read_sigvals(const std::string& path, uint& N, uint& M, Mat1D& S) {
  double val;
  Double1D V;
  std::ifstream fin(path);
  std::string line;
  getline(fin, line);
  // parse line #nsamples,nsnps
  size_t comma = line.find(',');
  N = std::stoi(line.substr(1, comma - 1));
  M = std::stoi(line.substr(comma + 1));
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

    if (input.pos != input.size)
      cao.error("Impossible: zstd only returns 0 when the input is completely consumed!");
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
