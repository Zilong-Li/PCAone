/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileCsv.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileCsv.hpp"

using namespace std;

// normalize count x of sample i according to --scale. libsize is only
// populated for the modes that need it (see csv_needs_libsize).
static inline double normalize_count(double x,
                                     int scale,
                                     const std::vector<double>& libsize,
                                     size_t i,
                                     double median_libsize,
                                     double scaleFactor) {
  if (!csv_needs_libsize(scale)) return x;
  const double total = libsize[i];
  if (total <= 0) return 0.0;  // empty sample: all its counts are zero, keep them at 0 instead of 0/0
  if (scale == 2) return log10(x * median_libsize / total + 1);
  if (scale == 3) return log1p(x / total * scaleFactor);
  return x / total * scaleFactor;  // scale == 4
}

// Offsets of the fields of one CSV line: field c is line[tidx[c], tidx[c+1]-1).
// Returns how many fields the line has. tidx holds ncol + 1 entries and nothing
// is written past them, so a line with too many fields is counted instead of
// overflowing tidx, as it did before the column count was checked.
static size_t csv_fields(const std::string& line, std::vector<size_t>& tidx, size_t ncol) {
  size_t nf = 1;
  tidx[0] = 0;
  for (size_t i = 0; i < line.size(); i++) {
    if (line[i] == ',') {
      if (nf < ncol) tidx[nf] = i + 1;
      nf++;
    }
  }
  tidx[ncol] = line.size() + 1;
  return nf;
}

// std::stof / std::stod on the field starting at line[b], without the throw: an
// exception inside an OpenMP loop cannot be caught and terminated the run
// ("what(): stof"). False if the field holds no finite number (NA, text, empty).
static inline bool csv_float(const std::string& line, size_t b, double& x) {
  const char* s = line.c_str() + b;
  char* end = nullptr;
  x = std::strtof(s, &end);
  return end != s && std::isfinite(x);
}
static inline bool csv_double(const std::string& line, size_t b, double& x) {
  const char* s = line.c_str() + b;
  char* end = nullptr;
  x = std::strtod(s, &end);
  return end != s && std::isfinite(x);
}

static void csv_check_ncol(size_t nf, size_t ncol, uint64 lineno) {
  if (nf != ncol)
    cao.error("line " + std::to_string(lineno) + " of the csv file has " + std::to_string(nf) + " columns, but " +
              std::to_string(ncol) + " were expected (from the first line, or --N)");
}

static void csv_bad_field(const std::string& line, const std::vector<size_t>& tidx, size_t col, uint64 lineno) {
  cao.error("line " + std::to_string(lineno) + ", column " + std::to_string(col + 1) +
            " of the csv file is not a number: '" + line.substr(tidx[col], tidx[col + 1] - tidx[col] - 1) + "'");
}

// parse one line of counts into column `col` of G
static void csv_parse_row(const std::string& line,
                          uint64 lineno,
                          std::vector<size_t>& tidx,
                          uint nsamples,
                          Mat2D& G,
                          Eigen::Index col,
                          int scale,
                          const std::vector<double>& libsize,
                          double median_libsize,
                          double scaleFactor) {
  csv_check_ncol(csv_fields(line, tidx, nsamples), nsamples, lineno);
  size_t bad = nsamples;
#pragma omp parallel for reduction(min : bad)
  for (size_t i = 0; i < nsamples; i++) {
    double entry;
    if (!csv_float(line, tidx[i], entry)) {
      if (i < bad) bad = i;
      continue;
    }
    G(i, col) = normalize_count(entry, scale, libsize, i, median_libsize, scaleFactor);
  }
  if (bad < nsamples) csv_bad_field(line, tidx, bad, lineno);
}

void FileCsv::read_all() {
  check_file_offset_first_var();

  auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
  auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
  size_t read, e, lastSNP = 0;
  zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
  std::string buffLine;
  G = Mat2D::Zero(nsamples, nsnps);
  while ((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin))) {
    ZSTD_inBuffer input = {buffIn, read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
      zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
      if (ZSTD_isError(zbuf.lastRet)) cao.error("Error: ZSTD decompression failed");
      buffCur += std::string((char*)buffOut, output.pos);
      while ((e = buffCur.find("\n")) != std::string::npos) {
        buffLine = buffCur.substr(0, e);
        buffCur.erase(0, e + 1);
        // G has nsnps columns, and --M can say fewer than the file has
        if (lastSNP >= nsnps) cao.error("the csv file has more than " + std::to_string(nsnps) + " lines (--M)");
        csv_parse_row(buffLine, lastSNP + 1, tidx, nsamples, G, lastSNP, params.scale, libsize, median_libsize,
                      params.scaleFactor);
        lastSNP++;
      }
    }
  }

  if (params.scale >= 1) standardize(G);  // standardization

  if (zbuf.lastRet != 0) cao.error("EOF before end of ZSTD_decompressStream.\n");

  // deal with the case there is no "\n" for the last line of file
  if (lastSNP != nsnps)
    cao.error("the csv file has " + std::to_string(lastSNP) + " lines, but " + std::to_string(nsnps) +
              " were expected (--M)");
}

void FileCsv::check_file_offset_first_var() {
  if (zbuf.fin == nullptr) {
    zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
  } else if (feof(zbuf.fin) || zbuf.lastRet == 0) {
    rewind(zbuf.fin);
  } else {
    rewind(zbuf.fin);
    if (params.verbose) cao.warn("confirm you are running the window-based RSVD (algorithm2)");
  }
  zbuf.lastRet = 1;
  buffCur = "";
}

void FileCsv::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {
  read_csvzstd_block(zbuf, buffCur, blocksize, start_idx, stop_idx, G, nsamples, libsize, tidx, median_libsize,
                     params.scale, params.scaleFactor);
}

void parse_csvzstd(ZstdDS& zbuf,
                   uint& nsamples,
                   uint& nsnps,
                   uint scale,
                   std::vector<double>& libsize,
                   std::vector<size_t>& tidx,
                   double& median_libsize) {
  auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
  auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
  size_t read, p, ncol = 0;
  int isEmpty = 1;
  const bool need_libsize = csv_needs_libsize(scale);
  nsnps = 0;
  std::string buffLine{""}, buffCur{""};
  while ((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin))) {
    isEmpty = 0;
    ZSTD_inBuffer input = {buffIn, read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
      zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
      if (ZSTD_isError(zbuf.lastRet)) cao.error("Error: ZSTD decompression failed");
      buffCur += std::string((char*)buffOut, output.pos);
      while ((p = buffCur.find("\n")) != std::string::npos) {
        nsnps++;
        buffLine = buffCur.substr(0, p);
        buffCur.erase(0, p + 1);
        // the first line sets the number of columns, and every line must match
        // it, including the second (the check used to start at the third line)
        if (nsnps == 1) {
          ncol = std::count(buffLine.begin(), buffLine.end(), ',') + 1;
          libsize.assign(ncol, 0.0);
          tidx.resize(ncol + 1);
        }
        csv_check_ncol(csv_fields(buffLine, tidx, ncol), ncol, nsnps);
        if (need_libsize)  // total counts per sample for --scale 2/3/4
        {
          size_t bad = ncol;
#pragma omp parallel for reduction(min : bad)
          for (size_t i = 0; i < ncol; i++) {
            double x;
            if (!csv_double(buffLine, tidx[i], x)) {
              if (i < bad) bad = i;
              continue;
            }
            libsize[i] += x;
          }
          if (bad < ncol) csv_bad_field(buffLine, tidx, bad, nsnps);
        }
      }
    }
  }

  if (isEmpty || nsnps == 0) cao.error("input file is empty.");
  if (zbuf.lastRet != 0) cao.error("EOF before end of ZSTD_decompressStream.");

  nsamples = ncol;
  zbuf.lastRet = 1;
  if (need_libsize) {
    std::vector<double> positive;
    positive.reserve(libsize.size());
    for (auto s : libsize) {
      if (s < 0) cao.error("negative total counts found. --scale 2, 3 and 4 expect non-negative counts.");
      if (s > 0) positive.push_back(s);
    }
    if (positive.empty()) cao.error("all samples have zero total counts. cannot normalize with --scale 2, 3 or 4.");
    if (positive.size() < libsize.size())
      cao.warn(std::to_string(libsize.size() - positive.size()) +
               " sample(s) have zero total counts. their normalized values are set to 0.");
    // median over non-empty samples, so empty samples do not drag the CPMED target to 0
    if (scale == 2) median_libsize = get_median(positive);
  }
}

void read_csvzstd_block(ZstdDS& zbuf,
                        std::string& buffCur,
                        uint blocksize,
                        uint64 start_idx,
                        uint64 stop_idx,
                        Mat2D& G,
                        uint nsamples,
                        std::vector<double>& libsize,
                        std::vector<size_t>& tidx,
                        double median_libsize,
                        uint scale,
                        double scaleFactor) {
  const uint actual_block_size = stop_idx - start_idx + 1;

  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = Mat2D::Zero(nsamples, actual_block_size);
  }
  auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
  auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
  size_t read, e, lastSNP = 0;
  std::string buffLine;
  if (buffCur != "") {
    while (lastSNP < actual_block_size && ((e = buffCur.find("\n")) != std::string::npos)) {
      buffLine = buffCur.substr(0, e);
      buffCur.erase(0, e + 1);
      csv_parse_row(buffLine, start_idx + lastSNP + 1, tidx, nsamples, G, lastSNP, scale, libsize, median_libsize,
                    scaleFactor);
      lastSNP++;
    }
  }

  if (zbuf.lastRet != 0 && lastSNP < actual_block_size) {
    while ((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin))) {
      ZSTD_inBuffer input = {buffIn, read, 0};
      while (input.pos < input.size) {
        ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
        zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
        if (ZSTD_isError(zbuf.lastRet)) cao.error("Error: ZSTD decompression failed");
        buffCur += std::string((char*)buffOut, output.pos);
        while (lastSNP < actual_block_size && ((e = buffCur.find("\n")) != std::string::npos)) {
          buffLine = buffCur.substr(0, e);
          buffCur.erase(0, e + 1);
          csv_parse_row(buffLine, start_idx + lastSNP + 1, tidx, nsamples, G, lastSNP, scale, libsize, median_libsize,
                        scaleFactor);
          lastSNP++;
        }
      }
      if (lastSNP >= actual_block_size) break;
    }
  }

  if (scale > 0) standardize(G);  // standardization
  if (lastSNP != actual_block_size) cao.error("something wrong when read_block_initial");
}

PermMat shuffle_csvzstd_to_bin(std::string& fin, std::string fout, uint gb, uint scale, double scaleFactor) {
  std::vector<size_t> tidx;
  std::vector<double> libsize;
  double median_libsize{0};
  uint nsnps, nsamples;
  const uint ibyte = 4;
  ZstdDS zbuf;
  {
    zbuf.fin = fopenOrDie(fin.c_str(), "rb");
    parse_csvzstd(zbuf, nsamples, nsnps, scale, libsize, tidx, median_libsize);
    fcloseOrDie(zbuf.fin);
  }
  uint64 bytes_per_snp = nsamples * ibyte;
  // in 64 bits: 1073741824 * gb wrapped to 0 at --buffer 4, then divided by it (SIGFPE)
  const uint64 bufsnps = (uint64)1073741824 * gb / bytes_per_snp;
  const uint blocksize = (uint)std::max<uint64>(1, std::min<uint64>(bufsnps, std::max<uint>(1, nsnps)));
  uint nblocks = (nsnps + blocksize - 1) / blocksize;
  std::ofstream ofs(fout + ".perm.bin", std::ios::binary);
  std::ofstream ofs2(fout + ".perm.txt");
  ofs.write((char*)&nsnps, ibyte);
  ofs.write((char*)&nsamples, ibyte);
  uint64 magic = ibyte * 2;
  zbuf.fin = fopenOrDie(fin.c_str(), "rb");
  zbuf.lastRet = 1;
  Mat2D G;
  std::vector<int> perm(nsnps);
  std::iota(perm.begin(), perm.end(), 0);
  auto rng = std::default_random_engine{};
  std::shuffle(perm.begin(), perm.end(), rng);
  Eigen::VectorXf fg;
  uint64 start_idx, stop_idx, idx;
  int ia{0}, ib{0};
  Eigen::VectorXi indices(nsnps);
  std::string buffCur{""};
  for (uint i = 0; i < nblocks; i++) {
    start_idx = i * blocksize;
    stop_idx = start_idx + blocksize - 1;
    stop_idx = stop_idx >= nsnps ? nsnps - 1 : stop_idx;
    read_csvzstd_block(zbuf, buffCur, blocksize, start_idx, stop_idx, G, nsamples, libsize, tidx, median_libsize, scale,
                       scaleFactor);
    for (Eigen::Index p = 0; p < G.cols(); p++, ia++) {
      ib = perm[ia];
      indices(ib) = ia;
      idx = magic + ib * bytes_per_snp;
      ofs.seekp(idx, std::ios_base::beg);
      fg = G.col(p).cast<float>();
      ofs.write((char*)fg.data(), bytes_per_snp);
    }
  }
  fin = fout + ".perm.bin";
  ofs2 << indices << "\n";
  // zstd_compress_file(fin, fin+".zst", 3);
  return PermMat(indices);
}
