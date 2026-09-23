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

void FileCsv::read_all() {
  check_file_offset_first_var();

  auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
  auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
  size_t read, i, j, e, lastSNP = 0;
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
        for (i = 0, j = 1; i < buffLine.size(); i++) {
          if (buffLine[i] == ',') {
            tidx[j] = i + 1;
            j++;
          }
        }
        tidx[nsamples] = buffLine.size() + 1;

#pragma omp parallel for
        for (size_t i = 0; i < nsamples; i++) {
          auto entry = std::stof(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
          G(i, lastSNP) = normalize_count(entry, params.scale, libsize, i, median_libsize, params.scaleFactor);
        }

        lastSNP++;
      }
    }
  }

  if (params.scale >= 1) standardize(G);  // standardization

  if (zbuf.lastRet != 0) cao.error("EOF before end of ZSTD_decompressStream.\n");

  // deal with the case there is no "\n" for the last line of file
  if (lastSNP != nsnps) cao.error("error when parsing csv file\n");
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
  size_t read, i, j, p, ncol = 0, lastCol = 0;
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
        lastCol = ncol;
        ncol = 1;
        for (i = 0, j = 1; i < buffLine.size(); i++) {
          if (buffLine[i] == ',') {
            ncol++;
            if (nsnps > 1) {
              tidx[j++] = i + 1;
            }
          }
        }
        // get ncol from the first line or header
        if (nsnps == 1) {
          libsize.assign(ncol, 0.0);
          tidx.resize(ncol + 1);
          for (i = 0, j = 1; i < buffLine.size(); i++) {
            if (buffLine[i] == ',') {
              tidx[j++] = i + 1;
            }
          }
        }

        tidx[ncol] = buffLine.size() + 1;
        if (need_libsize)  // total counts per sample for --scale 2/3/4
        {
#pragma omp parallel for
          for (size_t i = 0; i < ncol; i++) {
            libsize[i] += std::stod(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
          }
        }
        if (nsnps > 2 && (lastCol != ncol)) cao.error("the csv file has unaligned columns");
      }
    }
  }

  if (isEmpty) cao.error("input file is empty.");
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
  size_t read, i, j, e, lastSNP = 0;
  std::string buffLine;
  if (buffCur != "") {
    while (lastSNP < actual_block_size && ((e = buffCur.find("\n")) != std::string::npos)) {
      buffLine = buffCur.substr(0, e);
      buffCur.erase(0, e + 1);
      for (i = 0, j = 1; i < buffLine.size(); i++) {
        if (buffLine[i] == ',') {
          tidx[j++] = i + 1;
        }
      }
      tidx[nsamples] = buffLine.size() + 1;

#pragma omp parallel for
      for (size_t i = 0; i < nsamples; i++) {
        auto entry = std::stof(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
        G(i, lastSNP) = normalize_count(entry, scale, libsize, i, median_libsize, scaleFactor);
      }

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
          for (i = 0, j = 1; i < buffLine.size(); i++) {
            if (buffLine[i] == ',') {
              tidx[j] = i + 1;
              j++;
            }
          }
          tidx[nsamples] = buffLine.size() + 1;

#pragma omp parallel for
          for (size_t i = 0; i < nsamples; i++) {
            auto entry = std::stof(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
            G(i, lastSNP) = normalize_count(entry, scale, libsize, i, median_libsize, scaleFactor);
          }

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
  uint blocksize = 1073741824 * gb / bytes_per_snp;
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
