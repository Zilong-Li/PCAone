/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileCsv.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileCsv.hpp"

#include "Eigen/src/Core/Matrix.h"
#include "zstdpp/zstdpp.hpp"

using namespace std;

void FileCsv::read_all() {
  check_file_offset_first_var();

  size_t read;
  size_t i, j, e, lastSNP = 0;

  std::string buffLine;
  G = Mat2D::Zero(nsamples, nsnps);

  while (!isEmpty) {
    read = res.readFrom(ifs);
    isEmpty = read == 0;

    ZSTD_inBuffer input = {res.getRawInData(), read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {res.getRawOutData(), res.getToWrite(), 0};
      lastRet = ctx(input, output);  // perform decompression
      if (ZSTD_isError(lastRet)) cao.error("Error: ZSTD decompression failed. lastRet:", lastRet);
      res.writeTo(buffCur, output.pos);

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
          G(i, lastSNP) = std::stod(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
          if (params.scale == 1)
            G(i, lastSNP) = log10(G(i, lastSNP) + 0.01);
          else if (params.scale == 2)
            G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
        }

        lastSNP++;
      }
    }
  }

  if (lastRet != 0) cao.error("EOF before end of stream:", lastRet);

  // deal with the case there is no "\n" for the last line of file
  if (lastSNP != nsnps) cao.error("when parsing csv file: lastSNP=", lastSNP, ", M=", nsnps);

  G.rowwise() -= G.colwise().mean();  // do centering
}

void FileCsv::check_file_offset_first_var() {
  if (!(isEmpty || ifs.eof() || !ifs.tellg())) {
    cao.warn("confirm you are running winSVD with advanced settings");
  }
  if (lastRet != 0) cao.error("EOF before end of stream:", lastRet);
  // rewind stream
  ifs.clear();
  ifs.seekg(0, std::ios::beg);
  // reset
  buffCur = "";
  isEmpty = 0;
}

void FileCsv::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {
  read_csvzstd_block(ifs, res, ctx, buffCur, isEmpty, lastRet, blocksize, start_idx, stop_idx, nsamples, tidx,
                     G);
  // if (scale == 1)
  //   G(i, lastSNP) = log10(G(i, lastSNP) + 0.01);
  // else if (scale == 2)
  //   G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
  G.rowwise() -= G.colwise().mean();  // do centering
}

void parse_csvzstd(const std::string &fin, uint32_t &nsamples, uint32_t &nsnps, uint scale,
                   std::vector<double> &libsize, std::vector<size_t> &tidx, double &median_libsize) {
  Resources res{};
  Context ctx{};
  // size_t const toRead = res.getToRead();
  size_t read, i, j, p, ncol = 0, lastCol = 0;
  int isEmpty{0};     // check if ifstream return nothing, which means EOF
  size_t lastRet{0};  // check the return value of decompress function
  std::ifstream ifs(fin, std::ios::binary);

  nsnps = 0;
  std::string buffLine{""}, buffCur{""};
  while (!isEmpty) {
    read = res.readFrom(ifs);
    isEmpty = read == 0;
    ZSTD_inBuffer input = {res.getRawInData(), read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {res.getRawOutData(), res.getToWrite(), 0};
      lastRet = ctx(input, output);  // perform decompression
      if (ZSTD_isError(lastRet)) cao.error("Error: ZSTD decompression failed");
      res.writeTo(buffCur, output.pos);

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
          libsize.resize(ncol);
          tidx.resize(ncol + 1);
          for (i = 0, j = 1; i < buffLine.size(); i++) {
            if (buffLine[i] == ',') {
              tidx[j++] = i + 1;
            }
          }
        }

        tidx[ncol] = buffLine.size() + 1;
        if (scale == 2)  // cpmed
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

  /* The last return value from ZSTD_decompressStream did not end on a
   * frame, but we reached the end of the file! We assume this is an
   * error, and the input was truncated.
   */
  if (lastRet != 0) cao.error("EOF before end of stream:", lastRet);

  nsamples = ncol;
  if (scale == 2) median_libsize = get_median(libsize);
}

template <typename MatrixType>
void read_csvzstd_block(std::ifstream &ifs, Resources &res, Context &ctx, std::string &buffCur, int &isEmpty,
                        size_t &lastRet, uint blocksize, uint64_t start_idx, uint64_t stop_idx,
                        uint32_t nsamples, std::vector<size_t> &tidx, MatrixType &G) {
  const uint actual_block_size = stop_idx - start_idx + 1;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = MatrixType::Zero(nsamples, actual_block_size);
  }

  size_t i, j, e, lastSNP = 0;
  std::string buffLine;

  if (buffCur != "") {  // we got buffers left to be processed
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
        G(i, lastSNP) = parseFromString<MatrixType>(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
      }
      lastSNP++;
    }
  }

  size_t read;
  if (lastSNP < actual_block_size) {
    while (!isEmpty) {
      read = res.readFrom(ifs);
      isEmpty = read == 0;
      ZSTD_inBuffer input = {res.getRawInData(), read, 0};

      while (input.pos < input.size) {
        ZSTD_outBuffer output = {res.getRawOutData(), res.getToWrite(), 0};
        lastRet = ctx(input, output);  // perform decompression
        if (ZSTD_isError(lastRet)) cao.error("Error: ZSTD decompression failed.");
        res.writeTo(buffCur, output.pos);

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
            G(i, lastSNP) = parseFromString<MatrixType>(buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
          }
          lastSNP++;
        }
      }
      if (lastSNP >= actual_block_size) break;
    }
  }

  if (lastSNP != actual_block_size) cao.error("something wrong when read_block_initial");
}

PermMat normCSV2BIN(std::string &fin, std::string fout, uint gb, uint scale, bool shuffle = true) {
  cao.print(tick.date(), "convert and normalize CSV to binary file");
  std::vector<size_t> tidx;
  std::vector<double> libsize;
  double median_libsize;
  uint32_t nsnps, nsamples;
  const size_t ibyte = sizeof(uint32_t);

  tick.clock();
  parse_csvzstd(fin, nsamples, nsnps, scale, libsize, tidx, median_libsize);
  cao.print(tick.date(), "elapsed time of normalizing data:", tick.reltime(), " seconds");

  std::vector<int> perm(nsnps);
  std::iota(perm.begin(), perm.end(), 0);
  std::string outZstdName = fout + ".perm.zst";
  std::ofstream ofs(outZstdName, std::ios::binary);
  size_t magic = 2 * sizeof(uint32_t);
  const std::uint8_t compressLevel = 3;

  if (shuffle) {
    auto rng = std::default_random_engine{};
    std::shuffle(perm.begin(), perm.end(), rng);
    ofs.write((char *)&nsnps, ibyte);
    ofs.write((char *)&nsamples, ibyte);
  } else {
    std::vector<uint8_t> frameData(magic);
    std::memcpy(frameData.data(), &nsnps, ibyte);
    std::memcpy(frameData.data() + ibyte, &nsamples, ibyte);
    auto compressed = zstdpp::compress(frameData, compressLevel);
    ofs.write(reinterpret_cast<const char *>(compressed.data()), compressed.size());
  }

  size_t bytes_per_snp = nsamples * ibyte;
  uint blocksize = 1073741824 * gb / bytes_per_snp;
  uint nblocks = (nsnps + blocksize - 1) / blocksize;
  std::ofstream ofs2(fout + ".perm.txt");

  Resources res{};
  Context ctx{};
  std::string buffCur{""};
  int isEmpty{0};     // check if ifstream return nothing, which means EOF
  size_t lastRet{0};  // check the return value of decompress function
  std::ifstream ifs(fin, std::ios::binary);

  Eigen::MatrixXi G;
  Eigen::VectorXf fg;
  uint64 start_idx, stop_idx, idx;
  int ia{0}, ib{0};
  Eigen::VectorXi indices(nsnps);
  std::vector<uint8_t> frameData(bytes_per_snp);

  for (uint i = 0; i < nblocks; i++) {
    start_idx = i * blocksize;
    stop_idx = start_idx + blocksize - 1;
    stop_idx = stop_idx >= nsnps ? nsnps - 1 : stop_idx;
    read_csvzstd_block(ifs, res, ctx, buffCur, isEmpty, lastRet, blocksize, start_idx, stop_idx, nsamples,
                       tidx, G);
    if (!shuffle) {  // compress stream
      for (Eigen::Index p = 0; p < G.cols(); p++, ia++) {
        indices(perm[ia]) = ia;
        // fg = G.col(p).cast<float>();
        std::memcpy(frameData.data(), G.col(p).data(), bytes_per_snp);
        auto compressed = zstdpp::compress(frameData, compressLevel);
        ofs.write(reinterpret_cast<const char *>(compressed.data()), compressed.size());
      }
    } else {
      for (Eigen::Index p = 0; p < G.cols(); p++, ia++) {
        ib = perm[ia];
        indices(ib) = ia;
        idx = magic + ib * bytes_per_snp;
        ofs.seekp(idx, std::ios_base::beg);
        // fg = G.col(p).cast<float>();
        // ofs.write((char *)fg.data(), bytes_per_snp);
        ofs.write((char *)G.col(p).data(), bytes_per_snp);
      }
    }
  }
  if (lastRet != 0) cao.error("EOF before end of stream:", lastRet);
  ofs2 << indices << "\n";

  if (shuffle) {  // shuffle to plain binary followed by zstd compression
    std::string tmpZstdName = fout + ".perm.zst.zst";
    zstdpp::stream_compress(outZstdName, tmpZstdName, compressLevel);
    moveFile(tmpZstdName, outZstdName);
  }

  fin = outZstdName;  // redirect input to new file
  cao.print(tick.date(), "done compressing data with zstd");

  return PermMat(indices);
}
