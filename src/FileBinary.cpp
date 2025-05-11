#include "FileBinary.hpp"

using namespace std;

void FileBin::check_file_offset_first_var() {
  if (is_zstd) {  // rewind stream
    if (lastRet != 0) cao.error("EOF before end of stream:", lastRet);
    ifs.clear();
    ifs.seekg(0, std::ios::beg);
    return;
  }
  // magic += missing_points.size() * sizeof(uint64);
  long long offset = ibyte * 2 + nsnps * bytes_per_snp;
  if (ifs.tellg() == offset) {
    // reach the end of bed, reset the position to the first variant;
    ifs.seekg(ibyte * 2, std::ios_base::beg);
  } else if (ifs.tellg() == ibyte * 2) {
    ;
  } else {
    ifs.seekg(ibyte * 2, std::ios_base::beg);
    cao.warn("confirm you are running winSVD with advanced settings");
  }
}

void FileBin::read_all() {
  check_file_offset_first_var();
  if (!is_zstd) {
    G = Mat2D::Zero(nsamples, nsnps);
    Eigen::VectorXf fg(nsamples);
    for (Eigen::Index i = 0; i < G.cols(); i++) {
      ifs.read((char *)fg.data(), bytes_per_snp);
      G.col(i) = fg.cast<double>();
      G.col(i).array() -= G.col(i).mean();
    }
    return;
  }

  G = Mat2D::Zero(nsamples, nsnps);
  std::vector<float> snpLine(nsamples);  // store samples of current SNP
  size_t lastSNP = 0;
  bool no_magic = false;
  size_t read;
  while (!isEmpty) {
    read = res.readFrom(ifs);
    isEmpty = read == 0;
    ZSTD_inBuffer input = {res.getRawInData(), read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {res.getRawOutData(), res.getToWrite(), 0};
      lastRet = ctx(input, output);  // perform decompression
      if (ZSTD_isError(lastRet)) cao.error("Error: ZSTD decompression failed. lastRet:", lastRet);
      std::vector<char> _tmp((char *)res.getRawOutData(), (char *)res.getRawOutData() + output.pos);
      decompressedBytes.insert(decompressedBytes.end(), _tmp.begin(), _tmp.end());
      if (no_magic) {
        while (decompressedBytes.size() >= bytes_per_snp) {  // parse values
          decompressedBytes.erase(decompressedBytes.begin(), decompressedBytes.begin() + bytes_per_snp);
          std::memcpy(snpLine.data(), decompressedBytes.data(), bytes_per_snp);
#pragma omp parallel for
          for (uint32_t i = 0; i < nsamples; i++) {
            G(i, lastSNP) = snpLine[i];  // TODO:  maybe do scaling
          }
          lastSNP++;
        }
      } else {
        while (decompressedBytes.size() >= bytes_per_snp + 8) {  // parse values
          decompressedBytes.erase(decompressedBytes.begin() + 8,
                                  decompressedBytes.begin() + bytes_per_snp + 8);
          std::memcpy(snpLine.data(), decompressedBytes.data(), bytes_per_snp);
#pragma omp parallel for
          for (uint32_t i = 0; i < nsamples; i++) {
            G(i, lastSNP) = snpLine[i];  // TODO:  maybe do scaling
          }
          lastSNP++;
        }
        // remove the magic now
        decompressedBytes.erase(decompressedBytes.begin(), decompressedBytes.begin() + 8);
        no_magic = true;
      }
    }
  }

  if (lastRet != 0) cao.error("EOF before end of stream:", lastRet);
  if (lastSNP != nsnps) cao.error("when parsing binary zstd file: lastSNP=", lastSNP, ", M=", nsnps);
  G.rowwise() -= G.colwise().mean();  // do centering

  return;
}

// TODO : can standardize
void FileBin::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  G = Mat2D(nsamples, actual_block_size);
  if (!is_zstd) {
    // check where we are
    long long offset = ibyte * 2 + start_idx * bytes_per_snp;
    if (ifs.tellg() != offset) cao.error("something wrong with read_snp_block!\n");
    Eigen::VectorXf fg(nsamples);
    for (Eigen::Index i = 0; i < G.cols(); i++) {
      ifs.read((char *)fg.data(), bytes_per_snp);
      G.col(i) = fg.cast<double>();
      G.col(i).array() -= G.col(i).mean();
    }
    return;
  }

  size_t lastSNP = 0;
  std::vector<float> snpLine(nsamples);  // store samples of current SNP
  bool no_magic = false;
  if (start_idx != 0) no_magic = true;  // no magic if it's not the first block

  if (decompressedBytes.size()) {  // process the left bytes from previous reading
    while (lastSNP < actual_block_size && decompressedBytes.size() >= bytes_per_snp) {  // assure no magic
      decompressedBytes.erase(decompressedBytes.begin(), decompressedBytes.begin() + bytes_per_snp);
      std::memcpy(snpLine.data(), decompressedBytes.data(), bytes_per_snp);
#pragma omp parallel for
      for (uint32_t i = 0; i < nsamples; i++) {
        G(i, lastSNP) = snpLine[i];  // TODO:  maybe do scaling
      }
      lastSNP++;
    }
  }

  if (lastSNP < actual_block_size) {  // keep reading blocks

    size_t read;
    while (!isEmpty) {
      read = res.readFrom(ifs);
      isEmpty = read == 0;
      ZSTD_inBuffer input = {res.getRawInData(), read, 0};
      while (input.pos < input.size) {
        ZSTD_outBuffer output = {res.getRawOutData(), res.getToWrite(), 0};
        lastRet = ctx(input, output);  // perform decompression
        if (ZSTD_isError(lastRet)) cao.error("Error: ZSTD decompression failed. lastRet:", lastRet);
        std::vector<char> _tmp((char *)res.getRawOutData(), (char *)res.getRawOutData() + output.pos);
        decompressedBytes.insert(decompressedBytes.end(), _tmp.begin(), _tmp.end());
        if (no_magic) {
          while (lastSNP < actual_block_size && decompressedBytes.size() >= bytes_per_snp) {  // parse values
            decompressedBytes.erase(decompressedBytes.begin(), decompressedBytes.begin() + bytes_per_snp);
            std::memcpy(snpLine.data(), decompressedBytes.data(), bytes_per_snp);
#pragma omp parallel for
            for (uint32_t i = 0; i < nsamples; i++) {
              G(i, lastSNP) = snpLine[i];  // TODO:  maybe do scaling
            }
            lastSNP++;
          }
        } else {
          while (lastSNP < actual_block_size &&
                 decompressedBytes.size() >= bytes_per_snp + 8) {  // parse values
            decompressedBytes.erase(decompressedBytes.begin() + 8,
                                    decompressedBytes.begin() + bytes_per_snp + 8);
            std::memcpy(snpLine.data(), decompressedBytes.data(), bytes_per_snp);
#pragma omp parallel for
            for (uint32_t i = 0; i < nsamples; i++) {
              G(i, lastSNP) = snpLine[i];  // TODO:  maybe do scaling
            }
            lastSNP++;
          }
          // remove the magic now
          decompressedBytes.erase(decompressedBytes.begin(), decompressedBytes.begin() + 8);
          no_magic = true;
        }
      }
      if (lastSNP >= actual_block_size) break;
    }
  }

  if (lastSNP != actual_block_size) cao.error("something wrong when read_block_initial");
  G.rowwise() -= G.colwise().mean();  // do centering

  return;
}
