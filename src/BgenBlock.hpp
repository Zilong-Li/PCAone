/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/BgenBlock.hpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_BGENBLOCK_HPP
#define PCAONE_BGENBLOCK_HPP

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Prefetch.hpp"

namespace PCAone {

// Reads the minor allele dosages of BGEN variants, many at once.
//
// The bgen library reads a variant through one std::ifstream: it seeks to the
// variant, reads it, inflates it and decodes it, one variant after another,
// while the PCA waits on a single core. Here the offsets of all the variants
// are found once, by a pass over their headers, and then each OpenMP thread
// reads (pread), inflates and decodes variants of its own, in any order. That
// also lets the out-of-core winSVD read its permutation from the input itself,
// as PGEN does, instead of writing a shuffled copy of the file.
//
// The dosages are those of bgen's Variant::minor_allele_dosage(): the dosage
// of the first allele, swapped to 2 - dosage when bgen's sampling estimate
// calls the first allele the major one, and NaN for a missing sample.
// Biallelic variants only; layouts 1 and 2; no, zlib or zstd compression.
class BgenBlockReader {
 public:
  // layout, compression and nsamples from the BGEN header; first: offset of the
  // first variant; nthreads: decoding threads, by OpenMP thread number
  BgenBlockReader(const std::string& path,
                  int layout,
                  int compression,
                  uint32_t nsamples,
                  uint64_t first,
                  uint64_t nvariants,
                  int nthreads);
  ~BgenBlockReader();
  BgenBlockReader(const BgenBlockReader&) = delete;
  BgenBlockReader& operator=(const BgenBlockReader&) = delete;

  uint64_t variant_ct() const { return offsets_.size() - 1; }
  // bytes of the variants [first, first + count) in the file
  uint64_t bytes(uint64_t first, uint64_t count) const { return offsets_[first + count] - offsets_[first]; }

  // The minor allele dosages of variant vidx into dose[nsamples]. Thread-safe
  // for distinct thr.
  void minor_dosage(int thr, uint64_t vidx, float* dose);

  // Ask the kernel to read these variants (any order) into the page cache in
  // the background, e.g. the next block.
  void prefetch(const std::vector<uint32_t>& variants);
  void cancel();

 private:
  struct Thread;
  void index(uint64_t first, uint64_t nvariants);
  void read(uint64_t offset, uint64_t len, std::vector<unsigned char>& dst) const;
  const unsigned char* genotypes(Thread& t, uint64_t vidx, uint32_t& len);

  int fd_ = -1;
  int layout_, compression_;
  uint32_t nsamples_;
  std::vector<uint64_t> offsets_;  // of each variant, and the end of the last
  std::vector<std::unique_ptr<Thread>> threads_;
  std::unique_ptr<ReadAhead> readahead_;
  std::string path_;
};

}  // namespace PCAone

#endif  // PCAONE_BGENBLOCK_HPP
