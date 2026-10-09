/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PgenBlock.hpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_PGENBLOCK_HPP
#define PCAONE_PGENBLOCK_HPP

#include <atomic>
#include <cstdint>
#include <future>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include "pgenlib/include/pgenlib_read.h"

namespace PCAone {

// Out-of-core PGEN reads, with the records of the next block read ahead in
// the background. pgenlib itself is used as it is.
//
// Like PgenReader (pgenlibr), one pgenlib reader per thread reads its variants
// with fread(). The logical permutation scatters a block over the whole file,
// so on a cold disk every variant is a seek on the decoding thread, after the
// computation of the previous block. Here, while a block is decoded and
// multiplied, the records of the next one (with the LD base of an
// LD-compressed variant) are requested from the kernel in file order
// (POSIX_FADV_WILLNEED, F_RDADVISE on macOS), so that reading them overlaps the
// computation and the disk sees them sorted. Blocks already in the page cache
// (by a sample of their pages) are not requested. Decoding is pgenlib's
// own, so the values are exactly those of PgenReader, and the 2-bit calls are
// available to decode straight into G.
class PgenBlockReader {
  enum class Probe { None, NoWait, Mincore };

 public:
  // nthreads: decoding threads, by OpenMP thread number
  PgenBlockReader(const std::string& pgen, uint32_t nsamples, int nthreads);
  ~PgenBlockReader();
  PgenBlockReader(const PgenBlockReader&) = delete;
  PgenBlockReader& operator=(const PgenBlockReader&) = delete;

  uint32_t variant_ct() const { return pgfi_.raw_variant_ct; }

  // About to decode these variants (file indices, any order). Unless they were
  // requested ahead, they are requested now, sorted.
  void load(const std::vector<uint32_t>& variants);
  // Request these variants in the background, e.g. the next block.
  void prefetch(const std::vector<uint32_t>& variants);
  // Stop a request in the background, e.g. at the end of a run.
  void cancel();

  // Thread-safe for distinct thr. The same values as PgenReader::ReadHardcalls
  // and PgenReader::Read with allele_idx = 1: ALT counts 0, 1, 2 and -3 for a
  // missing call, or dosages in [0, 2].
  void hardcalls(int thr, uint32_t vidx, double* buf);
  void dosages(int thr, uint32_t vidx, double* buf);
  // ALT counts at 2 bits per sample, 3 = missing
  const uintptr_t* genovec(int thr, uint32_t vidx);
  // the same with the dosages: the dosage_ct samples set in *present have the
  // dosage (*main)[i] / 16384 instead of their call
  const uintptr_t* genovec_dosages(int thr, uint32_t vidx, const uintptr_t** present, const uint16_t** main,
                                   uint32_t* dosage_ct);

  uint64_t predicted() const { return hits_; }
  uint64_t requested() const { return requested_; }  // blocks not found in the page cache
  const char* probe_name() const {
    return probe_ == Probe::NoWait ? "RWF_NOWAIT" : probe_ == Probe::Mincore ? "mincore" : "none";
  }

 private:
  void init(const std::string& pgen, int nthreads);
  void release() noexcept;
  // the page runs pgenlib reads for these variants, sorted and merged
  std::vector<std::pair<uint64_t, uint64_t>> pages(const std::vector<uint32_t>& variants) const;
  void choose_probe(const std::string& pgen, const struct stat& st);
  bool page_cached(uint64_t offset) const;
  bool cached(const std::vector<uint32_t>& variants) const;
  void willneed(const std::vector<std::pair<uint64_t, uint64_t>>& runs, const std::atomic<bool>& stop);

  plink2::PgenFileInfo pgfi_;
  unsigned char* pgfi_alloc_ = nullptr;
  std::vector<uintptr_t> nonref_flags_;
  uint32_t nsamples_ = 0;
  int fd_ = -1;                    // for the read-ahead requests
  Probe probe_ = Probe::None;       // how page_cached() asks the page cache
  unsigned char* map_ = nullptr;   // the file, for mincore() only; never read
  uint64_t file_bytes_ = 0, page_ = 4096;
  std::vector<plink2::PgenReader*> pgr_;
  std::vector<unsigned char*> pgr_alloc_;
  std::vector<plink2::PgrSampleSubsetIndex> pssi_;
  std::vector<uintptr_t*> genovec_, dosage_present_, sample_include_;
  std::vector<uint16_t*> dosage_main_;
  std::vector<uint32_t> ahead_;                        // the variants requested ahead
  std::future<void> pending_;
  std::atomic<bool> stop_{false};
  const std::atomic<bool> never_stop_{false};
  uint64_t hits_ = 0;
  std::atomic<uint64_t> requested_{0};
};

}  // namespace PCAone

#endif  // PCAONE_PGENBLOCK_HPP
