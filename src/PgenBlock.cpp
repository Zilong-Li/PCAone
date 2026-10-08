/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PgenBlock.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "PgenBlock.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "pgenlib/pgenlib_ffi_support.h"

namespace PCAone {

namespace {
// ALT counts of the hard calls, as PgenReader returns them
alignas(16) const double kAltCountPairs[32] = PAIR_TABLE16(0.0, 1.0, 2.0, -3.0);

template <typename T>
T* aligned_buffer(uintptr_t bytes) {
  unsigned char* p = nullptr;
  if (plink2::cachealigned_malloc(std::max<uintptr_t>(bytes, plink2::kCacheline), &p)) throw std::bad_alloc();
  std::memset(p, 0, bytes);
  return reinterpret_cast<T*>(p);
}

[[noreturn]] void fail(const char* what, int code) {
  // as PgenReader: an error while decoding inside an OpenMP region cannot be thrown
  std::fprintf(stderr, "PGEN %s error %d\n", what, code);
  std::exit(EXIT_FAILURE);
}
}  // namespace

PgenBlockReader::PgenBlockReader(const std::string& pgen, uint32_t nsamples, int nthreads)
    : nsamples_(nsamples) {
  using namespace plink2;
  char errstr[kPglErrstrBufBlen];
  PreinitPgfi(&pgfi_);
  PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_cachelines = 0;
  if (PgfiInitPhase1(pgen.c_str(), nullptr, UINT32_MAX, nsamples, &header_ctrl, &pgfi_, &pgfi_cachelines, errstr) !=
      kPglRetSuccess)
    throw std::runtime_error(&errstr[7]);
  // the limits of PgenReader, which FilePgen uses for everything else
  if (header_ctrl & 0x30) throw std::runtime_error("PGEN with stored allele counts (multiallelic) is not supported");
  pgfi_.max_allele_ct = 2;
  if ((header_ctrl & 0xc0) == 0xc0) {
    nonref_flags_.assign(DivUp(pgfi_.raw_variant_ct, kBitsPerWord) + 1, 0);
    pgfi_.nonref_flags = nonref_flags_.data();
  }
  if (pgfi_cachelines && cachealigned_malloc(pgfi_cachelines * kCacheline, &pgfi_alloc_)) throw std::bad_alloc();
  uint32_t max_vrec_width = 0;
  uintptr_t pgr_cachelines = 0;
  // per-variant fread() mode, as PgenReader
  if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, pgfi_.raw_variant_ct, &max_vrec_width, &pgfi_, pgfi_alloc_,
                     &pgr_cachelines, errstr) != kPglRetSuccess)
    throw std::runtime_error(&errstr[7]);
  if (pgfi_.gflags & 4)
    throw std::runtime_error("PGEN with multiallelic variants and phase or dosages is not supported");

  fd_ = ::open(pgen.c_str(), O_RDONLY);
  if (fd_ < 0) throw std::runtime_error("cannot open " + pgen);
  struct stat st;
  if (::fstat(fd_, &st) == 0) file_bytes_ = (uint64_t)st.st_size;
  page_ = (uint64_t)::sysconf(_SC_PAGESIZE);

  const uintptr_t genovec_bytes = DivUp(nsamples, kNypsPerVec) * kBytesPerVec;
  const uintptr_t bitvec_bytes = DivUp(nsamples, kBitsPerVec) * kBytesPerVec;
  const uintptr_t dosage_bytes = DivUp(nsamples, 2 * kInt32PerVec) * kBytesPerVec;
  for (int t = 0; t < std::max(1, nthreads); ++t) {
    auto* pgr = static_cast<PgenReader*>(std::malloc(sizeof(PgenReader)));
    if (!pgr) throw std::bad_alloc();
    PreinitPgr(pgr);
    PgrSetFreadBuf(nullptr, pgr);
    pgr_.push_back(pgr);
    pgr_alloc_.push_back(aligned_buffer<unsigned char>(pgr_cachelines * kCacheline));
    // the first reader takes over the file pgfi has open, the others open their own
    const PglErr err = PgrInit(pgen.c_str(), max_vrec_width, &pgfi_, pgr, pgr_alloc_.back());
    if (err != kPglRetSuccess) throw std::runtime_error("PgrInit() error " + std::to_string((int)err));
    pssi_.emplace_back();
    PgrClearSampleSubsetIndex(pgr, &pssi_.back());
    genovec_.push_back(aligned_buffer<uintptr_t>(genovec_bytes));
    dosage_present_.push_back(aligned_buffer<uintptr_t>(bitvec_bytes));
    dosage_main_.push_back(aligned_buffer<uint16_t>(dosage_bytes));
    // all samples; pgenlib does not subset when sample_ct is the file's
    sample_include_.push_back(aligned_buffer<uintptr_t>(bitvec_bytes));
    std::memset(sample_include_.back(), 0xff, bitvec_bytes);
  }
}

PgenBlockReader::~PgenBlockReader() {
  cancel();
  plink2::PglErr err = plink2::kPglRetSuccess;
  for (auto* pgr : pgr_) {
    plink2::CleanupPgr(pgr, &err);
    std::free(pgr);
  }
  for (auto* p : pgr_alloc_) plink2::aligned_free(p);
  for (auto* p : genovec_) plink2::aligned_free(p);
  for (auto* p : dosage_present_) plink2::aligned_free(p);
  for (auto* p : dosage_main_) plink2::aligned_free(p);
  for (auto* p : sample_include_) plink2::aligned_free(p);
  plink2::CleanupPgfi(&pgfi_, &err);
  if (pgfi_alloc_) plink2::aligned_free(pgfi_alloc_);
  if (fd_ >= 0) ::close(fd_);
}

// pgenlib reads a variant's record, from the record of its LD base when it is
// LD-compressed (PgfiMultiread loads the same), so the pages of those bytes.
std::vector<std::pair<uint64_t, uint64_t>> PgenBlockReader::pages(const std::vector<uint32_t>& variants) const {
  std::vector<std::pair<uint64_t, uint64_t>> r;
  r.reserve(variants.size());
  for (const uint32_t v : variants) {
    const uint64_t a = plink2::GetPgfiLdbaseFpos(&pgfi_, v), b = plink2::GetPgfiFpos(&pgfi_, v + 1);
    r.emplace_back(a / page_ * page_, (b + page_ - 1) / page_ * page_);
  }
  std::sort(r.begin(), r.end());
  std::vector<std::pair<uint64_t, uint64_t>> runs;
  for (const auto& [a, b] : r) {
    if (!runs.empty() && a <= runs.back().second)
      runs.back().second = std::max(runs.back().second, b);
    else
      runs.emplace_back(a, b);
  }
  return runs;
}

void PgenBlockReader::willneed(const std::vector<std::pair<uint64_t, uint64_t>>& runs, const std::atomic<bool>& stop) {
  // asynchronous: the kernel reads the pages into the page cache, in file order
  for (const auto& [a, b] : runs) {
    if (stop.load(std::memory_order_relaxed)) return;
#if defined(POSIX_FADV_WILLNEED)
    ::posix_fadvise(fd_, (off_t)a, (off_t)(b - a), POSIX_FADV_WILLNEED);
#elif defined(F_RDADVISE)
    struct radvisory ra;
    ra.ra_offset = (off_t)a;
    ra.ra_count = (int)std::min<uint64_t>(b - a, INT32_MAX);
    ::fcntl(fd_, F_RDADVISE, &ra);
#endif
  }
}

void PgenBlockReader::cancel() {
  if (!pending_.valid()) return;
  stop_ = true;
  pending_.get();
  stop_ = false;
}

void PgenBlockReader::load(const std::vector<uint32_t>& variants) {
  if (!ahead_.empty() && ahead_ == variants) {
    ++hits_;  // requested while the previous block was used
  } else {
    cancel();
    willneed(pages(variants), never_stop_);
  }
  ahead_.clear();
}

void PgenBlockReader::prefetch(const std::vector<uint32_t>& variants) {
  cancel();
  ahead_ = variants;
  pending_ = std::async(std::launch::async, [this, runs = pages(variants)] { willneed(runs, stop_); });
}

const uintptr_t* PgenBlockReader::genovec(int thr, uint32_t vidx) {
  const plink2::PglErr err =
      plink2::PgrGet1(sample_include_[thr], pssi_[thr], nsamples_, vidx, 1, pgr_[thr], genovec_[thr]);
  if (err != plink2::kPglRetSuccess) fail("PgrGet1()", (int)err);
  return genovec_[thr];
}

void PgenBlockReader::hardcalls(int thr, uint32_t vidx, double* buf) {
  plink2::GenoarrLookup16x8bx2(genovec(thr, vidx), kAltCountPairs, nsamples_, buf);
}

const uintptr_t* PgenBlockReader::genovec_dosages(int thr, uint32_t vidx, const uintptr_t** present,
                                                  const uint16_t** main, uint32_t* dosage_ct) {
  const plink2::PglErr err = plink2::PgrGet1D(sample_include_[thr], pssi_[thr], nsamples_, vidx, 1, pgr_[thr],
                                              genovec_[thr], dosage_present_[thr], dosage_main_[thr], dosage_ct);
  if (err != plink2::kPglRetSuccess) fail("PgrGet1D()", (int)err);
  *present = dosage_present_[thr];
  *main = dosage_main_[thr];
  return genovec_[thr];
}

void PgenBlockReader::dosages(int thr, uint32_t vidx, double* buf) {
  const uintptr_t* present;
  const uint16_t* main;
  uint32_t dosage_ct = 0;
  const uintptr_t* geno = genovec_dosages(thr, vidx, &present, &main, &dosage_ct);
  plink2::Dosage16ToDoubles(kAltCountPairs, geno, present, main, nsamples_, dosage_ct, buf);
}

}  // namespace PCAone
