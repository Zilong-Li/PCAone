/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/PgenBlock.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "PgenBlock.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <unistd.h>

#include <cerrno>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "Prefetch.hpp"
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
  plink2::PreinitPgfi(&pgfi_);
  // The destructor does not run when a constructor throws, and FilePgen falls
  // back to PgenReader then: give back what was acquired so far (pgenlib's
  // files, the descriptor, the mapping, the readers and their buffers).
  try {
    init(pgen, nthreads);
  } catch (...) {
    release();
    throw;
  }
}

PgenBlockReader::~PgenBlockReader() { release(); }

void PgenBlockReader::init(const std::string& pgen, int nthreads) {
  using namespace plink2;
  char errstr[kPglErrstrBufBlen];
  PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_cachelines = 0;
  if (PgfiInitPhase1(pgen.c_str(), nullptr, UINT32_MAX, nsamples_, &header_ctrl, &pgfi_, &pgfi_cachelines, errstr) !=
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
  choose_probe(pgen, st);

  const uintptr_t genovec_bytes = DivUp(nsamples_, kNypsPerVec) * kBytesPerVec;
  const uintptr_t bitvec_bytes = DivUp(nsamples_, kBitsPerVec) * kBytesPerVec;
  const uintptr_t dosage_bytes = DivUp(nsamples_, 2 * kInt32PerVec) * kBytesPerVec;
  const size_t readers = (size_t)std::max(1, nthreads);
  // reserved, so that pushing a pointer cannot throw and lose it
  for (auto* v : {&genovec_, &dosage_present_, &sample_include_}) v->reserve(readers);
  pgr_.reserve(readers);
  pgr_alloc_.reserve(readers);
  pssi_.reserve(readers);
  dosage_main_.reserve(readers);
  for (size_t t = 0; t < readers; ++t) {
    auto* pgr = static_cast<PgenReader*>(std::malloc(sizeof(PgenReader)));
    if (!pgr) throw std::bad_alloc();
    PreinitPgr(pgr);  // no file: CleanupPgr() is a no-op until PgrInit() opens one
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

// Everything acquired, whether construction got to the end or not. Each piece
// starts out empty (null, -1, no file) and is reset, so this can run twice.
void PgenBlockReader::release() noexcept {
  try {
    cancel();
  } catch (...) {
  }
  plink2::PglErr err = plink2::kPglRetSuccess;
  for (auto* pgr : pgr_) {
    plink2::CleanupPgr(pgr, &err);  // closes the reader's file, if PgrInit() opened one
    std::free(pgr);
  }
  pgr_.clear();
  for (auto* p : pgr_alloc_) plink2::aligned_free(p);
  pgr_alloc_.clear();
  for (auto* v : {&genovec_, &dosage_present_, &sample_include_}) {
    for (auto* p : *v) plink2::aligned_free(p);
    v->clear();
  }
  for (auto* p : dosage_main_) plink2::aligned_free(p);
  dosage_main_.clear();
  pssi_.clear();
  plink2::CleanupPgfi(&pgfi_, &err);  // the file of PgfiInitPhase1(), unless a reader took it over
  if (pgfi_alloc_) plink2::aligned_free(pgfi_alloc_);
  pgfi_alloc_ = nullptr;
  if (map_) ::munmap(map_, file_bytes_);
  map_ = nullptr;
  if (fd_ >= 0) ::close(fd_);
  fd_ = -1;
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

// How to tell whether a page is in the page cache. Linux: a one-byte read that
// must not wait (preadv2 with RWF_NOWAIT), checked on the first page, read just
// before; file systems that do not support it (FUSE, virtiofs) fail the check.
// Then mincore() of the mapped file, the only way on macOS. Since Linux 5.0,
// mincore() reports every page as cached for a file the caller neither owns
// nor may write, such as shared read-only data, so it is not used for those,
// and every block is requested.
void PgenBlockReader::choose_probe(const std::string& pgen, const struct stat& st) {
  if (!file_bytes_) return;
  char byte;
  if (::pread(fd_, &byte, 1, 0) != 1) return;  // now in the page cache
#if defined(__linux__) && defined(RWF_NOWAIT)
  struct iovec v = {&byte, 1};
  if (::preadv2(fd_, &v, 1, 0, RWF_NOWAIT) == 1) {
    probe_ = Probe::NoWait;
    return;
  }
  const bool reliable = ::geteuid() == 0 || st.st_uid == ::geteuid() || ::access(pgen.c_str(), W_OK) == 0;
  if (!reliable) return;
#else
  (void)pgen;
  (void)st;
#endif
  void* m = ::mmap(nullptr, file_bytes_, PROT_READ, MAP_SHARED, fd_, 0);  // never read, only asked
  if (m == MAP_FAILED) return;
  map_ = static_cast<unsigned char*>(m);
  probe_ = Probe::Mincore;
}

bool PgenBlockReader::page_cached(uint64_t offset) const {
  switch (probe_) {
#if defined(__linux__) && defined(RWF_NOWAIT)
    case Probe::NoWait: {
      char byte;
      struct iovec v = {&byte, 1};
      return ::preadv2(fd_, &v, 1, (off_t)offset, RWF_NOWAIT) == 1;
    }
#endif
    case Probe::Mincore: {
#ifdef __APPLE__
      char in = 0;
#else
      unsigned char in = 0;
#endif
      const uint64_t page = offset / page_ * page_;
      return ::mincore(map_ + page, std::min<uint64_t>(page_, file_bytes_ - page), &in) == 0 && (in & 1);
    }
    default:
      return false;  // unknown: request it
  }
}

// Whether the page cache holds the records of these variants already, by all
// pages of up to 32 of them spread over the block. Cheap enough for the calling
// thread, so a block in the page cache costs neither a thread nor requests.
bool PgenBlockReader::cached(const std::vector<uint32_t>& variants) const {
  if (probe_ == Probe::None || variants.empty()) return false;
  const size_t probes = std::min<size_t>(32, variants.size());
  size_t hits = 0;
  for (size_t i = 0; i < probes; ++i) {
    const uint32_t v = variants[i * variants.size() / probes];
    const uint64_t a = plink2::GetPgfiLdbaseFpos(&pgfi_, v) / page_ * page_;
    const uint64_t b = std::min<uint64_t>(plink2::GetPgfiFpos(&pgfi_, v + 1), file_bytes_);
    bool all = true;
    for (uint64_t p = a; all && p < b; p += page_) all = page_cached(p);
    hits += all;
  }
  return hits * 10 >= probes * 9;
}

void PgenBlockReader::willneed(const std::vector<std::pair<uint64_t, uint64_t>>& runs, const std::atomic<bool>& stop) {
  ++requested_;
  // asynchronous: the kernel reads the pages into the page cache, in file order
  for (const auto& [a, b] : runs) {
    if (stop.load(std::memory_order_relaxed)) return;
    detail::willneed(fd_, a, b - a);
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
    ++hits_;  // seen, and requested if needed, while the previous block was used
  } else {
    cancel();
    if (!cached(variants)) willneed(pages(variants), never_stop_);
  }
  ahead_.clear();
}

void PgenBlockReader::prefetch(const std::vector<uint32_t>& variants) {
  cancel();
  ahead_ = variants;
  if (cached(variants)) return;
  // the runs are worked out in the background too: sorting them is not free
  pending_ = std::async(std::launch::async, [this, v = variants] { willneed(pages(v), stop_); });
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
