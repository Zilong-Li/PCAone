/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Prefetch.hpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_PREFETCH_HPP
#define PCAONE_PREFETCH_HPP

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <future>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace PCAone {

namespace detail {
// pread() [offset, offset + len) into dst in chunks, giving up early when stop
// is set. Returns false if it stopped, throws on a read error or a short file.
inline bool pread_all(int fd, unsigned char* dst, uint64_t len, uint64_t offset, const std::atomic<bool>& stop) {
  constexpr uint64_t chunk = 8ULL << 20;  // check for a cancel every 8 MiB
  uint64_t done = 0;
  while (done < len) {
    if (stop.load(std::memory_order_relaxed)) return false;
    const uint64_t want = std::min(chunk, len - done);
    const ssize_t got = ::pread(fd, dst + done, want, (off_t)(offset + done));
    if (got < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error(std::string("read error: ") + std::strerror(errno));
    }
    if (got == 0) throw std::runtime_error("unexpected end of file");
    done += (uint64_t)got;
  }
  return true;
}

inline int open_for_streaming(const std::string& path) {
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
#if defined(POSIX_FADV_SEQUENTIAL)
  ::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);  // larger kernel readahead
#endif
  return fd;
}

#if defined(POSIX_FADV_WILLNEED) || defined(F_RDADVISE)
constexpr bool kCanAdvise = true;
#else
constexpr bool kCanAdvise = false;
#endif

// Ask the kernel to read [offset, offset + len) into the page cache. It queues
// the reads and returns; the process holds no copy of the bytes. Only a hint:
// a file system that ignores it (or an error) leaves the reads as they were.
inline void willneed(int fd, uint64_t offset, uint64_t len) {
#if defined(POSIX_FADV_WILLNEED)
  ::posix_fadvise(fd, (off_t)offset, (off_t)len, POSIX_FADV_WILLNEED);
#elif defined(F_RDADVISE)
  while (len) {  // ra_count is an int
    struct radvisory ra;
    ra.ra_offset = (off_t)offset;
    ra.ra_count = (int)std::min<uint64_t>(len, 1ULL << 30);
    ::fcntl(fd, F_RDADVISE, &ra);
    offset += (uint64_t)ra.ra_count;
    len -= (uint64_t)ra.ra_count;
  }
#else
  (void)fd;
  (void)offset;
  (void)len;
#endif
}
}  // namespace detail

// Reads the fixed-width records of a file (BED SNPs, float columns of the binary
// format) block by block, and the next block in the background when that pays.
//
// The out-of-core PCA goes through the blocks in order, pass after pass. While
// the caller decodes and multiplies block i, a second thread reads block i + 1,
// or the first block again after the last one, so the reads overlap the
// computation instead of alternating with it. The bytes are the same as those
// of a plain read, so the results do not change. A block that was not
// predicted (e.g. LD reading blocks out of order) is read on the spot, and
// nothing is read ahead until the access is sequential again.
//
// A read in the background is a thread that wants a core while the PCA uses
// all of them, so it is worth it only when the reads take time. A block whose
// read took less than `min_share` of the computation between two calls (the
// file is in the page cache) has the next one read on the spot; a slower read
// switches the reading ahead back on.
//
// Memory: one block of raw records more than a plain read, only while reading
// ahead (for BED, 1/32 of the decoded block); freed when reading on the spot.
class RecordPrefetcher {
 public:
  RecordPrefetcher(const std::string& path, uint64_t data_offset, uint64_t record_bytes, uint64_t nrecords,
                   double min_share = 0.05)
      : fd_(detail::open_for_streaming(path)), offset_(data_offset), width_(record_bytes), n_(nrecords),
        min_share_(min_share) {}

  ~RecordPrefetcher() {
    cancel();
    if (fd_ >= 0) ::close(fd_);
  }

  RecordPrefetcher(const RecordPrefetcher&) = delete;
  RecordPrefetcher& operator=(const RecordPrefetcher&) = delete;

  // Read record order[r] of the file for record r (a logical permutation),
  // instead of record r. Consecutive file records are read together, and the
  // runs of a block are requested from the kernel first, so that a disk can
  // serve them in its own order.
  void set_order(std::vector<uint32_t> order) {
    cancel();
    for (auto r : order)
      if (r >= n_) throw std::out_of_range("RecordPrefetcher: a record beyond the end of the file");
    order_ = std::move(order);
  }
  // false: never read ahead in the background (--no-prefetch)
  void set_read_ahead(bool on) {
    cancel();
    allow_ahead_ = on;
  }

  // The records [first, first + count), valid until the next call.
  const unsigned char* get(uint64_t first, uint64_t count) {
    using clock = std::chrono::steady_clock;
    if (first + count > n_) throw std::out_of_range("RecordPrefetcher: records beyond the end of the file");
    const auto t0 = clock::now();
    // the caller's work on the previous block, which a read ahead can hide
    const double work = returned_ ? std::chrono::duration<double>(t0 - last_return_).count() : 0.0;
    double read = 0;
    if (pending_.valid() && next_first_ == first && next_count_ == count) {
      pending_.get();  // rethrows a read error
      std::swap(cur_, next_);
      read = ahead_seconds_;
      ++hits_;
    } else {
      cancel();
      read_into(cur_, first, count, never_stop_);
      read = std::chrono::duration<double>(clock::now() - t0).count();
    }
    wait_ += std::chrono::duration<double>(clock::now() - t0).count();
    if (returned_) ahead_ = read >= min_share_ * work;
    if (!ahead_ && next_.capacity()) std::vector<unsigned char>().swap(next_);  // not needed until it pays again

    // read ahead only along a sequential pass (or the start of one)
    const bool sequential = first == 0 || first == last_end_;
    last_end_ = first + count;
    if (first == 0) first_count_ = count;
    if (allow_ahead_ && ahead_ && sequential && count < n_) {
      uint64_t nf = first + count, nc = std::min(count, n_ - nf);
      if (nf >= n_) {  // the first block of the next pass
        nf = 0;
        nc = std::min(first_count_ ? first_count_ : count, n_);
      }
      next_first_ = nf;
      next_count_ = nc;
      pending_ = std::async(std::launch::async, [this, nf, nc] {
        const auto a = clock::now();
        read_into(next_, nf, nc, stop_);
        ahead_seconds_ = std::chrono::duration<double>(clock::now() - a).count();
      });
    }
    returned_ = true;
    last_return_ = clock::now();
    return cur_.data();
  }

  // Stop and wait for a read in the background, e.g. at the end of a run.
  void cancel() {
    if (!pending_.valid()) return;
    stop_ = true;
    try {
      pending_.get();
    } catch (...) {
    }
    stop_ = false;
  }

  double wait_seconds() const { return wait_; }   // time the caller was blocked on reads
  uint64_t predicted() const { return hits_; }    // blocks that came from the background
  uint64_t bytes_read() const { return bytes_; }  // call cancel() first
  uint64_t buffer_bytes() const { return cur_.capacity() + next_.capacity(); }

 private:
  void read_into(std::vector<unsigned char>& buf, uint64_t first, uint64_t count, const std::atomic<bool>& stop) {
    const uint64_t len = count * width_;
    if (buf.size() < len) buf.resize(len);
    if (order_.empty()) {
      if (detail::pread_all(fd_, buf.data(), len, offset_ + first * width_, stop)) bytes_ += len;
      return;
    }
    // runs of consecutive file records: (first logical record, count)
    std::vector<std::pair<uint64_t, uint64_t>> runs;
    for (uint64_t i = first, end = first + count; i < end;) {
      uint64_t j = i + 1;
      while (j < end && order_[j] == order_[j - 1] + 1) ++j;
      runs.emplace_back(i, j - i);
      i = j;
    }
    if (runs.size() > 1)
      for (const auto& [i, n] : runs) detail::willneed(fd_, offset_ + (uint64_t)order_[i] * width_, n * width_);
    for (const auto& [i, n] : runs) {
      if (!detail::pread_all(fd_, buf.data() + (i - first) * width_, n * width_, offset_ + (uint64_t)order_[i] * width_,
                             stop))
        return;
      bytes_ += n * width_;
    }
  }

  int fd_;
  uint64_t offset_, width_, n_;
  double min_share_;
  bool ahead_ = true, returned_ = false, allow_ahead_ = true;
  std::vector<uint32_t> order_;  // logical record -> file record; empty: the same
  std::chrono::steady_clock::time_point last_return_;
  double ahead_seconds_ = 0;  // written by the reading thread, read after it is joined
  std::vector<unsigned char> cur_, next_;
  std::future<void> pending_;
  uint64_t next_first_ = 0, next_count_ = 0, last_end_ = 0, first_count_ = 0, hits_ = 0;
  std::atomic<bool> stop_{false};
  const std::atomic<bool> never_stop_{false};
  std::atomic<uint64_t> bytes_{0};
  double wait_ = 0;
};

// Asks for byte ranges of a file to be read into the page cache in the
// background, so that a reader that reads the file itself (the binary copy of
// CSV, BGEN) then finds them there. The process holds no copy of the bytes, so
// the memory of -m is unchanged. A hint for later reads, it never changes what
// is read; where the kernel ignores it, the reads are as without it.
class ReadAhead {
 public:
  explicit ReadAhead(const std::string& path) : fd_(::open(path.c_str(), O_RDONLY)) {
    if (fd_ < 0) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
  }

  ~ReadAhead() {
    cancel();
    if (fd_ >= 0) ::close(fd_);
  }

  ReadAhead(const ReadAhead&) = delete;
  ReadAhead& operator=(const ReadAhead&) = delete;

  // Request the ranges (offset, length) in the background; they are sorted and
  // merged when less than `gap` bytes apart. A hint still running is stopped.
  void hint(std::vector<std::pair<uint64_t, uint64_t>> ranges, uint64_t gap = 4096) {
    cancel();
    if (ranges.empty()) return;
    std::sort(ranges.begin(), ranges.end());
    std::vector<std::pair<uint64_t, uint64_t>> merged;  // (start, end)
    for (const auto& r : ranges) {
      if (!r.second) continue;
      if (!merged.empty() && r.first <= merged.back().second + gap)
        merged.back().second = std::max(merged.back().second, r.first + r.second);
      else
        merged.emplace_back(r.first, r.first + r.second);
    }
    if (merged.empty()) return;
    ++hints_;
    pending_ = std::async(std::launch::async, [this, merged = std::move(merged)] { run(merged); });
  }

  void hint(uint64_t offset, uint64_t length) { hint({{offset, length}}); }

  // Wait until the ranges of the last hint are requested.
  void wait() {
    if (pending_.valid()) pending_.get();
  }

  void cancel() {
    if (!pending_.valid()) return;
    stop_ = true;
    try {
      pending_.get();
    } catch (...) {
    }
    stop_ = false;
  }

  uint64_t hints() const { return hints_; }

 private:
  void run(const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
    // in pieces, so that a cancel does not wait for a whole block to be queued
    constexpr uint64_t piece = 8ULL << 20;
#if !defined(POSIX_FADV_WILLNEED) && !defined(F_RDADVISE)
    std::vector<unsigned char> sink(1ULL << 20);  // no advice: read and discard
#endif
    for (const auto& [a, b] : ranges)
      for (uint64_t off = a; off < b; off += piece) {
        if (stop_.load(std::memory_order_relaxed)) return;
        const uint64_t len = std::min(piece, b - off);
#if defined(POSIX_FADV_WILLNEED) || defined(F_RDADVISE)
        detail::willneed(fd_, off, len);
#else
        for (uint64_t done = 0; done < len;) {
          const ssize_t got = ::pread(fd_, sink.data(), std::min<uint64_t>(sink.size(), len - done), (off_t)(off + done));
          if (got < 0 && errno == EINTR) continue;
          if (got <= 0) return;  // only a hint: the reader reports real read errors
          done += (uint64_t)got;
        }
#endif
      }
  }

  int fd_;
  uint64_t hints_ = 0;
  std::future<void> pending_;
  std::atomic<bool> stop_{false};
};

}  // namespace PCAone
#endif  // PCAONE_PREFETCH_HPP
