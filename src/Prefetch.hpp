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
#include <thread>
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
// Memory: two blocks of raw records.
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

    // read ahead only along a sequential pass (or the start of one)
    const bool sequential = first == 0 || first == last_end_;
    last_end_ = first + count;
    if (first == 0) first_count_ = count;
    if (ahead_ && sequential && count < n_) {
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

 private:
  void read_into(std::vector<unsigned char>& buf, uint64_t first, uint64_t count, const std::atomic<bool>& stop) {
    const uint64_t len = count * width_;
    if (buf.size() < len) buf.resize(len);
    if (detail::pread_all(fd_, buf.data(), len, offset_ + first * width_, stop)) bytes_ += len;
  }

  int fd_;
  uint64_t offset_, width_, n_;
  double min_share_;
  bool ahead_ = true, returned_ = false;
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

// Reads byte ranges of a file in the background and throws the bytes away, so
// that a library which reads the file itself (BGEN, PGEN) then finds them in
// the page cache. A hint for later reads, it never changes what is read.
class ReadAhead {
 public:
  explicit ReadAhead(const std::string& path, int threads = 1)
      : fd_(detail::open_for_streaming(path)), threads_(std::max(1, threads)) {}

  ~ReadAhead() {
    cancel();
    if (fd_ >= 0) ::close(fd_);
  }

  ReadAhead(const ReadAhead&) = delete;
  ReadAhead& operator=(const ReadAhead&) = delete;

  // Read the ranges (offset, length) in the background; they are sorted and
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
    pending_ = std::async(std::launch::async, [this, merged = std::move(merged)] { run(merged); });
  }

  void hint(uint64_t offset, uint64_t length) { hint({{offset, length}}); }

  // Wait until the ranges of the last hint are read.
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

 private:
  void run(const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
    // a few readers keep a disk queue busy when the ranges are scattered
    const int nt = (int)std::min<size_t>(threads_, ranges.size());
    std::atomic<size_t> next{0};
    auto worker = [&] {
      std::vector<unsigned char> sink(1ULL << 20);
      for (size_t i; (i = next++) < ranges.size() && !stop_;) {
        uint64_t off = ranges[i].first;
        const uint64_t end = ranges[i].second;
        while (off < end && !stop_) {
          const uint64_t len = std::min<uint64_t>(sink.size(), end - off);
          const ssize_t got = ::pread(fd_, sink.data(), len, (off_t)off);
          if (got <= 0) {
            if (got < 0 && errno == EINTR) continue;
            return;  // only a hint: the library reports real read errors
          }
          off += (uint64_t)got;
        }
      }
    };
    std::vector<std::thread> pool;
    for (int t = 1; t < nt; ++t) pool.emplace_back(worker);
    worker();
    for (auto& t : pool) t.join();
  }

  int fd_;
  int threads_;
  std::future<void> pending_;
  std::atomic<bool> stop_{false};
};

}  // namespace PCAone
#endif  // PCAONE_PREFETCH_HPP
