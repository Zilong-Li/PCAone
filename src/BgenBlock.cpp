/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/BgenBlock.cpp
 * @author      Zilong Li
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "BgenBlock.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <cerrno>
#include <cmath>
#include <cstring>
#include <stdexcept>

#include "zstd.h"

namespace PCAone {

namespace {
inline uint16_t rd16(const unsigned char* p) {
  uint16_t v;
  std::memcpy(&v, p, sizeof v);
  return v;
}
inline uint32_t rd32(const unsigned char* p) {
  uint32_t v;
  std::memcpy(&v, p, sizeof v);
  return v;
}
inline uint64_t rd64(const unsigned char* p) {
  uint64_t v;
  std::memcpy(&v, p, sizeof v);
  return v;
}

// pread up to len bytes; fewer only at the end of the file
uint64_t pread_some(int fd, unsigned char* dst, uint64_t len, uint64_t offset) {
  uint64_t done = 0;
  while (done < len) {
    const ssize_t got = ::pread(fd, dst + done, len - done, (off_t)(offset + done));
    if (got < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error(std::string("BGEN read error: ") + std::strerror(errno));
    }
    if (got == 0) break;
    done += (uint64_t)got;
  }
  return done;
}

// bgen's minor_certain(): does the confidence interval of freq miss 0.5?
bool minor_certain(double freq, int n_checked, double z) {
  double delta = (z * std::sqrt((freq * (1 - freq)) / n_checked));
  return !((freq - delta < 0.5) & (freq + delta > 0.5));
}

// bgen's Genotypes::find_minor_allele(): 0 if the first allele is the minor
// one, by the dosages of the first allele of evenly spaced samples
int find_minor_allele(const float* dose, uint32_t n_samples) {
  const uint32_t batchsize = 100;
  const uint32_t increment = std::max(n_samples / batchsize, (uint32_t)1);
  double total = 0;
  double freq = 0;
  for (uint32_t idx2 = 0; idx2 < increment; idx2++) {
    for (uint32_t n = idx2; n < n_samples; n += increment) total += dose[n];
    freq = total / (batchsize * (idx2 + 1) * 2);
    if (minor_certain(freq, batchsize * (idx2 + 1), 5.0)) break;
  }
  return freq <= 0.5 ? 0 : 1;
}

constexpr uint64_t kPad = 8;  // the bit reader loads 8 bytes at a time
}  // namespace

struct BgenBlockReader::Thread {
  std::vector<unsigned char> raw, unc;
  std::vector<uint32_t> missing;
  z_stream z{};
  bool zinit = false;
  ZSTD_DCtx* zd = nullptr;
  ~Thread() {
    if (zinit) inflateEnd(&z);
    if (zd) ZSTD_freeDCtx(zd);
  }
};

BgenBlockReader::BgenBlockReader(const std::string& path, int layout, int compression, uint32_t nsamples,
                                 uint64_t first, uint64_t nvariants, int nthreads)
    : layout_(layout), compression_(compression), nsamples_(nsamples), path_(path) {
  if (layout_ != 1 && layout_ != 2) throw std::invalid_argument("unsupported BGEN layout " + std::to_string(layout_));
  if (compression_ < 0 || compression_ > 2)
    throw std::invalid_argument("unsupported BGEN compression " + std::to_string(compression_));
  if (layout_ == 1 && compression_ == 2) throw std::invalid_argument("BGEN layout 1 cannot be zstd compressed");
  fd_ = detail::open_for_streaming(path);
  try {
    index(first, nvariants);
  } catch (...) {
    ::close(fd_);
    throw;
  }
  for (int i = 0; i < std::max(1, nthreads); ++i) threads_.emplace_back(new Thread);
}

BgenBlockReader::~BgenBlockReader() {
  readahead_.reset();
  if (fd_ >= 0) ::close(fd_);
}

void BgenBlockReader::read(uint64_t offset, uint64_t len, std::vector<unsigned char>& dst) const {
  if (dst.size() < len + kPad) dst.resize(len + kPad);
  if (pread_some(fd_, dst.data(), len, offset) != len) throw std::runtime_error("BGEN: unexpected end of file");
  std::memset(dst.data() + len, 0, kPad);
}

// one pass over the variant headers, skipping the genotypes: the offset of
// every variant
void BgenBlockReader::index(uint64_t first, uint64_t nvariants) {
  struct stat st;
  if (::fstat(fd_, &st) != 0) throw std::runtime_error("cannot stat " + path_);
  const uint64_t fsize = (uint64_t)st.st_size;

  // the bytes [at, at + n), from a window of the file read 64 KiB at a time:
  // headers of small variants come many to a window, those of large ones one
  std::vector<unsigned char> buf;
  uint64_t buf_off = 0, buf_len = 0;
  auto need = [&](uint64_t at, uint64_t n) -> const unsigned char* {
    if (at < buf_off || at + n > buf_off + buf_len) {
      const uint64_t want = std::max<uint64_t>(n, 1ULL << 16);
      if (buf.size() < want) buf.resize(want);
      buf_off = at;
      buf_len = pread_some(fd_, buf.data(), want, at);
      if (buf_len < n) throw std::runtime_error("BGEN: a variant runs past the end of the file");
    }
    return buf.data() + (at - buf_off);
  };

  offsets_.resize(nvariants + 1);
  uint64_t pos = first;
  for (uint64_t v = 0; v < nvariants; ++v) {
    if (pos >= fsize) throw std::runtime_error("BGEN: the file ends after " + std::to_string(v) + " variants");
    offsets_[v] = pos;
    uint64_t at = pos;
    if (layout_ == 1) at += 4;  // the number of samples
    for (int f = 0; f < 3; ++f) at += 2 + rd16(need(at, 2));  // variant id, rsid, chromosome
    at += 4;                                                    // position
    uint16_t nalleles = 2;
    if (layout_ == 2) {
      nalleles = rd16(need(at, 2));
      at += 2;
    }
    for (uint16_t a = 0; a < nalleles; ++a) at += 4 + rd32(need(at, 4));
    uint64_t length = (uint64_t)nsamples_ * 6;
    if (!(layout_ == 1 && compression_ == 0)) {
      length = rd32(need(at, 4));
      at += 4;
    }
    pos = at + length;
  }
  if (pos > fsize) throw std::runtime_error("BGEN: the last variant runs past the end of the file");
  offsets_[nvariants] = pos;
}

// read variant vidx into t.raw; returns its genotype data and its length
const unsigned char* BgenBlockReader::genotypes(Thread& t, uint64_t vidx, uint32_t& len) {
  const uint64_t off = offsets_[vidx], n = offsets_[vidx + 1] - off;
  read(off, n, t.raw);
  const unsigned char* p = t.raw.data();
  const unsigned char* end = p + n;
  auto check = [&](uint64_t k) {
    if (p + k > end) throw std::runtime_error("BGEN: malformed variant " + std::to_string(vidx));
  };
  if (layout_ == 1) {
    check(4);
    if (rd32(p) != nsamples_) throw std::runtime_error("BGEN: number of samples doesn't match");
    p += 4;
  }
  for (int f = 0; f < 3; ++f) {
    check(2);
    p += 2 + rd16(p);
  }
  p += 4;
  uint16_t nalleles = 2;
  if (layout_ == 2) {
    check(2);
    nalleles = rd16(p);
    p += 2;
  }
  if (nalleles != 2)
    throw std::runtime_error("BGEN: variant " + std::to_string(vidx) +
                             " is not biallelic; PCAone uses the dosages of biallelic variants only");
  for (uint16_t a = 0; a < nalleles; ++a) {
    check(4);
    p += 4 + rd32(p);
  }
  len = nsamples_ * 6;
  if (!(layout_ == 1 && compression_ == 0)) {
    check(4);
    len = rd32(p);
    p += 4;
  }
  check(len);
  return p;
}

void BgenBlockReader::minor_dosage(int thr, uint64_t vidx, float* dose) {
  Thread& t = *threads_.at(thr);
  uint32_t glen;
  const unsigned char* g = genotypes(t, vidx, glen);

  // inflate
  const unsigned char* d = g;
  uint32_t dlen = glen;
  if (compression_ != 0) {
    const unsigned char* src = g;
    uint32_t slen = glen;
    if (layout_ == 1) {
      dlen = nsamples_ * 6;
    } else {
      if (glen < 4) throw std::runtime_error("BGEN: malformed variant " + std::to_string(vidx));
      dlen = rd32(g);
      src += 4;
      slen -= 4;
    }
    if (t.unc.size() < (uint64_t)dlen + kPad) t.unc.resize((uint64_t)dlen + kPad);
    if (compression_ == 1) {
      if (!t.zinit) {
        if (inflateInit(&t.z) != Z_OK) throw std::runtime_error("zlib: inflateInit failed");
        t.zinit = true;
      } else {
        inflateReset(&t.z);
      }
      t.z.next_in = const_cast<Bytef*>(src);
      t.z.avail_in = slen;
      t.z.next_out = t.unc.data();
      t.z.avail_out = dlen;
      const int r = inflate(&t.z, Z_FINISH);
      if (r != Z_STREAM_END || t.z.total_out != dlen)
        throw std::runtime_error("BGEN: zlib decompression gave data of wrong length at variant " +
                                 std::to_string(vidx));
    } else {
      if (!t.zd) t.zd = ZSTD_createDCtx();
      const size_t got = ZSTD_decompressDCtx(t.zd, t.unc.data(), dlen, src, slen);
      if (ZSTD_isError(got) || got != dlen)
        throw std::runtime_error("BGEN: zstd decompression gave data of wrong length at variant " +
                                 std::to_string(vidx));
    }
    std::memset(t.unc.data() + dlen, 0, kPad);
    d = t.unc.data();
  }

  // the dosage of the first allele; missing samples are counted as 0 here, as
  // bgen does, until the minor allele is chosen. bgen multiplies by 1 / max;
  // the division makes a certain genotype exactly 0, 1 or 2 at any bit depth.
  const uint32_t N = nsamples_;
  t.missing.clear();
  if (layout_ == 1) {
    if (dlen < (uint64_t)N * 6) throw std::runtime_error("BGEN: malformed variant " + std::to_string(vidx));
    const float factor = 1.0f / 32768;
    for (uint32_t n = 0; n < N; ++n) {
      const uint32_t hom = rd16(d + 6 * n), het = rd16(d + 6 * n + 2), alt = rd16(d + 6 * n + 4);
      dose[n] = (hom * 2 + het) * factor;
      if ((hom == 0) & (het == 0) & (alt == 0)) t.missing.push_back(n);
    }
  } else {
    if (dlen < 10 + (uint64_t)N) throw std::runtime_error("BGEN: malformed variant " + std::to_string(vidx));
    if (rd32(d) != N) throw std::runtime_error("BGEN: number of samples doesn't match!");
    if (rd16(d + 4) != 2) throw std::runtime_error("BGEN: number of alleles doesn't match!");
    const uint8_t min_ploidy = d[6], max_ploidy = d[7];
    const unsigned char* ploidy = d + 8;
    const bool phased = d[8 + N];
    const uint32_t bits = d[9 + N];
    const unsigned char* probs = d + 10 + N;
    if (bits < 1 || bits > 32) throw std::runtime_error("BGEN: probabilities bit depth out of bounds");
    const uint64_t avail = (uint64_t)(dlen - 10 - N) * 8;  // bits of probabilities

    // missing samples: the top bit of their ploidy, 8 samples at a time
    for (uint32_t x = 0; x < N; x += 8) {
      if (x + 8 <= N && !(rd64(ploidy + x) & 0x8080808080808080ULL)) continue;
      for (uint32_t y = x; y < std::min(x + 8, N); ++y)
        if (ploidy[y] & 0x80) t.missing.push_back(y);
    }

    if (!phased && min_ploidy == 2 && max_ploidy == 2 && bits == 8) {
      // the common case: diploid, unphased, 8 bits
      if (avail < (uint64_t)N * 16) throw std::runtime_error("BGEN: malformed variant " + std::to_string(vidx));
      for (uint32_t n = 0; n < N; ++n) dose[n] = (float)(2 * (uint32_t)probs[2 * n] + probs[2 * n + 1]) / 255.0f;
    } else {
      // any ploidy and bit depth. Unphased, the stored probabilities of a
      // sample of ploidy p are those of 0..p-1 copies of the second allele, so
      // the first allele has p - k copies in the k-th; phased, one probability
      // of the first allele per haplotype.
      const float maxval = (float)(((uint64_t)1 << bits) - 1);
      const uint64_t mask = ~0ULL >> (64 - bits);
      uint64_t bit = 0;
      for (uint32_t n = 0; n < N; ++n) {
        const uint32_t p = ploidy[n] & 63;
        if (bit + (uint64_t)p * bits > avail)
          throw std::runtime_error("BGEN: malformed variant " + std::to_string(vidx));
        uint64_t sum = 0;
        for (uint32_t k = 0; k < p; ++k) {
          const uint64_t v = (rd64(probs + bit / 8) >> (bit % 8)) & mask;
          bit += bits;
          sum += phased ? v : v * (p - k);
        }
        dose[n] = (float)sum / maxval;
      }
    }
  }

  if (find_minor_allele(dose, N) != 0)
    for (uint32_t n = 0; n < N; ++n) dose[n] = 2.0f - dose[n];
  for (auto n : t.missing) dose[n] = std::nanf("1");
}

void BgenBlockReader::prefetch(const std::vector<uint32_t>& variants) {
  if (!readahead_) readahead_ = std::make_unique<ReadAhead>(path_);
  std::vector<std::pair<uint64_t, uint64_t>> ranges;
  ranges.reserve(variants.size());
  for (auto v : variants) ranges.emplace_back(offsets_[v], offsets_[v + 1] - offsets_[v]);
  readahead_->hint(std::move(ranges));
}

void BgenBlockReader::cancel() {
  if (readahead_) readahead_->cancel();
}

}  // namespace PCAone
