#ifndef PCAONE_BED_SHUFFLE_HPP
#define PCAONE_BED_SHUFFLE_HPP

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace PCAone {
// Rewrite the packed BED payload so that bucket b, the records
// [b * bucket, (b + 1) * bucket), holds the same SNPs as in the full shuffle
// `order` (order[destination] = source), kept in source order. On return,
// order holds the order actually written.
//
// Out-of-core winSVD (Halko.cpp) changes Omg only at multiples of bandFactor
// read blocks or at the end of a pass. With bucket = blocksize * bandFactor,
// the order within a bucket therefore does not change its updates in exact
// arithmetic, and the PCs equal those of the full shuffle.
//
// One sequential read of the input. The input is processed in chunks of
// budget / 2 bytes (the other half holds the regrouped chunk), and each chunk
// costs at most one seek and write per bucket. With about -w buckets the
// writes stay large however big the file is.
inline void rewrite_bed_buckets(std::istream& in, std::ostream& out, std::vector<uint32_t>& order,
                                uint64_t width, uint64_t budget, uint64_t bucket) {
  if (!width || budget / 2 < width || order.empty() || !bucket)
    throw std::invalid_argument("Invalid BED shuffle dimensions or buffer");
  const uint64_t n = order.size();
  const uint64_t capacity = std::min(n, budget / 2 / width);
  const uint64_t buckets = (n + bucket - 1) / bucket;
  const auto payload_start = out.tellp();  // preserve the caller's BED header
  if (payload_start == std::streampos(-1)) throw std::runtime_error("BED output must be seekable");
  std::vector<uint32_t> membership(n);
  for (uint64_t d = 0; d < n; ++d) membership[order[d]] = d / bucket;
  std::vector<char> input(capacity * width), output(capacity * width);
  std::vector<uint64_t> written(buckets, 0), counts(buckets), starts(buckets), next(buckets);
  for (uint64_t first = 0; first < n; first += capacity) {
    const uint64_t count = std::min(capacity, n - first);
    in.read(input.data(), count * width);
    if (!in) throw std::runtime_error("Short BED read during permutation");
    std::fill(counts.begin(), counts.end(), 0);
    for (uint64_t j = 0; j < count; ++j) ++counts[membership[first + j]];
    uint64_t offset = 0;
    for (uint64_t b = 0; b < buckets; ++b) {
      starts[b] = next[b] = offset;
      offset += counts[b];
    }
    for (uint64_t j = 0; j < count; ++j) {
      const uint64_t b = membership[first + j];
      const uint64_t slot = next[b]++;
      std::memcpy(output.data() + slot * width, input.data() + j * width, width);
      order[b * bucket + written[b] + slot - starts[b]] = first + j;
    }
    for (uint64_t b = 0; b < buckets; ++b) {
      if (!counts[b]) continue;
      out.seekp(payload_start + static_cast<std::streamoff>((b * bucket + written[b]) * width));
      out.write(output.data() + starts[b] * width, counts[b] * width);
      if (!out) throw std::runtime_error("Cannot write permuted BED file");
      written[b] += counts[b];
    }
  }
}
}  // namespace PCAone
#endif
