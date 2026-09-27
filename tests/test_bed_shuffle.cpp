#include "../src/BedShuffle.hpp"
#include <cassert>
#include <numeric>
#include <random>
#include <sstream>
#include <string>

int main() {
  // Odd record widths, a single chunk and many chunks, more buckets than
  // records per chunk, and a partial final bucket must all give the bytes of
  // the full shuffle regrouped stably within each bucket.
  const std::string header("\x6c\x1b\x01", 3);
  for (uint32_t n : {1u, 3u, 17u, 101u}) {
    for (uint64_t width : {1u, 3u, 19u}) {
      std::string data(n * width, '\0');
      for (size_t j = 0; j < data.size(); ++j) data[j] = (j * 73 + j / 7) % 256;
      std::vector<uint32_t> order(n);
      std::iota(order.begin(), order.end(), 0);
      std::mt19937 rng(42);  // any fixed order will do here
      std::shuffle(order.begin(), order.end(), rng);
      for (uint64_t bucket : {1u, 3u, 11u, 200u}) {
        auto stable = order;
        for (size_t first = 0; first < n; first += bucket)
          std::sort(stable.begin() + first, stable.begin() + std::min<uint64_t>(n, first + bucket));
        std::string expected = header;
        for (auto source : stable) expected += data.substr(source * width, width);
        for (uint64_t capacity : {1u, 2u, 7u, 200u}) {
          auto actual = order;
          std::istringstream in(data);
          // pre-sized: a stringstream, unlike a file, cannot seek past its end
          std::stringstream out(std::string(3 + data.size(), '\0'));
          out.write(header.data(), 3);
          PCAone::rewrite_bed_buckets(in, out, actual, width, 2 * capacity * width, bucket);
          assert(actual == stable);
          assert(out.str() == expected);
          assert(in.tellg() == static_cast<std::streamoff>(data.size()));
        }
      }
      std::istringstream truncated(data.substr(1));
      std::stringstream out(std::string(3 + data.size(), '\0'));
      auto actual = order;
      bool failed = false;
      try { PCAone::rewrite_bed_buckets(truncated, out, actual, width, 2 * width, 3); }
      catch (const std::runtime_error&) { failed = true; }
      assert(failed);
      failed = false;
      try { PCAone::rewrite_bed_buckets(truncated, out, actual, width, width, 3); }  // buffer < two records
      catch (const std::invalid_argument&) { failed = true; }
      assert(failed);
    }
  }
}
