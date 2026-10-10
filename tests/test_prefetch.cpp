// Prefetch.hpp: the records RecordPrefetcher returns are those of the file,
// whatever the order of the requests (passes, wrap-around, random access, a
// short last block) and with a read in flight; the second buffer is held only
// while reading ahead; ReadAhead never breaks reads.
#include "../src/Prefetch.hpp"

#include <unistd.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(cond)                                                       \
  do {                                                                    \
    if (!(cond)) {                                                        \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond); \
      ++failures;                                                         \
    }                                                                     \
  } while (0)

static unsigned char byte_of(uint64_t record, uint64_t b) { return (unsigned char)((record * 131 + b * 7 + 3) & 0xff); }

int main() {
  char name[] = "/tmp/pcaone-prefetch-XXXXXX";
  const int fd = mkstemp(name);
  CHECK(fd >= 0);
  const uint64_t header = 3, width = 37, n = 1001;  // a short last block
  {
    std::vector<unsigned char> bytes(header + n * width);
    for (uint64_t r = 0; r < n; ++r)
      for (uint64_t b = 0; b < width; ++b) bytes[header + r * width + b] = byte_of(r, b);
    std::ofstream(name, std::ios::binary).write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
  }
  auto same = [&](const unsigned char* p, uint64_t first, uint64_t count) {
    for (uint64_t r = 0; r < count; ++r)
      for (uint64_t b = 0; b < width; ++b)
        if (p[r * width + b] != byte_of(first + r, b)) return false;
    return true;
  };

  {
    // min_share 0: always read ahead, whatever the speed of the reads
    PCAone::RecordPrefetcher reader(name, header, width, n, 0.0);
    const uint64_t block = 64;
    for (int pass = 0; pass < 3; ++pass)
      for (uint64_t first = 0; first < n; first += block) {
        const uint64_t count = std::min(block, n - first);
        CHECK(same(reader.get(first, count), first, count));
      }
    // every block after the very first came from the background
    CHECK(reader.predicted() == 3 * ((n + block - 1) / block) - 1);
    // out of order, as LD reads blocks: still the right bytes
    std::mt19937 rng(7);
    for (int i = 0; i < 200; ++i) {
      const uint64_t first = rng() % n, count = 1 + rng() % std::min<uint64_t>(100, n - first);
      CHECK(same(reader.get(first, count), first, count));
    }
    // back to a pass after random access
    for (uint64_t first = 0; first < n; first += 100) {
      const uint64_t count = std::min<uint64_t>(100, n - first);
      CHECK(same(reader.get(first, count), first, count));
    }
    bool threw = false;
    try {
      reader.get(n - 1, 2);
    } catch (const std::out_of_range&) {
      threw = true;
    }
    CHECK(threw);
  }  // destroyed with the first block of the next pass in flight

  {
    // reads from the page cache and slow work between the calls: read on the spot
    PCAone::RecordPrefetcher reader(name, header, width, n);
    for (int pass = 0; pass < 2; ++pass)
      for (uint64_t first = 0; first < n; first += 100) {
        const uint64_t count = std::min<uint64_t>(100, n - first);
        CHECK(same(reader.get(first, count), first, count));
        usleep(2000);
      }
    CHECK(reader.predicted() <= 1);  // only the first block, before any timing
    // reading on the spot: the second buffer is given back
    CHECK(reader.buffer_bytes() <= 100 * width);
  }

  {
    // reading ahead holds two blocks
    PCAone::RecordPrefetcher reader(name, header, width, n, 0.0);
    for (uint64_t first = 0; first < 300; first += 100) CHECK(same(reader.get(first, 100), first, 100));
    reader.cancel();
    CHECK(reader.buffer_bytes() >= 2 * 100 * width);
  }

  {
    // a block as large as the file: nothing to read ahead
    PCAone::RecordPrefetcher reader(name, header, width, n);
    for (int pass = 0; pass < 2; ++pass) CHECK(same(reader.get(0, n), 0, n));
    CHECK(reader.predicted() == 0);
  }

  {
    // page cache requests: never an error, and the file reads as before
    PCAone::ReadAhead ahead(name);
    ahead.hint({{header, 100}, {header + 5000, 10}, {header + 90, 50}, {0, 0}});
    ahead.wait();
    ahead.hint(header, n * width);
    ahead.cancel();
    ahead.hint(n * width * 10, 4096);  // beyond the end: only a hint
    ahead.wait();
    ahead.hint({{0, 0}});  // nothing to ask for
    ahead.wait();
    CHECK(ahead.hints() == 3);
    PCAone::RecordPrefetcher reader(name, header, width, n);
    CHECK(same(reader.get(0, n), 0, n));
  }
  {
    // set_order: record r is record order[r] of the file (the BED permutation
    // read from the input), in passes with reads ahead, at random, and with
    // the reads ahead off; runs of consecutive records and single ones
    std::vector<uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::mt19937 rng(11);
    std::shuffle(order.begin(), order.end(), rng);
    for (uint64_t b = 0; b < n; b += 250) std::sort(order.begin() + b, order.begin() + std::min(n, b + 250));
    std::iota(order.begin() + 500, order.begin() + 600, 500);  // a run of 100
    auto same_order = [&](const unsigned char* p, uint64_t first, uint64_t count) {
      for (uint64_t r = 0; r < count; ++r)
        for (uint64_t b = 0; b < width; ++b)
          if (p[r * width + b] != byte_of(order[first + r], b)) return false;
      return true;
    };
    PCAone::RecordPrefetcher reader(name, header, width, n, 0.0);
    reader.set_order(order);
    for (int pass = 0; pass < 2; ++pass)
      for (uint64_t first = 0; first < n; first += 50) {
        const uint64_t count = std::min<uint64_t>(50, n - first);
        CHECK(same_order(reader.get(first, count), first, count));
      }
    CHECK(reader.predicted() > 0);
    for (int i = 0; i < 100; ++i) {
      const uint64_t first = rng() % n, count = 1 + rng() % std::min<uint64_t>(100, n - first);
      CHECK(same_order(reader.get(first, count), first, count));
    }
    PCAone::RecordPrefetcher quiet(name, header, width, n, 0.0);
    quiet.set_order(order);
    quiet.set_read_ahead(false);
    for (uint64_t first = 0; first < n; first += 100)
      CHECK(same_order(quiet.get(first, std::min<uint64_t>(100, n - first)), first, std::min<uint64_t>(100, n - first)));
    CHECK(quiet.predicted() == 0);
    bool threw = false;
    try {
      std::vector<uint32_t> bad(order);
      bad[3] = n;
      quiet.set_order(bad);
    } catch (const std::out_of_range&) {
      threw = true;
    }
    CHECK(threw);
  }

  {
    // set_order with records larger than a page
    char big[] = "/tmp/pcaone-prefetch-big-XXXXXX";
    const int bfd = mkstemp(big);
    CHECK(bfd >= 0);
    const uint64_t bw = 5000, bn = 120;
    std::vector<unsigned char> bytes(header + bn * bw);
    for (uint64_t r = 0; r < bn; ++r)
      for (uint64_t b = 0; b < bw; ++b) bytes[header + r * bw + b] = byte_of(r, b);
    std::ofstream(big, std::ios::binary).write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    std::vector<uint32_t> order(bn);
    for (uint64_t r = 0; r < bn; ++r) order[r] = (uint32_t)((r * 7) % bn);  // 7 and 120 are coprime
    PCAone::RecordPrefetcher reader(big, header, bw, bn, 0.0);
    reader.set_order(order);
    for (int pass = 0; pass < 2; ++pass)
      for (uint64_t first = 0; first < bn; first += 16) {
        const uint64_t count = std::min<uint64_t>(16, bn - first);
        const unsigned char* p = reader.get(first, count);
        bool ok = true;
        for (uint64_t r = 0; r < count && ok; ++r)
          for (uint64_t b = 0; b < bw; ++b)
            if (p[r * bw + b] != byte_of(order[first + r], b)) {
              ok = false;
              break;
            }
        CHECK(ok);
      }
    ::close(bfd);
    std::remove(big);
  }

  bool threw = false;
  try {
    PCAone::ReadAhead missing(std::string(name) + ".none");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);

  ::close(fd);
  std::remove(name);
  if (failures) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("Prefetch: records match the file in passes, at random, with a read in flight and through set_order\n");
  return 0;
}
