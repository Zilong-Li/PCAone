// Prefetch.hpp: the records RecordPrefetcher returns are those of the file,
// whatever the order of the requests (passes, wrap-around, random access, a
// short last block) and with a read in flight; ReadAhead never breaks reads.
#include "../src/Prefetch.hpp"

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
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
    PCAone::RecordPrefetcher reader(name, header, width, n);
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
    // a block as large as the file: nothing to read ahead
    PCAone::RecordPrefetcher reader(name, header, width, n);
    for (int pass = 0; pass < 2; ++pass) CHECK(same(reader.get(0, n), 0, n));
    CHECK(reader.predicted() == 0);
  }

  {
    PCAone::ReadAhead ahead(name, 3);
    ahead.hint({{header, 100}, {header + 5000, 10}, {header + 90, 50}, {0, 0}});
    ahead.wait();
    ahead.hint(header, n * width);
    ahead.cancel();
    ahead.hint(n * width * 10, 4096);  // beyond the end: only a hint
    ahead.wait();
  }

  ::close(fd);
  std::remove(name);
  if (failures) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("Prefetch: records match the file in passes, at random and with a read in flight\n");
  return 0;
}
