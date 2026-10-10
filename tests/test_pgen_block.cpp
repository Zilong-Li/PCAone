// PgenBlockReader (src/PgenBlock.cpp) on a PLINK 1 .bed, which pgenlib reads
// as a PGEN: the calls are those of the file, and a constructor that fails
// (wrong sample count, missing file, too few file descriptors for the readers)
// gives back every descriptor and mapping it took, because FilePgen goes on
// with PgenReader after catching the exception.
#include "../src/PgenBlock.hpp"

#include <dirent.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(cond, ...)                                              \
  do {                                                                \
    if (!(cond)) {                                                    \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed: ", __FILE__, __LINE__, #cond); \
      std::fprintf(stderr, __VA_ARGS__);                              \
      std::fprintf(stderr, "\n");                                     \
      ++failures;                                                     \
    }                                                                 \
  } while (0)

// open descriptors, and the highest one (the listing's own is counted each time)
static int open_fds(int* highest = nullptr) {
  DIR* d = opendir("/dev/fd");
  if (!d) return -1;
  int n = 0, hi = -1;
  while (dirent* e = readdir(d)) {
    if (e->d_name[0] == '.') continue;
    ++n;
    hi = std::max(hi, std::atoi(e->d_name));
  }
  closedir(d);
  if (highest) *highest = hi;
  return n;
}

static std::string bed;  // the test file

// mappings of the test file (Linux), or -1
static int mappings() {
  std::ifstream maps("/proc/self/maps");
  if (!maps) return -1;
  int n = 0;
  for (std::string line; std::getline(maps, line);) n += line.find(bed) != std::string::npos;
  return n;
}

// constructing a reader throws, and leaves the process as it found it
template <typename F>
static void throws_cleanly(const char* what, F make) {
  const int fds = open_fds(), maps = mappings();
  bool threw = false;
  try {
    make();
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw, "%s: no exception", what);
  const int fds_after = open_fds(), maps_after = mappings();
  CHECK(fds_after == fds, "%s: %d open descriptors before, %d after", what, fds, fds_after);
  if (maps >= 0) CHECK(maps_after == maps, "%s: %d mappings of the file before, %d after", what, maps, maps_after);
}

int main() {
  const uint32_t N = 301, M = 257;  // N not a multiple of 4: a partial last byte
  const uint32_t bytes = (N + 3) / 4;
  const char* tmp = std::getenv("TMPDIR");
  std::string name = std::string(tmp && *tmp ? tmp : "/tmp") + "/pcaone-pgenblock-XXXXXX";
  std::vector<char> path(name.begin(), name.end());
  path.push_back('\0');
  const int fd = mkstemp(path.data());
  CHECK(fd >= 0, "mkstemp");
  ::close(fd);
  bed = path.data();

  // PLINK 1 codes: 00 hom A1 (2 ALT), 01 missing, 10 het, 11 hom A2 (0 ALT)
  std::mt19937 rng(11);
  std::vector<std::vector<int>> alt(M, std::vector<int>(N));
  {
    std::vector<unsigned char> file = {0x6c, 0x1b, 0x01};
    for (uint32_t v = 0; v < M; ++v) {
      std::vector<unsigned char> rec(bytes, 0);
      for (uint32_t s = 0; s < N; ++s) {
        const unsigned code = (v % 7 == 0) ? 3u : rng() % 4;  // some monomorphic variants
        alt[v][s] = code == 0 ? 2 : code == 1 ? -3 : code == 2 ? 1 : 0;
        rec[s / 4] |= (unsigned char)(code << (2 * (s % 4)));
      }
      file.insert(file.end(), rec.begin(), rec.end());
    }
    std::ofstream(bed, std::ios::binary).write(reinterpret_cast<const char*>(file.data()), (std::streamsize)file.size());
  }

  const int fds0 = open_fds();
  std::string probe;
  {
    PCAone::PgenBlockReader reader(bed, N, 3);
    probe = reader.probe_name();
    CHECK(reader.variant_ct() == M, "variant_ct %u", reader.variant_ct());
    std::vector<uint32_t> all(M);
    for (uint32_t v = 0; v < M; ++v) all[v] = M - 1 - v;  // any order
    reader.load(all);
    reader.prefetch(all);  // the next "block", then stopped
    reader.cancel();
    std::vector<double> buf(N);
    for (uint32_t v = 0; v < M; ++v) {
      const int thr = (int)(v % 3);
      reader.hardcalls(thr, v, buf.data());
      int bad = 0;
      for (uint32_t s = 0; s < N; ++s) bad += buf[s] != alt[v][s];
      CHECK(!bad, "hard calls of variant %u: %d samples differ", v, bad);
      reader.dosages(thr, v, buf.data());
      bad = 0;
      for (uint32_t s = 0; s < N; ++s) bad += buf[s] != alt[v][s];
      CHECK(!bad, "dosages of variant %u: %d samples differ", v, bad);
      const uintptr_t* g = reader.genovec(thr, v);
      bad = 0;
      for (uint32_t s = 0; s < N; ++s) {
        const int code = (int)((g[s / (4 * sizeof(uintptr_t))] >> (2 * (s % (4 * sizeof(uintptr_t))))) & 3);
        bad += code != (alt[v][s] < 0 ? 3 : alt[v][s]);
      }
      CHECK(!bad, "genovec of variant %u: %d samples differ", v, bad);
    }
  }
  CHECK(open_fds() == fds0, "a reader that was used leaves %d descriptors, not %d", open_fds(), fds0);

  // failures in pgenlib's header checks: its file is already open then
  throws_cleanly("wrong sample count", [&] { PCAone::PgenBlockReader r(bed, N + 4, 2); });
  throws_cleanly("missing file", [&] { PCAone::PgenBlockReader r(bed + ".none", N, 2); });

  // failure while opening the readers' files, after the descriptor, the page
  // cache probe and the first readers were set up
  {
    int highest = -1;
    open_fds(&highest);
    struct rlimit old;
    CHECK(getrlimit(RLIMIT_NOFILE, &old) == 0, "getrlimit");
    const int fds = open_fds(), maps = mappings();
    // RLIMIT_NOFILE bounds the descriptor numbers, and a new descriptor takes
    // the lowest free one. Fill the free numbers below the highest open
    // descriptor (a CI runner can leave e.g. 145 open), so that only the few
    // above it are left for the readers.
    std::vector<int> fillers;
    for (int f; (f = ::open("/dev/null", O_RDONLY)) >= 0;) {
      if (f > highest) {
        ::close(f);
        break;
      }
      fillers.push_back(f);
    }
    struct rlimit low = old;
    low.rlim_cur = (rlim_t)highest + 1 + 4;  // room for a few files, not for 32 readers
    bool threw = false;
    if (setrlimit(RLIMIT_NOFILE, &low) == 0) {
      try {
        PCAone::PgenBlockReader r(bed, N, 32);
      } catch (const std::exception&) {
        threw = true;
      }
      setrlimit(RLIMIT_NOFILE, &old);
      for (int f : fillers) ::close(f);
      fillers.clear();
      CHECK(threw, "32 readers under a limit of %ld descriptors: no exception", (long)low.rlim_cur);
      const int fds_after = open_fds(), maps_after = mappings();
      CHECK(fds_after == fds, "descriptor limit: %d open descriptors before, %d after", fds, fds_after);
      if (maps >= 0)
        CHECK(maps_after == maps, "descriptor limit: %d mappings of the file before, %d after", maps, maps_after);
    } else {
      for (int f : fillers) ::close(f);
      std::printf("skip the descriptor limit case: setrlimit failed\n");
    }
  }

  // and a reader can still be made afterwards
  {
    PCAone::PgenBlockReader reader(bed, N, 4);
    std::vector<double> buf(N);
    reader.hardcalls(3, 1, buf.data());
    CHECK(buf[0] == alt[1][0], "after the failures");
  }
  CHECK(open_fds() == fds0, "at the end %d descriptors, not %d", open_fds(), fds0);

  std::remove(bed.c_str());
  if (failures) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("PgenBlockReader: calls of a .bed match; failed constructions leave no descriptor or mapping "
              "open (page cache probe: %s)\n",
              probe.c_str());
  return 0;
}
