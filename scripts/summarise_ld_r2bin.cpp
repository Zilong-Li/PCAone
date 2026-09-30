// -*- compile-command: "g++ summarise_ld_r2bin.cpp -o summarise_ld_r2bin -O3 -std=c++17 -lz"; -*-
//
// Bin the SNP pairs of an LD file by distance, per chromosome, for
// scripts/plot-ld-decay.R. It streams the file, so it handles .ld.gz files too
// large to read into R, and writes the same bins and table as
// `plot-ld-decay.R --save-bins`:
//
//   summarise_ld_r2bin -i adj.ld.gz --max 1000000 -o adj.decay.tsv
//   Rscript plot-ld-decay.R adj.decay.tsv std.decay.tsv
//
// Input: PCAone -R/--print-r2 output, or plink --ld (use --ld-window-r2 0),
// gzipped or not. Columns are found by name: CHR_A, BP_A, BP_B, R2 and, if
// present, CHR_B; pairs on different chromosomes go to one row "inter".
//
// Bins: [0, --min), then --bins log-spaced bins up to --max; with --linear,
// --bins equal bins from 0 to --max. Set --max to the --ld-bp of the run.
//
// Output (tsv): chr lo hi n mean_dist mean_r2 sd_r2

#include <zlib.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct Acc {
  long long n = 0;
  long double dist = 0, r2 = 0, r2sq = 0;
  void add(double d, double r) {
    n++;
    dist += d;
    r2 += r;
    r2sq += (long double)r * r;
  }
};

static void usage() {
  std::cerr << "Bin LD pairs by distance for scripts/plot-ld-decay.R\n\n"
               "Usage: summarise_ld_r2bin -i FILE [options]\n\n"
               "  -i FILE      .ld.gz/.ld from PCAone -R or plink --ld (- for stdin)\n"
               "  -o FILE      output table [stdout]\n"
               "  --bins N     log-spaced bins between --min and --max [40]\n"
               "  --min D      lower edge of the first log bin, bp [1000]\n"
               "  --max D      upper edge of the last bin, bp; use the run's --ld-bp [1000000]\n"
               "  --linear     N equal bins from 0 to --max instead\n";
}

// the same edges as ld_edges() in plot-ld-decay.R (R's seq(from, to, length.out))
static std::vector<double> make_edges(int bins, double min, double max, bool linear) {
  std::vector<double> e;
  if (linear) {
    e.push_back(0);
    for (int i = 1; i < bins; i++) e.push_back(i * (max / bins));
    e.push_back(max);
    return e;
  }
  const double a = std::log10(min), b = std::log10(max), by = (b - a) / bins;
  e.push_back(0);
  e.push_back(std::pow(10.0, a));
  for (int i = 1; i < bins; i++) e.push_back(std::pow(10.0, a + i * by));
  e.push_back(std::pow(10.0, b));
  return e;
}

// split on blanks and tabs, in place; returns the number of fields
static int split(char* s, std::vector<char*>& f) {
  f.clear();
  char* p = s;
  while (*p) {
    while (*p == ' ' || *p == '\t') p++;
    if (!*p || *p == '\n' || *p == '\r') break;
    f.push_back(p);
    while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
    if (*p) *p++ = '\0';
  }
  return (int)f.size();
}

static std::string strip_chr(const char* c) {
  if ((c[0] == 'c' || c[0] == 'C') && (c[1] == 'h' || c[1] == 'H') && (c[2] == 'r' || c[2] == 'R')) return c + 3;
  return c;
}

int main(int argc, char** argv) {
  std::string in, out;
  int bins = 40;
  double min = 1000, max = 1000000;
  bool linear = false;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto val = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "error: " << a << " needs a value\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "-h" || a == "--help") {
      usage();
      return 0;
    } else if (a == "-i") in = val();
    else if (a == "-o") out = val();
    else if (a == "--bins") bins = std::atoi(val().c_str());
    else if (a == "--min") min = std::atof(val().c_str());
    else if (a == "--max") max = std::atof(val().c_str());
    else if (a == "--linear") linear = true;
    else {
      std::cerr << "error: unknown option " << a << "\n\n";
      usage();
      return 1;
    }
  }
  if (in.empty()) {
    usage();
    return 1;
  }
  if (bins < 1 || max <= 0 || (!linear && !(min > 0 && min < max))) {
    std::cerr << "error: need --bins >= 1 and 0 < --min < --max\n";
    return 1;
  }
  const std::vector<double> e = make_edges(bins, min, max, linear);
  const size_t nbin = e.size() - 1;

  gzFile fp = in == "-" ? gzdopen(0, "rb") : gzopen(in.c_str(), "rb");
  if (!fp) {
    std::cerr << "error: cannot open " << in << "\n";
    return 1;
  }
  gzbuffer(fp, 1 << 20);

  std::vector<char> buf(1 << 16);
  std::string line;
  std::vector<char*> f;
  auto getline = [&]() -> bool {  // whole line, however long
    line.clear();
    while (gzgets(fp, buf.data(), (int)buf.size())) {
      line += buf.data();
      if (!line.empty() && line.back() == '\n') return true;
    }
    return !line.empty();
  };

  if (!getline()) {
    std::cerr << "error: " << in << " is empty\n";
    return 1;
  }
  int cA = -1, pA = -1, cB = -1, pB = -1, cR = -1;
  int nf = split(&line[0], f);
  for (int k = 0; k < nf; k++) {
    if (!std::strcmp(f[k], "CHR_A")) cA = k;
    else if (!std::strcmp(f[k], "BP_A")) pA = k;
    else if (!std::strcmp(f[k], "CHR_B")) cB = k;
    else if (!std::strcmp(f[k], "BP_B")) pB = k;
    else if (!std::strcmp(f[k], "R2")) cR = k;
  }
  if (cA < 0 || pA < 0 || pB < 0 || cR < 0) {
    std::cerr << "error: " << in << " needs the columns CHR_A, BP_A, BP_B and R2\n";
    return 1;
  }
  const int need = std::max({cA, pA, cB, pB, cR}) + 1;

  std::vector<std::string> chrs;  // in order of appearance
  std::unordered_map<std::string, size_t> idx;
  std::vector<std::vector<Acc>> acc;
  Acc inter;
  long long far = 0, bad = 0, lines = 0;
  std::string prev;
  size_t cur = 0;
  while (getline()) {
    lines++;
    if (split(&line[0], f) < need) {
      bad++;
      continue;
    }
    const char* ca = f[cA];
    const double r2 = std::strtod(f[cR], nullptr);
    const double d = std::fabs(std::strtod(f[pB], nullptr) - std::strtod(f[pA], nullptr));
    if (cB >= 0 && std::strcmp(ca, f[cB])) {
      inter.add(d, r2);
      continue;
    }
    if (prev != ca) {  // chromosome change: look it up once
      prev = ca;
      std::string c = strip_chr(ca);
      auto it = idx.find(c);
      if (it == idx.end()) {
        it = idx.emplace(c, chrs.size()).first;
        chrs.push_back(c);
        acc.emplace_back(nbin);
      }
      cur = it->second;
    }
    size_t b = std::upper_bound(e.begin(), e.end(), d) - e.begin();  // edges <= d
    if (b > nbin) {
      if (d > e.back()) {
        far++;
        continue;
      }
      b = nbin;  // d == max: last bin is closed
    }
    acc[cur][b - 1].add(d, r2);
  }
  gzclose(fp);
  if (bad) std::cerr << "warning: " << bad << " lines with too few columns skipped\n";
  if (far) std::cerr << "warning: " << far << " pairs beyond --max " << max << " ignored; raise --max\n";

  FILE* o = out.empty() ? stdout : std::fopen(out.c_str(), "w");
  if (!o) {
    std::cerr << "error: cannot write " << out << "\n";
    return 1;
  }
  auto row = [&](const std::string& c, const char* lo, const char* hi, const Acc& a) {
    const long double n = a.n;
    const double sd = std::sqrt((double)std::max<long double>(0, (a.r2sq - a.r2 * a.r2 / n) / std::max<long double>(n - 1, 1)));
    std::fprintf(o, "%s\t%s\t%s\t%lld\t%.7g\t%.7g\t%.7g\n", c.c_str(), lo, hi, a.n, (double)(a.dist / n),
                 (double)(a.r2 / n), sd);
  };
  std::fprintf(o, "chr\tlo\thi\tn\tmean_dist\tmean_r2\tsd_r2\n");
  char lo[32], hi[32];
  long long total = inter.n;
  for (size_t c = 0; c < chrs.size(); c++) {
    for (size_t b = 0; b < nbin; b++) {
      if (!acc[c][b].n) continue;
      std::snprintf(lo, sizeof lo, "%.15g", e[b]);
      std::snprintf(hi, sizeof hi, "%.15g", e[b + 1]);
      row(chrs[c], lo, hi, acc[c][b]);
      total += acc[c][b].n;
    }
  }
  if (inter.n) row("inter", "NA", "NA", inter);
  if (o != stdout) std::fclose(o);
  std::cerr << "binned " << total << " of " << lines << " pairs\n";
  return 0;
}
