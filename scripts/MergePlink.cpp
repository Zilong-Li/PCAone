#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using uchar = unsigned char;
using uint64 = std::uint64_t;

namespace {

struct InputSet {
  std::string prefix;
  uint64 nsnps = 0;
};

struct Options {
  std::string out = "merged";
  std::string list_file;
  std::vector<std::string> prefixes;
  uint32_t bands = 64;
  bool permute = true;
  bool force = false;
};

void usage(std::ostream& out) {
  out << "Merge PLINK BED files and optionally apply PCAone's SNP permutation.\n\n"
      << "Usage:\n"
      << "  PCAone-merge-plink -o OUT PREFIX [PREFIX ...]\n"
      << "  PCAone-merge-plink -o merged plink.chr{1..22}\n\n"
      << "Options:\n"
      << "  -o, --out OUT       Output PLINK prefix. Required unless using default [merged].\n"
      << "  -l, --list FILE     Text file with one PLINK prefix per line.\n"
      << "  -w, --batches N     Number of PCAone batches/bands for permutation [64].\n"
      << "      --no-permute    Only merge; do not reorder SNP records.\n"
      << "      --force         Overwrite OUT.bed/OUT.bim/OUT.fam if present.\n"
      << "  -h, --help          Print this help message.\n\n"
      << "The input files must be SNP-major PLINK 1 .bed files with identical .fam files.\n"
      << "The default output is ready for: PCAone -b OUT -S ...\n";
}

Options parse_args(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto require_value = [&](const std::string& name) {
      if (i + 1 >= argc) throw std::invalid_argument(name + " requires an argument");
      return std::string(argv[++i]);
    };
    if (arg == "-h" || arg == "--help") {
      usage(std::cout);
      std::exit(EXIT_SUCCESS);
    } else if (arg == "-o" || arg == "--out") {
      opt.out = require_value(arg);
    } else if (arg == "-l" || arg == "--list") {
      opt.list_file = require_value(arg);
    } else if (arg == "-w" || arg == "--batches") {
      auto value = require_value(arg);
      std::size_t pos = 0;
      unsigned long parsed = std::stoul(value, &pos);
      if (pos != value.size() || parsed == 0) throw std::invalid_argument("--batches must be a positive integer");
      opt.bands = static_cast<uint32_t>(parsed);
    } else if (arg == "--no-permute") {
      opt.permute = false;
    } else if (arg == "--force") {
      opt.force = true;
    } else if (!arg.empty() && arg[0] == '-') {
      throw std::invalid_argument("unknown option: " + arg);
    } else {
      opt.prefixes.push_back(arg);
    }
  }
  if (!opt.list_file.empty()) {
    std::ifstream in(opt.list_file);
    if (!in.is_open()) throw std::invalid_argument("cannot open list file: " + opt.list_file);
    std::string line;
    while (std::getline(in, line)) {
      auto first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos || line[first] == '#') continue;
      auto last = line.find_last_not_of(" \t\r\n");
      opt.prefixes.push_back(line.substr(first, last - first + 1));
    }
  }
  if (opt.prefixes.empty()) throw std::invalid_argument("at least one input prefix is required");
  return opt;
}

uint64 count_lines(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open()) throw std::invalid_argument("cannot open: " + path);
  uint64 count = 0;
  std::string line;
  while (std::getline(in, line)) ++count;
  return count;
}

std::string read_text_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) throw std::invalid_argument("cannot open: " + path);
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

void ensure_can_write(const Options& opt) {
  if (opt.force) return;
  for (const auto& ext : {".bed", ".bim", ".fam"}) {
    auto path = opt.out + ext;
    if (std::filesystem::exists(path)) throw std::invalid_argument(path + " exists; use --force to overwrite");
  }
}

std::vector<uint64> band_starts(uint64 nsnps, uint32_t bands) {
  if (bands > nsnps && nsnps > 0) bands = static_cast<uint32_t>(nsnps);
  std::vector<uint64> starts(bands, 0);
  uint64 mod = nsnps % bands;
  uint64 bandsize = (nsnps + bands - 1) / bands;
  for (uint32_t i = 0; i < bands; ++i) {
    if (mod == 0 || i < mod)
      starts[i] = i * bandsize;
    else
      starts[i] = mod * bandsize + (bandsize - 1) * (i - mod);
  }
  return starts;
}

uint64 permuted_index(uint64 original_index, uint64 nsnps, uint32_t bands, const std::vector<uint64>& starts) {
  if (bands > nsnps && nsnps > 0) bands = static_cast<uint32_t>(nsnps);
  return starts[original_index % bands] + original_index / bands;
}

void check_bed_header(std::ifstream& bed, const std::string& path) {
  uchar header[3];
  bed.read(reinterpret_cast<char*>(header), 3);
  if (!bed || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    throw std::invalid_argument(path + " is not a SNP-major PLINK .bed file");
}

std::vector<InputSet> inspect_inputs(const Options& opt, uint64& nsamples, std::string& fam_text) {
  std::vector<InputSet> inputs;
  fam_text = read_text_file(opt.prefixes.front() + ".fam");
  nsamples = count_lines(opt.prefixes.front() + ".fam");
  if (nsamples == 0) throw std::invalid_argument(opt.prefixes.front() + ".fam has no samples");

  for (const auto& prefix : opt.prefixes) {
    auto fam = read_text_file(prefix + ".fam");
    if (fam != fam_text) throw std::invalid_argument(prefix + ".fam differs from the first .fam file");
    uint64 nsnps = count_lines(prefix + ".bim");
    uint64 bytes_per_snp = (nsamples + 3) >> 2;
    auto expected_bed_size = 3 + nsnps * bytes_per_snp;
    auto observed_bed_size = std::filesystem::file_size(prefix + ".bed");
    if (observed_bed_size != expected_bed_size) {
      throw std::invalid_argument(prefix + ".bed size does not match its .bim/.fam dimensions");
    }
    inputs.push_back({prefix, nsnps});
  }
  return inputs;
}

void merge_plink(const Options& opt) {
  ensure_can_write(opt);
  uint64 nsamples = 0;
  std::string fam_text;
  auto inputs = inspect_inputs(opt, nsamples, fam_text);
  uint64 total_snps = 0;
  for (const auto& input : inputs) total_snps += input.nsnps;
  if (total_snps == 0) throw std::invalid_argument("no SNPs found in input .bim files");

  const uint64 bytes_per_snp = (nsamples + 3) >> 2;
  const auto starts = opt.permute ? band_starts(total_snps, opt.bands) : std::vector<uint64>{};
  const uint32_t active_bands =
      opt.permute && opt.bands > total_snps ? static_cast<uint32_t>(total_snps) : opt.bands;

  std::ofstream out_bed(opt.out + ".bed", std::ios::binary);
  std::ofstream out_bim(opt.out + ".bim");
  std::ofstream out_fam(opt.out + ".fam", std::ios::binary);
  if (!out_bed.is_open() || !out_bim.is_open() || !out_fam.is_open())
    throw std::invalid_argument("cannot open output files for prefix: " + opt.out);

  uchar header[3] = {0x6c, 0x1b, 0x01};
  out_bed.write(reinterpret_cast<char*>(header), 3);
  out_bed.seekp(3 + total_snps * bytes_per_snp - 1);
  char zero = 0;
  out_bed.write(&zero, 1);
  out_fam << fam_text;

  std::vector<std::string> bim_lines(total_snps);
  std::vector<char> record(bytes_per_snp);
  uint64 original_index = 0;
  for (const auto& input : inputs) {
    std::ifstream bed(input.prefix + ".bed", std::ios::binary);
    std::ifstream bim(input.prefix + ".bim");
    if (!bed.is_open() || !bim.is_open()) throw std::invalid_argument("cannot open input prefix: " + input.prefix);
    check_bed_header(bed, input.prefix + ".bed");

    std::string bim_line;
    for (uint64 i = 0; i < input.nsnps; ++i, ++original_index) {
      if (!std::getline(bim, bim_line)) throw std::invalid_argument(input.prefix + ".bim ended early");
      bed.read(record.data(), record.size());
      if (!bed) throw std::invalid_argument(input.prefix + ".bed ended early");
      uint64 out_index = opt.permute ? permuted_index(original_index, total_snps, active_bands, starts) : original_index;
      out_bed.seekp(3 + out_index * bytes_per_snp);
      out_bed.write(record.data(), record.size());
      bim_lines[out_index] = bim_line;
    }
  }
  for (const auto& line : bim_lines) out_bim << line << '\n';

  std::cerr << "merged " << inputs.size() << " PLINK prefixes: " << nsamples << " samples, " << total_snps
            << " SNPs";
  if (opt.permute) std::cerr << ", permuted into " << active_bands << " bands";
  std::cerr << "\noutput prefix: " << opt.out << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    auto opt = parse_args(argc, argv);
    merge_plink(opt);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    usage(std::cerr);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
