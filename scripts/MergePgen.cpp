#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using uchar = unsigned char;
using uint64 = std::uint64_t;

namespace {

struct Options {
  std::string out = "merged";
  std::string list_file;
  std::string plink2 = "plink2";
  std::vector<std::string> prefixes;
  uint32_t bands = 64;
  bool permute = true;
  bool force = false;
  bool keep_temp = false;
};

void usage(std::ostream& out) {
  out << "Merge PLINK2 PGEN files and optionally apply PCAone's SNP permutation.\n\n"
      << "Usage:\n"
      << "  PCAone-merge-pgen -o OUT PREFIX [PREFIX ...]\n"
      << "  PCAone-merge-pgen -o merged pgen.chr{1..22}\n\n"
      << "Options:\n"
      << "  -o, --out OUT       Output PLINK2 prefix. Required unless using default [merged].\n"
      << "  -l, --list FILE     Text file with one PLINK2 prefix per line.\n"
      << "  -w, --batches N     Number of PCAone batches/bands for permutation [64].\n"
      << "      --plink2 PATH   Path to plink2 executable [plink2].\n"
      << "      --no-permute    Only merge; do not reorder variant records.\n"
      << "      --force         Overwrite output/temp files if present.\n"
      << "      --keep-temp     Keep intermediate merged files.\n"
      << "  -h, --help          Print this help message.\n\n"
      << "This wrapper uses plink2 to write fixed-width hardcall PGEN format=2 files.\n"
      << "Dosages and multiallelic variants are not preserved by this helper.\n"
      << "The default output is ready for: PCAone -p OUT -S ...\n";
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
    } else if (arg == "--plink2") {
      opt.plink2 = require_value(arg);
    } else if (arg == "--no-permute") {
      opt.permute = false;
    } else if (arg == "--force") {
      opt.force = true;
    } else if (arg == "--keep-temp") {
      opt.keep_temp = true;
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

std::string shell_quote(const std::string& s) {
  std::string out = "'";
  for (char c : s) {
    if (c == '\'')
      out += "'\\''";
    else
      out += c;
  }
  out += "'";
  return out;
}

void run_command(const std::string& cmd) {
  std::cerr << cmd << "\n";
  int rc = std::system(cmd.c_str());
  if (rc != 0) throw std::runtime_error("command failed with status " + std::to_string(rc));
}

void remove_fileset(const std::string& prefix) {
  for (const auto& ext : {".pgen", ".pvar", ".psam", ".log", ".pgen.pgi"}) std::filesystem::remove(prefix + ext);
}

void remove_pmerge_sidecars(const std::string& prefix) {
  for (const auto& ext : {".pgen", ".pvar", ".psam"}) std::filesystem::remove(prefix + "-merge" + ext);
}

void ensure_can_write(const Options& opt, const std::string& prefix) {
  if (opt.force) {
    remove_fileset(prefix);
    remove_pmerge_sidecars(prefix);
    return;
  }
  for (const auto& ext : {".pgen", ".pvar", ".psam"}) {
    auto path = prefix + ext;
    if (std::filesystem::exists(path)) throw std::invalid_argument(path + " exists; use --force to overwrite");
  }
}

uint64 count_psam_samples(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open()) throw std::invalid_argument("cannot open: " + path);
  uint64 count = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] != '#') ++count;
  }
  return count;
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

void write_pmerge_list(const Options& opt, const std::string& path) {
  std::ofstream out(path);
  if (!out.is_open()) throw std::invalid_argument("cannot write: " + path);
  for (const auto& prefix : opt.prefixes) out << prefix << '\n';
}

void merge_with_plink2(const Options& opt, const std::string& out_prefix) {
  std::string list_path = out_prefix + ".pmerge-list";
  write_pmerge_list(opt, list_path);
  std::string cmd = shell_quote(opt.plink2) + " --pmerge-list " + shell_quote(list_path) +
                    " pfile --indiv-sort none --make-pgen format=2 --out " + shell_quote(out_prefix);
  run_command(cmd);
  if (!opt.keep_temp) std::filesystem::remove(list_path);
}

void permute_fixed_pgen(const Options& opt, const std::string& in_prefix, const std::string& out_prefix) {
  ensure_can_write(opt, out_prefix);
  uint64 nsamples = count_psam_samples(in_prefix + ".psam");
  if (nsamples == 0) throw std::invalid_argument(in_prefix + ".psam has no samples");
  uint64 bytes_per_variant = (nsamples + 3) >> 2;

  std::ifstream in_pvar(in_prefix + ".pvar");
  if (!in_pvar.is_open()) throw std::invalid_argument("cannot open: " + in_prefix + ".pvar");
  std::vector<std::string> header_lines;
  std::vector<std::string> variant_lines;
  std::string line;
  while (std::getline(in_pvar, line)) {
    if (!line.empty() && line[0] == '#')
      header_lines.push_back(line);
    else
      variant_lines.push_back(line);
  }
  uint64 nsnps = variant_lines.size();
  if (nsnps == 0) throw std::invalid_argument(in_prefix + ".pvar has no variants");

  std::ifstream in_pgen(in_prefix + ".pgen", std::ios::binary);
  std::ofstream out_pgen(out_prefix + ".pgen", std::ios::binary);
  std::ofstream out_pvar(out_prefix + ".pvar");
  if (!in_pgen.is_open() || !out_pgen.is_open() || !out_pvar.is_open())
    throw std::invalid_argument("cannot open PGEN/PVAR files for permutation");

  uchar header[12];
  in_pgen.read(reinterpret_cast<char*>(header), sizeof(header));
  if (!in_pgen) throw std::invalid_argument(in_prefix + ".pgen is too small for fixed-width PGEN format=2");
  out_pgen.write(reinterpret_cast<char*>(header), sizeof(header));

  auto expected_size = sizeof(header) + nsnps * bytes_per_variant;
  auto observed_size = std::filesystem::file_size(in_prefix + ".pgen");
  if (observed_size != expected_size)
    throw std::invalid_argument(in_prefix + ".pgen is not fixed-width PGEN format=2 with expected dimensions");

  const auto starts = band_starts(nsnps, opt.bands);
  uint32_t active_bands = opt.bands > nsnps ? static_cast<uint32_t>(nsnps) : opt.bands;
  std::vector<std::string> out_variants(nsnps);
  std::vector<char> record(bytes_per_variant);
  for (uint64 i = 0; i < nsnps; ++i) {
    in_pgen.read(record.data(), record.size());
    if (!in_pgen) throw std::invalid_argument(in_prefix + ".pgen ended early");
    uint64 out_index = permuted_index(i, nsnps, active_bands, starts);
    out_pgen.seekp(sizeof(header) + out_index * bytes_per_variant);
    out_pgen.write(record.data(), record.size());
    out_variants[out_index] = variant_lines[i];
  }

  for (const auto& h : header_lines) out_pvar << h << '\n';
  for (const auto& v : out_variants) out_pvar << v << '\n';
  std::filesystem::copy_file(in_prefix + ".psam", out_prefix + ".psam", std::filesystem::copy_options::overwrite_existing);
  std::cerr << "permuted " << nsnps << " variants into " << active_bands << " bands\n";
}

void merge_pgen(const Options& opt) {
  if (!opt.permute) {
    ensure_can_write(opt, opt.out);
    merge_with_plink2(opt, opt.out);
    std::cerr << "output prefix: " << opt.out << "\n";
    return;
  }

  std::string tmp_prefix = opt.out + ".merge-pgen.tmp";
  ensure_can_write(opt, tmp_prefix);
  merge_with_plink2(opt, tmp_prefix);
  permute_fixed_pgen(opt, tmp_prefix, opt.out);
  if (!opt.keep_temp) {
    remove_fileset(tmp_prefix);
    remove_pmerge_sidecars(tmp_prefix);
  }
  std::cerr << "output prefix: " << opt.out << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    auto opt = parse_args(argc, argv);
    merge_pgen(opt);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    usage(std::cerr);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
