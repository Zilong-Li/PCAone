#ifndef BGEN_GENOTYPES_H_
#define BGEN_GENOTYPES_H_

#include <cstdint>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <sstream>

namespace bgen {

class Genotypes {
  std::ifstream* handle;
  int layout;
  int compression;
  int n_alleles;
  std::uint32_t n_samples;
  std::uint64_t offset;
  std::uint32_t length;
  std::uint32_t bit_depth;
  char * uncompressed;
  float * probs;
  float * alt_dose;
  float * minor_dose;
  bool is_decompressed = false;
  bool has_ploidy = false;
  bool probs_parsed = false;
  bool minor_dosage_parsed = false;
  bool alt_dosage_parsed = false;
  std::vector<int> missing;
public:
  Genotypes(std::ifstream* handle, int lay, int compr, int n_alleles, std::uint32_t n_samples, std::uint64_t offset, std::uint32_t length) :
     handle(handle), layout(lay), compression(compr), n_alleles(n_alleles), n_samples(n_samples), offset(offset), length(length) {};
  Genotypes() {};
  ~Genotypes() { clear_probs(); };
  void parse_preamble(char * uncompressed, std::uint32_t & idx);
  void parse_ploidy(char * uncompressed, std::uint32_t & idx);
  float * parse_layout1(char *, std::uint32_t & idx);
  float * parse_layout2(char *, std::uint32_t & idx);
  void decompress();
  float * probabilities();
  int find_minor_allele(float * dose);
  float * get_allele_dosage(bool use_alt=true, bool use_minor=false);
  void ref_dosage_fast(char * uncompressed, std::uint32_t & idx, float * dose);
  void ref_dosage_fast_fallback(char * uncompressed, std::uint32_t & idx, float * dose);
  void swap_allele_dosage(float * dose);
  void ref_dosage_slow(char * uncompressed, std::uint32_t & idx, float * dose);
  void clear_probs();
  bool phased;
  std::uint32_t max_probs = 0;
  bool constant_ploidy;
  int min_ploidy;
  int max_ploidy;
  int minor_idx;
  std::uint8_t * ploidy;
};

std::uint32_t get_max_probs(int &max_ploidy, int &n_alleles, bool &phased);

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
