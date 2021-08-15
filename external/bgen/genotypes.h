#ifndef BGEN_GENOTYPES_H_
#define BGEN_GENOTYPES_H_

#include <cstdint>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <sstream>

namespace bgen {

typedef unsigned int uint;

class Genotypes {
  std::uint64_t offset;
  std::ifstream* handle;
  int layout;
  int compression;
  int n_alleles;
  uint n_samples;
  uint bit_depth;
  char * uncompressed;
  float * probs;
  float * dose;
  bool is_decompressed = false;
  bool has_ploidy = false;
  bool probs_parsed = false;
  bool dosage_parsed = false;
  std::vector<int> missing;
public:
  Genotypes(std::ifstream* handle, int lay, int compr, int n_alleles, uint n_samples) :
    handle(handle), layout(lay), compression(compr), n_alleles(n_alleles), n_samples(n_samples) {
      std::uint32_t length;
      handle->read(reinterpret_cast<char*>(&length), sizeof(length));
      offset = handle->tellg();
      next_var_offset = offset + length;
  };
  Genotypes() {};
  ~Genotypes() { clear_probs(); };
  void parse_preamble(char * uncompressed, uint & idx);
  void parse_ploidy(char * uncompressed, uint & idx);
  float * parse_layout1(char *, uint & idx);
  float * parse_layout2(char *, uint & idx);
  void decompress();
  float * probabilities();
  int find_minor_allele(float * dose);
  float * minor_allele_dosage();
  void ref_dosage_fast(char * uncompressed, uint & idx);
  void alt_dosage();
  void ref_dosage_slow(char * uncompressed, uint & idx);
  void clear_probs();
  bool phased;
  uint max_probs = 0;
  bool constant_ploidy;
  int min_ploidy;
  int max_ploidy;
  int minor_idx;
  std::uint8_t * ploidy;
  std::uint64_t next_var_offset;
};

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
