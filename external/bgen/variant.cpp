
#include <stdexcept>
#include <cmath>

#include "variant.h"

namespace bgen {

/// initialise a single variant with chrom, pos, rsID identifiers
///
/// This starts a Genotypes object, but this doesn't parse the genotypes until
/// required, just starts it so we can get the offset of the next variant, so as
/// to parse the bgen variants at speed.
///
///  @param handle std::ifstream for bgen file
///  @param offset start byte for variant in bgen file
///  @param layout bgen layout version (1 or 2)
///  @param compression compression scheme (0=no compression, 1=zlib, 2=zstd)
///  @param expected_n number of samples for variant
Variant::Variant(std::ifstream & handle, std::uint64_t & varoffset, int layout, int compression, int expected_n) {
  offset = varoffset;
  handle.seekg(offset);
  if (layout == 1) {
    handle.read(reinterpret_cast<char*>(&n_samples), sizeof(n_samples));
  } else {
    n_samples = expected_n;
  }
  
  if ((int) n_samples != expected_n) {
    throw std::invalid_argument("number of samples doesn't match");
  }
  
  // get the variant ID (first need to know how long the field is)
  std::uint16_t item_len;
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(varid));
  }
  
  // get the rsID (first need to know how long the field is)
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(rsid));
  }
  
  // get the chromosome (first need to know how long the field is)
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(chrom));
  }
  
  handle.read(reinterpret_cast<char*>(&pos), sizeof(std::uint32_t));
  if (layout == 1) {
    n_alleles = 2;
  } else {
    handle.read(reinterpret_cast<char*>(&n_alleles), sizeof(std::uint16_t));
  }
  
  for (int x=0; x < n_alleles; x++) {
    std::uint32_t allele_len;
    std::string allele;
    handle.read(reinterpret_cast<char*>(&allele_len), sizeof(std::uint32_t));
    std::copy_n(std::istream_iterator<char>(handle), allele_len, std::back_inserter(allele));
    alleles.push_back(allele);
  }
  
  geno = Genotypes(&handle, layout, compression, n_alleles, n_samples);
}

/// uses the genotypes object to find the offset of the next variant
std::uint64_t Variant::next_variant_offset() {
  return geno.next_var_offset;
}

int Variant::probs_per_sample() {
  return geno.max_probs;
}

bool Variant::phased() {
  if (geno.max_probs == 0) {
    throw std::invalid_argument("unknown phase, run variant.probabilities() first");
  }
  return geno.phased;
}

std::uint8_t * Variant::ploidy() {
  if (geno.max_probs == 0) {
    throw std::invalid_argument("unknown ploidy, run variant.probabilities() first");
  }
  return geno.ploidy;
}

/// get genotype probabilities for the variant as a 1-dimensional vector
///
/// This makes it easy to pass the data via cython into a numpy array, which can
/// be reshaped to a 2-D array.
float * Variant::probs_1d() {
  return geno.probabilities();
}

/// get dosage of the minor allele (only works for biallelic variants)
float * Variant::minor_allele_dosage() {
  float * dose = geno.minor_allele_dosage();
  minor_allele = alleles[geno.minor_idx];
  return dose;
}

} // namespace bgen
