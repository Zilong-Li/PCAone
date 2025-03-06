
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
///  @param _handle std::ifstream for bgen file
///  @param varoffset start byte for variant in bgen file
///  @param layout bgen layout version (1 or 2)
///  @param compression compression scheme (0=no compression, 1=zlib, 2=zstd)
///  @param expected_n number of samples for variant
Variant::Variant(std::ifstream * _handle, std::uint64_t & varoffset, int layout, int compression, int expected_n, std::uint64_t fsize) : handle(_handle) {
  if (varoffset >= fsize) {
    throw std::out_of_range("reached end of file");
  }
  offset = varoffset;
  handle->seekg(offset);
  if (layout == 1) {
    handle->read(reinterpret_cast<char*>(&n_samples), sizeof(n_samples));
  } else {
    n_samples = expected_n;
  }
  
  if ((int) n_samples != expected_n) {
    throw std::invalid_argument("number of samples doesn't match");
  }
  
  // get the variant ID (first need to know how long the field is)
  std::uint16_t item_len;
  handle->read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    varid.resize(item_len);
    handle->read(&varid[0], item_len);
  }
  
  // get the rsID (first need to know how long the field is)
  handle->read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    rsid.resize(item_len);
    handle->read(&rsid[0], item_len);
  }
  
  // get the chromosome (first need to know how long the field is)
  handle->read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    chrom.resize(item_len);
    handle->read(&chrom[0], item_len);
  }
  
  handle->read(reinterpret_cast<char*>(&pos), sizeof(std::uint32_t));
  if (layout == 1) {
    n_alleles = 2;
  } else {
    handle->read(reinterpret_cast<char*>(&n_alleles), sizeof(std::uint16_t));
  }
  
  for (int x=0; x < n_alleles; x++) {
    std::uint32_t allele_len;
    std::string allele;
    handle->read(reinterpret_cast<char*>(&allele_len), sizeof(std::uint32_t));
    allele.resize(allele_len);
    handle->read(&allele[0], allele_len);
    alleles.push_back(allele);
  }
  
  std::uint32_t length;
  if ((layout == 1) && (compression == 0)) {
    length = n_samples * 6;
  } else {
    handle->read(reinterpret_cast<char *>(&length), sizeof(length));
  }
  std::uint64_t geno_offset = (std::uint64_t) handle->tellg();
  geno = Genotypes(handle, layout, compression, n_alleles, n_samples, geno_offset, length);
  next_variant_offset = geno_offset + length;
}

int Variant::probs_per_sample() {
  if (geno.max_probs == 0) {
    geno.load_data_and_parse_header();
  }
  return geno.max_probs;
}

bool Variant::phased() {
  if (geno.max_probs == 0) {
    geno.load_data_and_parse_header();
  }
  return geno.phased;
}

std::uint8_t * Variant::ploidy() {
  if (geno.max_probs == 0) {
    geno.load_data_and_parse_header();
  }
  return geno.ploidy;
}

/// get genotype probabilities for the variant as a 1-dimensional vector
///
/// This makes it easy to pass the data via cython into a numpy array, which can
/// be reshaped to a 2-D array.
void Variant::probs_1d(float * probs) {
  geno.probabilities(probs);
}

/// get dosage of the alt allele (only works for biallelic variants)
void Variant::alt_dosage(float * dose) {
  geno.get_allele_dosage(dose, true, false);
  minor_allele = alleles[geno.minor_idx];
}

/// get dosage of the minor allele (only works for biallelic variants)
void Variant::minor_allele_dosage(float * dose) {
  geno.get_allele_dosage(dose, false, true);
  minor_allele = alleles[geno.minor_idx];
}

std::vector<std::uint8_t> Variant::copy_data() {
  std::uint32_t length = next_variant_offset - offset;
  std::vector<std::uint8_t> data(length);
  handle->seekg(offset);
  handle->read(reinterpret_cast<char *>(&data[0]), length);
  return data;
}

} // namespace bgen
