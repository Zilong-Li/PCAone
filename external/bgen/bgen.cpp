
#include "bgen.h"

namespace bgen {

Bgen::Bgen(std::string path, std::string sample_path, bool delay_parsing) {
  handle.open(path, std::ios::binary);
  if (!handle) {
    throw std::invalid_argument("error reading from '" + path + "'");
  }
  fsize = handle.tellg();
  header = Header(handle);
  if (header.has_sample_ids) {
    samples = Samples(handle, header.nsamples);
  } else if (sample_path.size() > 0) {
    samples = Samples(sample_path, header.nsamples);
  } else {
    samples = Samples(header.nsamples);
  }
  
  // figure out the file length, so we don't go beyond it
  handle.seekg(0, std::ios::end);
  fsize = (std::uint64_t) handle.tellg() - fsize;
  
  offset = header.offset + 4;
  if (!delay_parsing) {
    parse_all_variants();
  }
}

Variant Bgen::next_var() {
  if (handle.eof() | (offset >= fsize)) {
    throw std::out_of_range("reached end of file");
  }
  Variant var(handle, offset, header.layout, header.compression, header.nsamples);
  offset = var.next_variant_offset();
  return var;
}

/// load all variants into memory at once
void Bgen::parse_all_variants() {
  offset = header.offset + 4;
  variants.clear();
  variants.resize(header.nvariants);
  int idx = 0;
  while (true) {
    try {
      variants[idx] = next_var();
      idx += 1;
    } catch (const std::out_of_range & e) {
      break;
    }
  }
  // finally reset the offset position to the first variant, so we can iterate
  // over variants more easily in python
  offset = header.offset + 4;
}

/// drop a subset of variants passed in by indexes
void Bgen::drop_variants(std::vector<int> indices) {
  // sort indices in descending order, so dropping elemtns doesn't affect later items
  std::sort(indices.rbegin(), indices.rend());
  
  auto it = std::unique(indices.begin(), indices.end());
  if (it != indices.end()) {
    throw std::invalid_argument("can't drop variants with duplicate indices");
  }
  
  for (auto idx : indices) {
    variants[idx] = variants.back();
    variants.pop_back();
  }
  variants.shrink_to_fit();
  
  // and sort the variants again afterward
  std::sort(variants.begin(), variants.end(),
          [] (Variant const& a, Variant const& b) { return a.pos < b.pos; });
}

/// get all the IDs for the variants in the bgen file
std::vector<std::string> Bgen::varids() {
  std::vector<std::string> varid(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    varid[x] = variants[x].varid;
  }
  return varid;
}

/// get all the rsIDs for the variants in the bgen file
std::vector<std::string> Bgen::rsids() {
  std::vector<std::string> rsid(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    rsid[x] = variants[x].rsid;
  }
  return rsid;
}

/// get all the chroms for the variants in the bgen file
std::vector<std::string> Bgen::chroms() {
  std::vector<std::string> chrom(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    chrom[x] = variants[x].chrom;
  }
  return chrom;
}

/// get all the positions for the variants in the bgen file
std::vector<std::uint32_t> Bgen::positions() {
  std::vector<std::uint32_t> position(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    position[x] = variants[x].pos;
  }
  return position;
}

} // namespace bgen
