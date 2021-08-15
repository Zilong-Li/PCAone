#ifndef BGEN_BGEN_H_
#define BGEN_BGEN_H_

#include <fstream>
#include <stdexcept>
#include <vector>

#include "header.h"
#include "samples.h"
#include "variant.h"

namespace bgen {

class Bgen {
  std::ifstream handle;
  std::uint64_t fsize;
  std::uint64_t offset;
public:
  Bgen(std::string path, std::string sample_path="", bool delay_parsing=false);
  void parse_all_variants();
  Variant next_var();
  void drop_variants(std::vector<int> indices);
  std::vector<std::string> varids();
  std::vector<std::string> rsids();
  std::vector<std::string> chroms();
  std::vector<std::uint32_t> positions();
  Variant & operator[](std::size_t idx) { return variants[idx]; }
  Variant & get(std::size_t idx) { return variants[idx]; };
  std::vector<Variant> variants;
  Header header;
  Samples samples;
};

} // namespace bgen

#endif  // BGEN_BGEN_H_
