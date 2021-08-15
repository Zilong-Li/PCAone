
#include <iterator>
#include <algorithm>
#include <cassert>
#include <bitset>

#include "header.h"

namespace bgen {

Header::Header(std::ifstream & handle) {
  handle.seekg(0);
  char buff[20];
  handle.read(&buff[0], 20);
  
  offset = *reinterpret_cast<const std::uint32_t*>(&buff[0]);
  header_length = *reinterpret_cast<const std::uint32_t*>(&buff[4]);
  nvariants = *reinterpret_cast<const std::uint32_t*>(&buff[8]);
  nsamples = *reinterpret_cast<const std::uint32_t*>(&buff[12]);
  magic = std::string(&buff[16], 4);
  
  // make sure we are reading a bgen file
  if ((magic != "bgen") & ((int) (magic[0] & magic[1] & magic[2] & magic[3]) != 0)) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  
  // read any extra data contained in the header
  int size = header_length - 20;
  if (size > 0) {
    std::copy_n(std::istream_iterator<char>(handle), size, std::back_inserter(extra));
  }
  
  // read flags data
  std::bitset<32> flags;
  handle.read(reinterpret_cast<char*>(&flags), sizeof(std::uint32_t));
  
  std::bitset<32> compr_mask(0b000000000000000000000000000011);
  std::bitset<32> layout_mask(0b000000000000000000000000111100);
  compression = (int) (flags & compr_mask).to_ulong();
  layout = (int) ((flags & layout_mask) >> 2).to_ulong();
  has_sample_ids = flags[31];
}

} // namespace bgen
