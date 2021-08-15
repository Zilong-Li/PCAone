
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>

#include "samples.h"
#include "utils.h"

namespace bgen {

Samples::Samples(std::ifstream & handle, int n_samples) {
  /* initialize sample list if present in the bgen file
  */
  std::uint32_t sample_header_length;
  handle.read(reinterpret_cast<char*>(&sample_header_length), sizeof(std::uint32_t));
  
  std::uint32_t sample_n_check;
  handle.read(reinterpret_cast<char*>(&sample_n_check), sizeof(std::uint32_t));
  if (n_samples != (int) sample_n_check) {
    throw std::invalid_argument("inconsistent number of samples");
  }
  
  samples.resize(n_samples);
  std::uint16_t id_len;
  for (int i=0; i<n_samples; i++) {
    handle.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
    std::string sample_id;
    std::copy_n(std::istream_iterator<char>(handle), id_len, std::back_inserter(sample_id));
    samples[i] = sample_id;
  }
}

Samples::Samples(std::string path, int n_samples) {
  /* initialize from external sample file
  */
  std::ifstream handle(path, std::ios::in);
  if (!handle) {
    throw std::invalid_argument("error with sample file: '" + path + "'");
  }
  
  std::string header;
  std::getline(handle, header);
  std::string types;
  std::getline(handle, types);
  
  // find the file length post header, then read it all in to memory
  auto pos = handle.tellg();
  handle.seekg(0, std::ios::end);
  auto fsize = (std::uint64_t) handle.tellg() - pos;
  std::string lines(fsize, '\0');
  handle.seekg(pos);
  handle.read(&lines[0], fsize);
  
  samples.resize(n_samples);
  std::string sample_id;
  std::istringstream iss(lines);
  
  // run through all lines and gte the first column as sample_id
  int idx = 0;
  std::string line;
  while (std::getline(iss, line, '\n')) {
    samples[idx] = line.substr(0, line.find(' '));
    idx += 1;
  }
  
  if (idx != n_samples) {
    throw std::invalid_argument("inconsistent number of samples");
  }
  
  if (n_samples != (int)samples.size()) {
    throw std::invalid_argument("inconsistent number of samples");
  }
}

Samples::Samples(int n_samples) {
  /* initialize with integer IDs if no sample list available
  */
  samples.resize(n_samples);
  for (int i=0; i<n_samples; i++) {
    samples[i] = std::to_string(i);
  }
}

} // namespace bgen
