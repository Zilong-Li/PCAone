#ifndef BGEN_SAMPLES_H_
#define BGEN_SAMPLES_H_

#include <fstream>
#include <vector>
#include <string>

namespace bgen {

class Samples {
public:
  Samples(std::ifstream & handle, int n_samples);
  Samples(std::string path, int n_samples);
  Samples(int n_samples);
  Samples() {};
  std::vector<std::string> samples;
};

} // namespace bgen

#endif  // BGEN_SAMPLES_H_
