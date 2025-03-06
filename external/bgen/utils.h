#ifndef BGEN_UTILS_H_
#define BGEN_UTILS_H_

#include <bitset>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace bgen {

struct Range {
    std::uint8_t _min;
    std::uint8_t _max;
};

std::uint32_t n_choose_k(int n, int k);
bool minor_certain(double freq, int n_checked, double z);
std::uint64_t fast_ploidy_sum(std::uint8_t * x, std::uint32_t & size);
Range fast_range(std::uint8_t * x, std::uint32_t & size);

} // namespace bgen

#endif  // BGEN_UTILS_H_
