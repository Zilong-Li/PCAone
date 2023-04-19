#ifndef BGEN_UTILS_H_
#define BGEN_UTILS_H_

#include <bitset>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace bgen {

std::uint32_t n_choose_k(int n, int k);
bool minor_certain(double freq, int n_checked, double z);

} // namespace bgen

#endif  // BGEN_UTILS_H_
