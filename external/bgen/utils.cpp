
#include "utils.h"

namespace bgen {

// Returns value of Binomial Coefficient C(n, k)
uint n_choose_k(int n, int k) {
  uint res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k ) {
    k = n - k;
  }

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (uint i = 0; i < (uint) k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

/// check if the minor allele is certain (to 99.9999999999999& confidence)
///
///  Take the frequency, and number of individuals checked so far, and see if the
///  99.99..(fifteen nines) confidence interval overlaps 0.5. If not, then we can
///  be sure we've identified the minor allele, even without checking the full
///  population.
///
///  @param freq estimated minor allele frequency
///  @param n_checked number of individsuals checked so far
///  @param z standard normal deviate (eg 1.96 for 95% CI, here we use 10.0 for
///    stronger confidence, and the fact the normal approximation for confidence
///    intervals isn't perfect)
///  @return True/False for whether to halt the permuations
bool minor_certain(double freq, int n_checked, double z) {
    double delta = (z * std::sqrt((freq * (1 - freq)) / n_checked));
    
    // check if the confidence interval overlaps 0.5
    return !((freq - delta < 0.5) & (freq + delta > 0.5));
}


} // namespace bgen
