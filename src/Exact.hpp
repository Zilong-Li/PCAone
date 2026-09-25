#ifndef PCAONE_EXACT_H
#define PCAONE_EXACT_H

#include "Cmd.hpp"
#include "Data.hpp"

// --svd 3. With params.out_of_core set, the N x N GRM is accumulated block by
// block, so G is never held whole: memory is N^2 doubles plus one block. Main
// sets it for N <= M (unless --maf, which only the in-core path supports).
// Otherwise the in-core path decomposes GG' (N <= M) or G'G (N > M).
void run_pca_exact(Data* data, const Param& params);

#endif  // PCAONE_EXACT_H
