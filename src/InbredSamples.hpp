#ifndef PCAONE_INBREDSAMPLES_
#define PCAONE_INBREDSAMPLES_

#include "Cmd.hpp"
#include "Data.hpp"

// --inbreed 2: the inbreeding coefficient of each sample accounting for
// population structure, i.e. with the individual allele frequencies pi of a
// reference PCA (-P/--USV) in place of one frequency per site. See
// InbredSamples.cpp for the estimator.
//
// Pi is the reference (a FileUSV). Writes <out>.inbred, one row per sample.
void run_inbred_samples(Data* Pi, const Param& params);

#endif  // PCAONE_INBREDSAMPLES_
