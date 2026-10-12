/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Variants.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 *
 * The variants of the input (.bim, .pvar, BEAGLE markers) and the .mbim of a
 * reference PCA, which is a .bim with the allele frequency in column 7.
 ******************************************************************************/
#ifndef PCAONE_VARIANTS_
#define PCAONE_VARIANTS_

#include <istream>

#include "Common.hpp"

class Param;

struct BimMatch {
  bool identical = false;  // the same variants in the same order
  Int1D bim_indices;       // the input's variants found in the .mbim
  Int1D mbim_indices;      // and where
  std::vector<bool> flip;  // parallel to bim_indices: true if ref/alt are swapped vs mbim
};

// The variants of the input (.bim, .pvar or BEAGLE) matched to those of the
// .mbim params.filebim by chromosome, position and alleles, either way round.
// Stops with an error when none match.
BimMatch match_variants_to_mbim(const Param& params);

// a .pvar line as the .bim line of the same variant: CHR ID 0 POS ALT REF, with
// the first ALT only. The PGEN reader counts ALT, so ALT is A1, as BED2GENO
// counts A1 of a .bim.
std::string pvar_line_to_bim_line(const std::string& line, const std::string& path);

// the next variant of a .bim, or of a .pvar (pvar = true) as its .bim line,
// the header lines skipped. false at the end of the file
bool read_bim_line(std::istream& is, bool pvar, const std::string& path, std::string& line);

#endif  // PCAONE_VARIANTS_
