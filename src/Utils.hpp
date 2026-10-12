/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Utils.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2026. Use of this code is governed by the LICENSE file.
 *
 * The logger and timer every module shares, and small helpers for strings and
 * files. The rest of what lived here has its own module:
 *   LinAlg.hpp    threaded products, flip_UV, convergence measures
 *   Stats.hpp     chi-squared distribution, median, allele frequencies from GLs
 *   TextIO.hpp    fast text output, sample IDs, the files of a PCA
 *   Variants.hpp  .bim/.pvar/BEAGLE variants and their match to a .mbim
 ******************************************************************************/
#ifndef PCAONE_UTILES_
#define PCAONE_UTILES_

#include "Common.hpp"
#include "Logger.hpp"
#include "Timer.hpp"
#include "zlib.h"

// MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
Logger cao;  // logger
Timer tick;  // Timer
#else
extern Timer tick;
extern Logger cao;
#endif

/// get the machine information
std::string get_machine();

// split s at every one of the separators; consecutive separators give empty fields
String1D split_string(const std::string& s, const std::string& separators);

size_t count_lines(const std::string& fpath);

/* @brief gets a line of gzfile
 * @param gz   file hander returned by gzopen
 * @param buf  buffer used for storing data
 * @param size buffer size for realloc buffer
 * @return extended buf length
 */
int tgets(gzFile gz, char** buf, uint64* size);

FILE* fopenOrDie(const char* filename, const char* instruction);

size_t freadOrDie(void* buffer, size_t sizeToRead, FILE* file);

void fcloseOrDie(FILE* file);

bool isZstdCompressed(const char* filename);

// zstd decompression buffer
struct ZstdDS {
  ZstdDS() {
    buffInTmp.reserve(buffInSize);
    buffOutTmp.reserve(buffOutSize);
  }
  ~ZstdDS() {
    ZSTD_freeDCtx(dctx);
    fcloseOrDie(fin);
  }
  FILE* fin = nullptr;
  size_t const buffInSize = ZSTD_DStreamInSize();
  size_t const buffOutSize = ZSTD_DStreamOutSize();
  ZSTD_DCtx* const dctx = ZSTD_createDCtx();
  size_t lastRet = 1;
  std::string buffInTmp, buffOutTmp;
};

// the EM-PCA loops ran out of --maxiter without a word
template <typename P>
void warn_em_not_converged(const P& params, double diff) {
  cao.warn("EM-PCA did not converge in --maxiter " + std::to_string(params.maxiter) +
           " iterations (diff = " + std::to_string(diff) + ", --tol-em = " + std::to_string(params.tolem) + ")");
}

#endif  // PCAONE_UTILES_
