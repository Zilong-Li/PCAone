#ifndef PCAONE_FILEBEAGLE_
#define PCAONE_FILEBEAGLE_

#include <zlib.h>

#include "Data.hpp"

int tgets(gzFile gz, char **buf, uint64 *l);

class FileBeagle : public Data {
 public:
  FileBeagle(Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing BEAGLE format");
    original = buffer = (char *)calloc(bufsize, sizeof(char));
    if (params.nsnps > 0 && params.nsamples > 0) {
      cao.print(tick.date(), "use nsamples and nsnps given by user");
      nsamples = params.nsamples;
      nsnps = params.nsnps;
    } else {
      fp = gzopen(params.filein.c_str(), "r");
      tgets(fp, &buffer, &bufsize);
      int nCol = 1;
      if (buffer != original) original = buffer;
      strtok_r(buffer, delims, &buffer);
      while (strtok_r(NULL, delims, &buffer)) nCol++;
      if (nCol % 3) cao.error("Number of columns should be a multiple of 3.");
      nsamples = nCol / 3 - 1;
      // continue getting the number of sites
      // FIXME:  assume the columns are aligne. maybe check it first.
      buffer = original;
      nsnps = 0;
      while (tgets(fp, &buffer, &bufsize)) nsnps++;
      gzclose(fp);
    }

    P = Mat2D::Zero(nsamples * 2, nsnps);  // Mat2D is column major
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps);
  }

  ~FileBeagle() {}

  virtual void read_all();
  // below are for blockwise, remain for future.
  virtual void check_file_offset_first_var() {}

  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {}

  virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D &U, const Mat1D &svals,
                                 const Mat2D &VT, bool standardize) {}

 private:
  gzFile fp = nullptr;
  char *original, *buffer;
  uint64 bufsize = (uint64)128 * 1024 * 1024;
  const char *delims = "\t \n";
};

#endif  // PCAONE_FILEBEAGLE_
