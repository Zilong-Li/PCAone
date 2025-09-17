#ifndef PCAONE_FILEBEAGLE_
#define PCAONE_FILEBEAGLE_

#include "Data.hpp"

class FileBeagle : public Data {
 public:
  FileBeagle(const Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing BEAGLE format");
    // TODO: ask user to use --pcangsd explicitly
    impute = true; // always imputing imcomplete information
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
      const char *delims = "\t \n";
      strtok_r(buffer, delims, &buffer);
      while (strtok_r(NULL, delims, &buffer)) nCol++;
      if (nCol % 3) cao.error("Number of columns should be a multiple of 3.");
      nsamples = nCol / 3 - 1;
      // NOTE:  assume the columns are alignd. maybe check it first.
      buffer = original;
      nsnps = 0;
      // continue getting the number of sites
      while (tgets(fp, &buffer, &bufsize)) {
        nsnps++;
      } 
    }

    if (params.pca) {  // initial F
      F = Mat1D::Zero(nsnps);
    }
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps);
  }

  ~FileBeagle() override {
    free(buffer);
    if (fp) gzclose(fp);
  }

  void read_all() final;
  
  // below are for blockwise
  void check_file_offset_first_var() final;

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final {}

 private:
  gzFile fp = nullptr;
  char *original, *buffer;
  uint64 bufsize = (uint64)128 * 1024 * 1024;
};

#endif  // PCAONE_FILEBEAGLE_
