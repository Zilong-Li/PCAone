#ifndef PCAONE_FILEBGEN_
#define PCAONE_FILEBGEN_

#include "Data.hpp"
#include "bgen/reader.h"
#include "bgen/writer.h"

// const double GENOTYPE_THRESHOLD = 0.9;
// const double BGEN_MISSING_VALUE = -9;
// const double BGEN2GENO[4] = {0, 0.5, 1, BGEN_MISSING_VALUE};


class FileBgen : public Data {
 public:
  // using Data::Data;
  FileBgen(const Param& params_) : Data(params_) {
    cao.print(tick.date(), "start parsing BGEN format");
    bg = new bgen::CppBgenReader(params.filein, "", true);
    nsamples = bg->header.nsamples;
    nsnps = bg->header.nvariants;
    dosages.resize(nsamples);
    F = Mat1D::Zero(nsnps);  // initial F
    cao.print(tick.date(), "N(#samples) =", nsamples, ", M(#SNPs) =", nsnps);
    cao.print(tick.date(), "the layout is", bg->header.layout, ", compressed by",
              bg->header.compression == 2 ? "zstd" : "zlib");
  }

  ~FileBgen() { delete bg; }

  virtual void read_all();
  // for blockwise
  virtual void check_file_offset_first_var() { bg->offset = bg->header.offset + 4; }

  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);

  virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals,
                                 const Mat2D& VT, bool standardize) {}

 private:
  bgen::CppBgenReader* bg;
  std::vector<float> dosages, probs1d;
  bool frequency_was_estimated = false;
};

void permute_bgen_thread(std::vector<int> idx, std::string fin, std::string fout, int ithread);

PermMat permute_bgen(std::string& fin, std::string fout, int nthreads);

#endif  // PCAONE_FILEBGEN_
