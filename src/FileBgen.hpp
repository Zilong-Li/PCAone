#ifndef PCAONE_FILEBGEN_
#define PCAONE_FILEBGEN_

#include "Data.hpp"
#include "bgen/reader.h"
#include "bgen/writer.h"

// const double GENOTYPE_THRESHOLD = 0.9;
// const double BGEN_MISSING_VALUE = -9;
// const double BGEN2GENO[4] = {0, 0.5, 1, BGEN_MISSING_VALUE};

void read_bgen_block(Mat2D& G, Mat1D& F, bgen::CppBgenReader* bg, float* dosages,
                     bool& frequency_was_estimated, uint64 nsamples, uint64 nsnps, uint blocksize,
                     uint64 start_idx, uint64 stop_idx, bool standardize);

int shuffle_bgen_to_bin(std::string& fin, std::string fout, uint gb, bool standardize);

void permute_bgen_thread(std::vector<int> idx, std::string fin, std::string fout, int ithread);

PermMat permute_bgen(std::string& fin, std::string fout, int nthreads);

class FileBgen : public Data {
 public:
  // using Data::Data;
  FileBgen(const Param& params_) : Data(params_) {
    cao.print(tick.date(), "start parsing BGEN format");
    bg = new bgen::CppBgenReader(params.filein, "", true);
    nsamples = bg->header.nsamples;
    nsnps = bg->header.nvariants;
    dosages.resize(nsamples);
    F = Mat1D::Zero(nsnps); // initial F
    cao.print(tick.date(), "the layout is", bg->header.layout, ", compression is", bg->header.compression , ". N (#samples):", nsamples,
              ". M (#SNPs):", nsnps);
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

#endif  // PCAONE_FILEBGEN_
