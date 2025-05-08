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

  ~FileBgen() override { delete bg; }

  void read_all() final;

  // for blockwise
  void check_file_offset_first_var() final { bg->offset = bg->header.offset + 4; }

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final {}

 private:
  bgen::CppBgenReader* bg;
  std::vector<float> dosages, probs1d;
  bool frequency_was_estimated = false;
};

void permute_bgen_thread(std::vector<int> idx, std::string fin, std::string fout, int ithread);

PermMat permute_bgen(std::string& fin, std::string fout, int nthreads);

#endif  // PCAONE_FILEBGEN_
