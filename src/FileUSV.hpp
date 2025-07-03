#ifndef PCAONE_FILEUSV
#define PCAONE_FILEUSV

#include "Data.hpp"

// To get \Pi mainly
class FileUSV : public Data {
 public:
  FileUSV(const Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing U:", params.fileU, ", S:", params.fileS, ", V:", params.fileV);
    read_sigvals(params.fileS, nsamples, nsnps, S); // could not structual bindings
    if (S.size() != params.k) cao.warn("the value of -k not equal the number of rows in " + params.fileS);
    K = fmin(S.size(), params.k);
    cao.print(tick.date(), "start parsing mbim and read allele frequency of SNPs from", params.filebim);
    F = read_frq(params.filebim);
    if (F.size() != nsnps) cao.error("the number of sites in mbim not matching the header line of .sigvals");
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps);
    V = read_eigvecs(params.fileV, nsnps, K);
    U = read_eigvecs(params.fileU, nsamples, K);
  }

  ~FileUSV() override = default;

  void read_all() final;
  
  // for blockwise
  void check_file_offset_first_var() final {}

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final {}

 private:
  int K;
  Mat1D S;
  Mat2D U, V;
};

#endif  // PCAONE_FILEUSV
