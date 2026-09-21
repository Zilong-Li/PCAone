#ifndef PCAONE_FILEUSV
#define PCAONE_FILEUSV

#include <cmath>

#include "Data.hpp"
#include "Utils.hpp"

// To get \Pi mainly
class FileUSV : public Data {
 public:
  FileUSV(const Param& params_)
      : Data(params_) {
    cao.print(tick.date(), "start parsing U:", params.fileU, ", S:", params.fileS, ", V:", params.fileV);
    read_sigvals(params.fileS, nsamples, nsnps, S, &usv);  // could not structual bindings
    if (S.size() != params.k) cao.warn("the value of -k not equal the number of rows in " + params.fileS);
    K = fmin(S.size(), params.k);
    cao.print(tick.date(), "start parsing mbim and read allele frequency of SNPs from", params.filebim);
    F = read_frq(params.filebim);
    if (F.size() != nsnps) cao.error("the number of sites in mbim not matching the header line of .sigvals");
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps);
    V = read_eigvecs(params.fileV, nsnps, K);
    U = read_eigvecs(params.fileU, nsamples, K);
    if (params.inbreed) check_transform();
  }

  ~FileUSV() override = default;

  void read_all() final;

  // for blockwise
  void check_file_offset_first_var() final {}

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D&, const Mat1D&, const Mat2D&, bool) final {}

  // factor mapping the reconstruction back onto the 0..1 allele-frequency
  // scale: pi = U*S*V' * inv_scale(f) + f. See check_transform().
  double inv_scale(double f) const {
    if (usv.gscale == 2) return 0.5;                              // centred dosages, 0..2
    if (usv.scale != SCALE_STANDARDIZE_GENETIC) return 1.0;       // centred only, 0..1
    const double sd = std::sqrt(f * (1.0 - f));                   // standardized, 0..1
    return (sd > VAR_TOL) ? sd / std::sqrt((double)usv.ploidy) : 1.0;
  }

 private:
  void check_transform();

  int K;
  Mat1D S;
  Mat2D U, V;
  UsvTransform usv;
};

#endif  // PCAONE_FILEUSV
