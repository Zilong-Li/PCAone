#ifndef PCAONE_ARNOLDI_
#define PCAONE_ARNOLDI_

#include "Data.hpp"

class ArnoldiOpData {
 public:
  ArnoldiOpData(Data* data_) : data(data_), n(data_->nsamples) {
    data->nops = 1;
  }

  ~ArnoldiOpData() {}

  // The line below is new for spectra v1.0.0
  using Scalar = double;

  inline uint64 rows() const { return n; }
  inline uint64 cols() const { return n; }
  // y = G * G' * x ; data.G is n x m;
  void perform_op(const double* x_in, double* y_out) const;
  inline void setFlags(bool is_update, bool is_standardize, bool is_pcangsd) {
    update = is_update;
    standardize = is_standardize;
    pcangsd = is_pcangsd;
  }

  Mat2D U, VT;
  Mat1D S;

 private:
  Data* data;
  const uint64 n;
  bool update = false, standardize = false, pcangsd = false;
};

void run_pca_with_arnoldi(Data* data, const Param& params);

#endif  // PCAONE_ARNOLDI_
