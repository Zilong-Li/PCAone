#ifndef PCAONE_HALKO_
#define PCAONE_HALKO_

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData {
 public:
  Data* data;
  using Index = Eigen::Index;
  bool update = false, standardize = false, reset = true;
  Mat2D U, Omg, Omg2;  // nsamples x nk
  Mat2D V;       // nsnps x nk
  Mat1D S;       // nk x 1

 public:
  RsvdOpData(Data* data_) : data(data_) {}

  virtual ~RsvdOpData() {}

  virtual Index rows() const = 0;
  virtual Index cols() const = 0;
  virtual Index ranks() const = 0;
  virtual Index oversamples() const = 0;
  inline Index size() const { return ranks() + oversamples(); }

  virtual void computeGandH(Mat2D& G, Mat2D& H, int pi) = 0;

  /// update: whether to update the data in place
  /// standardize: whether to standardize the data in place
  /// reset: whether to reset the random matrix Omg
  inline void setFlags(bool is_update, bool is_standardize, bool reset_omg = true) {
    update = is_update;
    standardize = is_standardize;
    reset = reset_omg;
  }

  void computeUSV(int p, double tol);

  void initOmg();

  Mat2D computeU(const Mat2D& G, const Mat2D& H);
};

class NormalRsvdOpData : public RsvdOpData {
 private:
  const Index nk, os, size;
  uint64 actual_block_size, start_idx, stop_idx;

 public:
  NormalRsvdOpData(Data* data_, int k_, int os_ = 10) : RsvdOpData(data_), nk(k_), os(os_), size(k_ + os_) {
    initOmg();
  }

  ~NormalRsvdOpData() {}

  Index rows() const { return data->nsnps; }  // for snpmajor input
  Index cols() const { return data->nsamples; }
  Index ranks() const { return nk; }
  Index oversamples() const { return os; }

  void computeGandH(Mat2D& G, Mat2D& H, int pi = 0);
};

class FancyRsvdOpData : public RsvdOpData {
 private:
  const Index nk, os, size;
  uint64 bandsize, blocksize, actual_block_size, start_idx, stop_idx;
  Mat2D H1, H2;

 public:
  FancyRsvdOpData(Data* data_, int k_, int os_ = 10) : RsvdOpData(data_), nk(k_), os(os_), size(k_ + os_) {
    initOmg();
    H1 = Mat2D::Zero(cols(), size);
    H2 = Mat2D::Zero(cols(), size);
  }

  ~FancyRsvdOpData() {}

  Index rows() const { return data->nsnps; }
  Index cols() const { return data->nsamples; }
  Index ranks() const { return nk; }
  Index oversamples() const { return os; }

  void computeGandH(Mat2D& G, Mat2D& H, int pi = 0);
};

void run_pca_with_halko(Data* data, const Param& params);

#endif  // PCAONE_HALKO_
