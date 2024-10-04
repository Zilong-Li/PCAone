#ifndef PCAONE_HALKO_
#define PCAONE_HALKO_

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData {
 public:
  Data* data;
  using Index = Eigen::Index;
  bool update = false, standardize = false, verbose = false;
  Mat2D U, V;
  Mat1D S;

 public:
  RsvdOpData(Data* data_) : data(data_) {}

  virtual ~RsvdOpData() {}

  virtual Index rows() const = 0;
  virtual Index cols() const = 0;
  virtual Index ranks() const = 0;
  virtual Index oversamples() const = 0;

  virtual void computeGandH(Mat2D& G, Mat2D& H, int pi) = 0;

  /// update, standardize, verbose
  inline void setFlags(bool is_update, bool is_standardize, bool is_verbose) {
    update = is_update;
    standardize = is_standardize;
    verbose = is_verbose;
  }

  void computeUSV(int p, double tol);
};

class NormalRsvdOpData : public RsvdOpData {
 private:
  const Index nk, os, size;
  uint64 actual_block_size, start_idx, stop_idx;
  Mat2D Omg;

 public:
  NormalRsvdOpData(Data* data_, int k_, int os_ = 10)
      : RsvdOpData(data_), nk(k_), os(os_), size(k_ + os_) {}

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
  uint64 band, blocksize, actual_block_size, start_idx, stop_idx;
  Mat2D Omg, Omg2, H1, H2;

 public:
  FancyRsvdOpData(Data* data_, int k_, int os_ = 10)
      : RsvdOpData(data_), nk(k_), os(os_), size(k_ + os_) {
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

// void print_summary_table(const Mat2D& Upre, const Mat2D& Ucur);
void run_pca_with_halko(Data* data, const Param& params);

#endif  // PCAONE_HALKO_
