#ifndef PCAONE_HALKO_
#define PCAONE_HALKO_

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData {
 public:
  Data* data;
  using Index = Eigen::Index;
  bool update = false, standardize = false, verbose = false;
  MyMatrix U, V;
  MyVector S;

 public:
  RsvdOpData(Data* data_) : data(data_) {}

  virtual ~RsvdOpData() {}

  virtual Index rows() const = 0;
  virtual Index cols() const = 0;
  virtual Index ranks() const = 0;
  virtual Index oversamples() const = 0;

  virtual void computeGandH(MyMatrix& G, MyMatrix& H, int pi) = 0;

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
  MyMatrix Omg;

 public:
  NormalRsvdOpData(Data* data_, int k_, int os_ = 10)
      : RsvdOpData(data_), nk(k_), os(os_), size(k_ + os_) {}

  ~NormalRsvdOpData() {}

  Index rows() const { return data->nsnps; }  // for snpmajor input
  Index cols() const { return data->nsamples; }
  Index ranks() const { return nk; }
  Index oversamples() const { return os; }

  void computeGandH(MyMatrix& G, MyMatrix& H, int pi = 0);
};

class FancyRsvdOpData : public RsvdOpData {
 private:
  const Index nk, os, size;
  uint64 band, blocksize, actual_block_size, start_idx, stop_idx;
  MyMatrix Omg, Omg2, H1, H2;

 public:
  FancyRsvdOpData(Data* data_, int k_, int os_ = 10)
      : RsvdOpData(data_), nk(k_), os(os_), size(k_ + os_) {
    H1 = MyMatrix::Zero(cols(), size);
    H2 = MyMatrix::Zero(cols(), size);
  }

  ~FancyRsvdOpData() {}

  Index rows() const { return data->nsnps; }
  Index cols() const { return data->nsamples; }
  Index ranks() const { return nk; }
  Index oversamples() const { return os; }

  void computeGandH(MyMatrix& G, MyMatrix& H, int pi = 0);
};

// void print_summary_table(const MyMatrix& Upre, const MyMatrix& Ucur);
void run_pca_with_halko(Data* data, const Param& params);

#endif  // PCAONE_HALKO_
