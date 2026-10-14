#ifndef PCAONE_HALKO_
#define PCAONE_HALKO_

#include "Data.hpp"
#include "LinAlg.hpp"

class RsvdOpData {
 public:
  Data* data;
  using Index = Eigen::Index;
  bool update = false, standardize = false;
  Mat2D U, Omg, Omg2;  // nsamples x nk
  Mat2D V;             // nsnps x nk
  Mat1D S;             // nk x 1

 public:
  RsvdOpData(Data* data_)
      : data(data_) {}

  virtual ~RsvdOpData() {}

  virtual Index rows() const = 0;
  virtual Index cols() const = 0;
  virtual Index ranks() const = 0;
  virtual Index oversamples() const = 0;
  inline Index size() const { return ranks() + oversamples(); }

  virtual void computeGandH(Mat2D& G, Mat2D& H, int pi) = 0;

  /// update: whether to update the data in place
  /// standardize: whether to standardize the data in place
  /// reset_omg: whether to reset the random matrix Omg
  inline void setFlags(bool is_update, bool is_standardize) {
    update = is_update;
    standardize = is_standardize;
  }

  void computeUSV(int p, double tol);

  // winEM (--winem): EMU's EM run inside the window updates of winSVD. Each
  // block of sites is imputed by Data::estep_reg() from the newest em_k PCs
  // (Uem) just before it enters the products, and each update of Omg (the
  // M-step) gives new PCs. Runs passes over the data until the PCs of two
  // passes agree to tol (1 - mev), but no fewer than log2(-w) + 3; U is then
  // the em_k PCs, and the object is left set for the final decomposition.
  void runWinEM(uint maxpass, double tol);
  // 0: the EM of --emu (fit_with_pi, read_block_update); 1: the passes of
  // winEM; 2: a decomposition of the matrix imputed from the fixed PCs Uem
  int winem = 0;
  Mat2D Uem;  // nsamples x em_k, orthonormal; empty before the first update

  void initOmg();
  void warmOmg();  // EM updates: start from the previous U instead of at random

 protected:
  // out-of-core: read the block of sites [start_idx, stop_idx], and add its
  // products: G's rows of the block = X_b' * Omg, then Hb += X_b * those rows
  void add_block(uint64 start_idx, uint64 stop_idx, Mat2D& G, Mat2D& Hb);

  // winEM: Uem = Q * (the em_k leading left singular vectors of R), the PCs of
  // H = QR just factored for the new Omg
  void update_em_basis(const TallQR& qr);
  // winEM: impute the block of sites [start_idx, start_idx + Z.cols()) of the
  // data, then standardize it if asked (out-of-core; in-core G is standardized whole)
  void estep_block(uint64 start_idx, Eigen::Ref<Mat2D> Z);

  int pi_offset = 0;  // > 0: a warm-started winSVD runs full bands from the first epoch
};

class NormalRsvdOpData : public RsvdOpData {
 private:
  const Index nk, os, size;

 public:
  NormalRsvdOpData(Data* data_, int k_, int os_ = 10)
      : RsvdOpData(data_),
        nk(k_),
        os(os_),
        size(k_ + os_) {
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
  FancyRsvdOpData(Data* data_, int k_, int os_ = 10)
      : RsvdOpData(data_),
        nk(k_),
        os(os_),
        size(k_ + os_) {
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
