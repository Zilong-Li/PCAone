/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/RSVD.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#ifndef PCAONE_RSVD_
#define PCAONE_RSVD_

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>
#include <stdexcept>

namespace PCAone {

template <typename MatrixType>
using RealType = typename Eigen::NumTraits<typename MatrixType::Scalar>::Real;

template <typename MatrixType, typename ScalarType, typename RandomEngineType>
struct StandardNormalRandomHelper {
  static inline MatrixType generate(Eigen::Index numRows, Eigen::Index numCols, RandomEngineType& engine);
};

template <typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, RealType<MatrixType>, RandomEngineType> {
  static inline MatrixType generate(const Eigen::Index numRows, const Eigen::Index numCols,
                                    RandomEngineType& engine) {
    std::normal_distribution<RealType<MatrixType>> distribution{0, 1};
    const auto normal{[&](typename MatrixType::Scalar) { return distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, normal);
  }
};

template <typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, std::complex<RealType<MatrixType>>, RandomEngineType> {
  static inline MatrixType generate(const Eigen::Index numRows, const Eigen::Index numCols,
                                    RandomEngineType& engine) {
    constexpr RealType<MatrixType> stdDev{0.707106781186547};
    std::normal_distribution<RealType<MatrixType>> distribution{0, stdDev};
    const auto complexNormal{[&](typename MatrixType::Scalar) {
      return std::complex<RealType<MatrixType>>(distribution(engine), distribution(engine));
    }};
    return MatrixType::NullaryExpr(numRows, numCols, complexNormal);
  }
};

template <typename MatrixType, typename RandomEngineType>
inline MatrixType UniformRandom(const Eigen::Index numRows, const Eigen::Index numCols,
                                RandomEngineType& engine) {
  std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{-1, 1};
  const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
  return MatrixType::NullaryExpr(numRows, numCols, uniform);
}

template <typename MatrixType, typename RandomEngineType>
inline MatrixType StandardNormalRandom(const Eigen::Index numRows, const Eigen::Index numCols,
                                       RandomEngineType& engine) {
  return StandardNormalRandomHelper<MatrixType, typename MatrixType::Scalar, RandomEngineType>::generate(
      numRows, numCols, engine);
}

template <typename MatrixType>
inline void permute_matrix(MatrixType& G, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P,
                           bool bycol = true) {
  if (bycol) {
    P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(G.cols());
    P.setIdentity();
    auto rng = std::default_random_engine{};
    std::shuffle(P.indices().data(), P.indices().data() + P.indices().size(), rng);
    G = G * P;  // permute columns in-place
  } else {
    P = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(G.rows());
    P.setIdentity();
    auto rng = std::default_random_engine{};
    std::shuffle(P.indices().data(), P.indices().data() + P.indices().size(), rng);
    G = P * G;  // permute rows in-place
  }
}

template <typename MatrixType>
void flipOmg(MatrixType& Omg2, MatrixType& Omg) {
  for (Eigen::Index i = 0; i < Omg.cols(); ++i) {
    // if signs of half of values are flipped then correct signs.
    if ((Omg2.col(i) - Omg.col(i)).array().abs().sum() > 2 * (Omg2.col(i) + Omg.col(i)).array().abs().sum()) {
      Omg.col(i) *= -1;
    }
  }
  Omg2 = Omg;
}

// RsvdOp: rows, cols, ranks, oversamples, computeGandH()
template <typename MatrixType>
class RsvdOpOnePass {
 private:
  using Index = Eigen::Index;
  using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

  ConstGenericMatrix mat;
  Index nrow, ncol;
  MatrixType Omg;
  int k, os, size, rand;
  bool trans;  // if matrix is wide then flip the matrix dimension

 public:
  int finder = 1;
  RsvdOpOnePass(ConstGenericMatrix& mat_, int k_, int os_ = 10, int rand_ = 1)
      : mat(mat_), k(k_), os(os_), size(k_ + os_), rand(rand_) {
    if (mat.rows() >= mat.cols()) {
      trans = false;
      nrow = mat.rows();
      ncol = mat.cols();
    } else {
      trans = true;
      nrow = mat.cols();
      ncol = mat.rows();
    }

    auto randomEngine = std::default_random_engine{};
    if (rand == 1) {
      Omg = StandardNormalRandom<MatrixType, std::default_random_engine>(ncol, size, randomEngine);
    } else {
      Omg = UniformRandom<MatrixType, std::default_random_engine>(ncol, size, randomEngine);
    }
  }

  virtual ~RsvdOpOnePass() {}

  Index rows() const { return nrow; }
  Index cols() const { return ncol; }
  Index ranks() const { return k; }
  Index oversamples() const { return os; }

  void computeGandH(MatrixType& G, MatrixType& H, uint32_t p) {
    if (trans) {
      G.noalias() = mat.transpose() * Omg;
      H.noalias() = mat * G;
    } else {
      G.noalias() = mat * Omg;
      H.noalias() = mat.transpose() * G;
    }
    if (p > 0) {
      for (uint32_t pi = 0; pi < p; pi++) {
        if (finder == 1) {
          Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(H);
          Omg.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
        } else if (finder == 2) {
          Eigen::FullPivLU<Eigen::Ref<MatrixType>> lu(H);
          Omg.setIdentity(cols(), size);
          Omg.template triangularView<Eigen::StrictlyLower>() = lu.matrixLU();
        } else {
          throw std::invalid_argument("finder must be 1 or 2");
        }
        if (trans) {
          G.noalias() = mat.transpose() * Omg;
          H.noalias() = mat * G;
        } else {
          G.noalias() = mat * Omg;
          H.noalias() = mat.transpose() * G;
        }
      }
    }
  }

  void computeGandH(MatrixType& G, MatrixType& H, uint32_t p, uint32_t windows) {
    if (windows % 2 != 0) throw std::runtime_error("windows must be a power of 2, ie. windows=2^x.\n");
    if (std::pow(2, p) < windows) throw std::runtime_error("pow(2, p) >= windows has to be met\n");
    uint32_t blocksize = (unsigned int)std::ceil((double)nrow / windows);
    if (blocksize < windows)
      throw std::runtime_error(
          "window size is smaller than number of windows because given matrix "
          "is "
          "too small. please consider other methods or adjust "
          "parameter windows.\n");
    if (trans) {
      G.noalias() = mat.transpose() * Omg;
      H.noalias() = mat * G;
    } else {
      G.noalias() = mat * Omg;
      H.noalias() = mat.transpose() * G;
    }
    uint32_t start_idx, stop_idx, actual_block_size;
    MatrixType H1 = MatrixType::Zero(ncol, size);
    MatrixType H2 = MatrixType::Zero(ncol, size);
    MatrixType Omg2;
    uint32_t band = 1;
    for (uint32_t pi = 0; pi <= p; pi++) {
      if (pi == 0) Omg2 = Omg;
      if (std::pow(2, pi) >= windows) {
        // reset H1, H2 to zero
        H1.setZero();
        H2.setZero();
      }
      band = std::fmin(band * 2, windows);
      for (uint32_t b = 0, i = 1, j = 1; b < windows; ++b, ++i, ++j) {
        start_idx = b * blocksize;
        stop_idx = (b + 1) * blocksize >= nrow ? nrow - 1 : (b + 1) * blocksize - 1;
        actual_block_size = stop_idx - start_idx + 1;
        if (trans)
          G.middleRows(start_idx, actual_block_size).noalias() =
              mat.middleCols(start_idx, actual_block_size).transpose() * Omg;
        else
          G.middleRows(start_idx, actual_block_size).noalias() =
              mat.middleRows(start_idx, actual_block_size) * Omg;
        if (pi > 0 && j <= std::pow(2, pi - 1) && std::pow(2, pi) < windows) {
          if (trans)
            H1.noalias() +=
                mat.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
          else
            H1.noalias() += mat.middleRows(start_idx, actual_block_size).transpose() *
                            G.middleRows(start_idx, actual_block_size);
          // additional complementary power iteration for last read
          if (j == std::pow(2, pi - 1)) {
            H = H1 + H2;
            Eigen::HouseholderQR<MatrixType> qr(H);
            Omg.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
            H2.setZero();
            flipOmg(Omg2, Omg);
          }
        } else if (i <= band / 2) {
          if (trans)
            H1.noalias() +=
                mat.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
          else
            H1.noalias() += mat.middleRows(start_idx, actual_block_size).transpose() *
                            G.middleRows(start_idx, actual_block_size);
        } else if (i > band / 2 && i <= band) {
          if (trans)
            H2.noalias() +=
                mat.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
          else
            H2.noalias() += mat.middleRows(start_idx, actual_block_size).transpose() *
                            G.middleRows(start_idx, actual_block_size);
        }
        if ((b + 1) >= band) {
          if (i == band) {
            H = H1 + H2;
            Eigen::HouseholderQR<MatrixType> qr(H);
            Omg.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
            flipOmg(Omg2, Omg);
            H1.setZero();
            i = 0;
          } else if (i == band / 2) {
            H = H1 + H2;
            Eigen::HouseholderQR<MatrixType> qr(H);
            Omg.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
            flipOmg(Omg2, Omg);
            H2.setZero();
          }
        }
      }
    }
  }
};

template <typename MatrixType, typename RsvdOp>
class RsvdOnePass {
 public:
  RsvdOnePass(RsvdOp& op) : b_op(op) {}

  ~RsvdOnePass() {}

  inline MatrixType singularValues() const { return b_singularValues; }

  inline MatrixType matrixU(bool trans = false) const {
    if (trans) {
      return b_rightSingularVectors;
    } else {
      return b_leftSingularVectors;
    }
  }

  inline MatrixType matrixV(bool trans = false) const {
    if (trans) {
      return b_leftSingularVectors;
    } else {
      return b_rightSingularVectors;
    }
  }

  // G = D * Omega; H = D.transpose() * G;
  void computeUSV(uint32_t p, uint32_t windows = 0) {
    const Eigen::Index nrow{b_op.rows()};
    const Eigen::Index ncol{b_op.cols()};
    const Eigen::Index size{b_op.ranks() + b_op.oversamples()};
    const Eigen::Index k{b_op.ranks()};
    MatrixType H(ncol, size), G(nrow, size), R(size, size), Rt(size, size);

    if (windows > 0)
      b_op.computeGandH(G, H, p, windows);
    else
      b_op.computeGandH(G, H, p);

    {
      Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(G);
      R.noalias() = MatrixType::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
      G.noalias() = qr.householderQ() * MatrixType::Identity(nrow, size);
    }
    {
      Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(G);
      Rt.noalias() = MatrixType::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
      G.noalias() = qr.householderQ() * MatrixType::Identity(nrow, size);
    }

    R = Rt * R;

    // R.T * B = H.T => lapack dtrtrs()
    MatrixType B = R.transpose().fullPivHouseholderQr().solve(H.transpose());
    Eigen::JacobiSVD<MatrixType> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    b_leftSingularVectors.noalias() = G * svd.matrixU().leftCols(k);
    b_rightSingularVectors = svd.matrixV().leftCols(k);
    b_singularValues = svd.singularValues().head(k);
  }

 private:
  RsvdOp& b_op;
  MatrixType b_leftSingularVectors{};
  MatrixType b_singularValues{};
  MatrixType b_rightSingularVectors{};
};

template <typename MatrixType>
class RsvdOne {
 private:
  using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

  ConstGenericMatrix mat;
  uint32_t k, os, rand;
  bool trans;

  RsvdOpOnePass<MatrixType>* op;
  RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>* rsvd;

 public:
  RsvdOne(ConstGenericMatrix& mat_, uint32_t k_, uint32_t os_ = 10, uint32_t rand_ = 1)
      : mat(mat_), k(k_), os(os_), rand(rand_) {
    if (mat.rows() >= mat.cols())
      trans = false;
    else
      trans = true;
    op = new RsvdOpOnePass<MatrixType>(mat, k, os, rand);
    rsvd = new RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>(*op);
  }

  ~RsvdOne() {
    delete op;
    delete rsvd;
  }

  inline void setRangeFinder(int flag) { op->finder = flag; }

  inline void compute(uint32_t p, uint32_t windows = 0) { rsvd->computeUSV(p, windows); }

  inline MatrixType matrixU() const { return rsvd->matrixU(trans); }

  inline MatrixType matrixV() const { return rsvd->matrixV(trans); }

  inline MatrixType singularValues() const { return rsvd->singularValues(); }
};

}  // namespace PCAone

#endif  // PCAONE_RSVD_
