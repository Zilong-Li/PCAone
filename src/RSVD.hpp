#ifndef __RSVD__
#define __RSVD__

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>

template <typename MatrixType>
using RealType = typename Eigen::NumTraits<typename MatrixType::Scalar>::Real;

template <typename MatrixType, typename ScalarType, typename RandomEngineType>
struct StandardNormalRandomHelper {
  static inline MatrixType generate(Eigen::Index numRows, Eigen::Index numCols,
                                    RandomEngineType &engine);
};

template <typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, RealType<MatrixType>, RandomEngineType> {
  static inline MatrixType generate(const Eigen::Index numRows, const Eigen::Index numCols,
                                    RandomEngineType &engine) {
    std::normal_distribution<RealType<MatrixType>> distribution{0, 1};
    const auto normal{[&](typename MatrixType::Scalar) { return distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, normal);
  }
};

template <typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, std::complex<RealType<MatrixType>>,
                                  RandomEngineType> {
  static inline MatrixType generate(const Eigen::Index numRows, const Eigen::Index numCols,
                                    RandomEngineType &engine) {
    constexpr RealType<MatrixType> stdDev{0.707106781186547};
    std::normal_distribution<RealType<MatrixType>> distribution{0, stdDev};
    const auto complexNormal{[&](typename MatrixType::Scalar) {
      return std::complex<RealType<MatrixType>>(distribution(engine), distribution(engine));
    }};
    return MatrixType::NullaryExpr(numRows, numCols, complexNormal);
  }
};

// @see: https://github.com/mp4096/rsvd/blob/master/include/rsvd/StandardNormalRandom.hpp
template <typename MatrixType, typename RandomEngineType>
inline MatrixType StandardNormalRandom(const Eigen::Index numRows, const Eigen::Index numCols, RandomEngineType &engine) {
  return StandardNormalRandomHelper<MatrixType, typename MatrixType::Scalar,
                                    RandomEngineType>::generate(numRows, numCols, engine);
};


// RsvdOp: rows, cols, ranks, oversamples, computeGandH()
template <typename MatrixType>
class RsvdOpOnePass
{
private:
    using Index = Eigen::Index;

public:

    virtual Index rows() const = 0;
    virtual Index cols() const = 0;
    virtual Index ranks() const = 0;
    virtual Index oversamples() const = 0;

    virtual void computeGandH(MatrixType& G, MatrixType& H, int p) = 0;

    virtual ~RsvdOpOnePass() {}
};


template <typename MatrixType>
class RsvdOpTallMat : public RsvdOpOnePass<MatrixType>
{
private:
    using Index = Eigen::Index;
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix mat;
    const Index nrow, ncol, nk, os, size;
    MatrixType Omg;

public:
    RsvdOpTallMat(ConstGenericMatrix& mat_, int k_, int os_ = 10) :
        mat(mat_), nrow(mat_.rows()), ncol(mat_.cols()), nk(k_), os(os_), size(k_ + os_)
        {
            std::mt19937_64 randomEngine{};
            randomEngine.seed(111);
            Omg = StandardNormalRandom<MatrixType, std::mt19937_64>(ncol, size, randomEngine);
        }

    ~RsvdOpTallMat() {}

    Index rows() const { return nrow; }
    Index cols() const { return ncol; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MatrixType& G, MatrixType& H, int p)
        {
            G.noalias() = mat * Omg;
            H.noalias() = mat.transpose() * G;
            if (p > 0)
            {
                for(int i=0; i<p; ++i)
                {
                    // Eigen::ColPivHouseholderQR<Eigen::Ref<MatrixType>> qr(H);
                    Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(H);
                    H.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
                    G.noalias() = mat * H;
                    H.noalias() = mat.transpose() * G;
                }
            }
        }

};

template <typename MatrixType>
class RsvdOpWideMat : public RsvdOpOnePass<MatrixType>
{
private:
    using Index = Eigen::Index;
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix mat;
    const Index nrow, ncol, nk, os, size;
    MatrixType Omg;

public:
    RsvdOpWideMat(ConstGenericMatrix& mat_, int k_, int os_ = 10) :
        mat(mat_), nrow(mat_.cols()), ncol(mat_.rows()), nk(k_), os(os_), size(k_ + os_)
        {
            std::mt19937_64 randomEngine{};
            randomEngine.seed(111);
            Omg = StandardNormalRandom<MatrixType, std::mt19937_64>(ncol, size, randomEngine);
        }

    ~RsvdOpWideMat() {}

    Index rows() const { return nrow; }
    Index cols() const { return ncol; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MatrixType& G, MatrixType& H, int p)
        {
            G.noalias() = mat.transpose() * Omg;
            H.noalias() = mat * G;
            if (p > 0)
            {
                for(int i=0; i<p; ++i)
                {
                    // Eigen::ColPivHouseholderQR<Eigen::Ref<MatrixType>> qr(H);
                    Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(H);
                    H.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
                    G.noalias() = mat.transpose() * H;
                    H.noalias() = mat * G;
                }
            }
        }

};


template <typename MatrixType, typename RsvdOp>
class RsvdOnePass
{
public:

    RsvdOnePass(RsvdOp& op) : b_op(op) {}

    ~RsvdOnePass() {}

    inline MatrixType singularValues() const { return b_singularValues; }

    inline MatrixType matrixU() const { return b_leftSingularVectors; }

    inline MatrixType matrixV() const { return b_rightSingularVectors; }

    // G = D * Omega; H = D.transpose() * G;
    void computeUSV(int p = 1)
        {
            const Eigen::Index nrow{b_op.rows()};
            const Eigen::Index ncol{b_op.cols()};
            const Eigen::Index size{b_op.ranks() + b_op.oversamples()};
            const Eigen::Index k{b_op.ranks()};
            MatrixType G(nrow, size), H(ncol, size), Q(nrow, size), B(size, ncol), R(size, size), Rt(size, size);
            b_op.computeGandH(G, H, p);

            Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(G);
            Q.noalias() = qr.householderQ() * MatrixType::Identity(nrow, size);
            R.noalias() = MatrixType::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
            Eigen::HouseholderQR<MatrixType> qr2(Q);
            Q.noalias() = qr2.householderQ() * MatrixType::Identity(nrow, size);
            Rt.noalias() = MatrixType::Identity(size, nrow) * qr2.matrixQR().template triangularView<Eigen::Upper>();
            R = Rt * R;
            // R.T * B = H.T
            B.noalias() = R.transpose().householderQr().solve(H.transpose());
            Eigen::JacobiSVD<MatrixType> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
            b_leftSingularVectors.noalias() = Q * svd.matrixU().leftCols(k);
            b_rightSingularVectors = svd.matrixV().leftCols(k);
            b_singularValues = svd.singularValues().head(k);
        }

    // for fancy in future
    void computeUSV(MatrixType& G, MatrixType& H, int p = 1)
        {
            b_op.computeGandH(G, H, p);

            const Eigen::Index nrow{G.rows()};
            const Eigen::Index ncol{H.cols()};
            const Eigen::Index size{b_op.ranks() + b_op.oversamples()};
            const Eigen::Index k{b_op.ranks()};
            MatrixType Q(nrow, size), B(size, ncol), R(size, size), Rt(size, size);

            Eigen::HouseholderQR<MatrixType> qr(G);  // shouldn't use in place qr
            Q.noalias() = qr.householderQ() * MatrixType::Identity(nrow, size);
            R.noalias() = MatrixType::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
            Eigen::HouseholderQR<MatrixType> qr2(Q);
            Q.noalias() = qr2.householderQ() * MatrixType::Identity(nrow, size);
            Rt.noalias() = MatrixType::Identity(size, nrow) * qr2.matrixQR().template triangularView<Eigen::Upper>();
            R = Rt * R;

            B.noalias() = R.transpose().householderQr().solve(H.transpose());
            Eigen::JacobiSVD<MatrixType> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
            b_leftSingularVectors.noalias() = Q * svd.matrixU().leftCols(k);
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
class Rsvd
{
private:
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix mat;
    int  k, p, os;
    bool flip = false;

    RsvdOpOnePass<MatrixType>* op;
    RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>* rsvd;

public:
    Rsvd(ConstGenericMatrix& mat_, int k_, int p_, int os_ = 10) :
        mat(mat_), k(k_), p(p_), os(os_)
        {
            if (mat.rows() >= mat.cols())
            {
                op = new RsvdOpTallMat<MatrixType>(mat, k, os);
            } else {
                op = new RsvdOpWideMat<MatrixType>(mat, k, os);
                flip = true;
            }
            rsvd = new RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>(*op);
            rsvd->computeUSV(p);
        }

    ~Rsvd()
        {
            delete op;
            delete rsvd;
        }

    inline MatrixType matrixU() const { if (flip) { return rsvd->matrixV();} else {return rsvd->matrixU();} }

    inline MatrixType matrixV() const { if (flip) { return rsvd->matrixU();} else {return rsvd->matrixV();} }

    inline MatrixType singularValues() const { return rsvd->singularValues(); }
};


#endif
