#ifndef PCAONE_RSVD_
#define PCAONE_RSVD_

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>
#include <stdexcept>

template <typename MatrixType>
using RealType = typename Eigen::NumTraits<typename MatrixType::Scalar>::Real;

template <typename MatrixType, typename ScalarType, typename RandomEngineType>
struct StandardNormalRandomHelper
{
    static inline MatrixType generate(Eigen::Index numRows, Eigen::Index numCols, RandomEngineType& engine);
};

template <typename MatrixType, typename RandomEngineType>
inline MatrixType UniformRandom(const Eigen::Index numRows, const Eigen::Index numCols, RandomEngineType& engine)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{-1, 1};
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
};

template <typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, RealType<MatrixType>, RandomEngineType>
{
    static inline MatrixType generate(const Eigen::Index numRows, const Eigen::Index numCols, RandomEngineType& engine)
    {
        std::normal_distribution<RealType<MatrixType>> distribution{0, 1};
        const auto normal{[&](typename MatrixType::Scalar) { return distribution(engine); }};
        return MatrixType::NullaryExpr(numRows, numCols, normal);
    }
};

template <typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, std::complex<RealType<MatrixType>>, RandomEngineType>
{
    static inline MatrixType generate(const Eigen::Index numRows, const Eigen::Index numCols, RandomEngineType& engine)
    {
        constexpr RealType<MatrixType> stdDev{0.707106781186547};
        std::normal_distribution<RealType<MatrixType>> distribution{0, stdDev};
        const auto complexNormal{[&](typename MatrixType::Scalar) { return std::complex<RealType<MatrixType>>(distribution(engine), distribution(engine)); }};
        return MatrixType::NullaryExpr(numRows, numCols, complexNormal);
    }
};

template <typename MatrixType, typename RandomEngineType>
inline MatrixType StandardNormalRandom(const Eigen::Index numRows, const Eigen::Index numCols, RandomEngineType& engine)
{
    return StandardNormalRandomHelper<MatrixType, typename MatrixType::Scalar, RandomEngineType>::generate(numRows, numCols, engine);
};


// RsvdOp: rows, cols, ranks, oversamples, computeGandH()
template <typename MatrixType>
class RsvdOpOnePass
{
private:
    using Index = Eigen::Index;
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix mat;
    Index nrow, ncol;
    MatrixType Omg;
    int k, os, size, rand;
    bool trans; // if matrix is wide then flip the matrix dimension

public:
    RsvdOpOnePass(ConstGenericMatrix& mat_, int k_, int os_ = 10, int rand_ = 1) : mat(mat_), k(k_), os(os_), size(k_ + os_), rand(rand_)
    {
        if (mat.rows() >= mat.cols())
        {
            trans = false;
            nrow = mat.rows();
            ncol = mat.cols();
        }
        else
        {
            trans = true;
            nrow = mat.cols();
            ncol = mat.rows();
        }

        auto randomEngine = std::default_random_engine{};
        if (rand == 1)
        {
            Omg = StandardNormalRandom<MatrixType, std::default_random_engine>(ncol, size, randomEngine);
        }
        else
        {
            Omg = UniformRandom<MatrixType, std::default_random_engine>(ncol, size, randomEngine);
        }
    }

    virtual ~RsvdOpOnePass()
    {
    }

    Index rows() const
    {
        return nrow;
    }
    Index cols() const
    {
        return ncol;
    }
    Index ranks() const
    {
        return k;
    }
    Index oversamples() const
    {
        return os;
    }

    void computeGandH(MatrixType& G, MatrixType& H, int p)
    {
        if (trans)
        {
            G.noalias() = mat.transpose() * Omg;
            H.noalias() = mat * G;
        }
        else
        {
            G.noalias() = mat * Omg;
            H.noalias() = mat.transpose() * G;
        }
        if (p > 0)
        {
            for (int i = 0; i < p; i++)
            {
                Eigen::HouseholderQR<Eigen::Ref<MatrixType>> qr(H);
                H.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
                if (trans)
                {
                    G.noalias() = mat.transpose() * Omg;
                    H.noalias() = mat * G;
                }
                else
                {
                    G.noalias() = mat * Omg;
                    H.noalias() = mat.transpose() * G;
                }
            }
        }
    }

    void computeGandH(MatrixType& G, MatrixType& H, int p, int windows)
    {
        if (windows % 2 != 0)
            throw std::runtime_error("windows must be a power of 2, ie. windows=2^x.\n");
        uint blocksize = (unsigned int)std::ceil((double)nrow / windows);
        if (blocksize < windows)
            throw std::runtime_error("window size is smaller than number of windows because given matrix is too small. please consider other methods or adjust "
                                     "parameter windows.\n");
        uint start_idx, stop_idx, actual_block_size;
        MatrixType H1, H2;
        uint band = 2;
        for (int pi = 0; pi <= p; pi++)
        {
            band = std::fmin(band * 2, windows);
            H1 = MatrixType::Zero(ncol, size);
            H2 = MatrixType::Zero(ncol, size);
            for (uint b = 0, i = 1; b < windows; ++b, ++i)
            {
                start_idx = b * blocksize;
                stop_idx = (b + 1) * blocksize >= nrow ? nrow - 1 : (b + 1) * blocksize - 1;
                actual_block_size = stop_idx - start_idx + 1;
                if (trans)
                    G.middleRows(start_idx, actual_block_size).noalias() = mat.middleCols(start_idx, actual_block_size).transpose() * Omg;
                else
                    G.middleRows(start_idx, actual_block_size).noalias() = mat.middleRows(start_idx, actual_block_size) * Omg;
                if (i <= band / 2)
                {
                    if (trans)
                        H1.noalias() += mat.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                    else
                        H1.noalias() += mat.middleRows(start_idx, actual_block_size).transpose() * G.middleRows(start_idx, actual_block_size);
                }
                else if (i > band / 2 && i <= band)
                {
                    if (trans)
                        H2.noalias() += mat.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                    else
                        H2.noalias() += mat.middleRows(start_idx, actual_block_size).transpose() * G.middleRows(start_idx, actual_block_size);
                }
                if ((b + 1) >= band)
                {
                    if (i == band)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixType> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
                        H1 = MatrixType::Zero(cols(), size);
                        i = 0;
                    }
                    else if (i == band / 2)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixType> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixType::Identity(cols(), size);
                        H2 = MatrixType::Zero(cols(), size);
                    }
                }
            }
        }
    }
};

template <typename MatrixType, typename RsvdOp>
class RsvdOnePass
{
public:
    RsvdOnePass(RsvdOp& op) : b_op(op)
    {
    }

    ~RsvdOnePass()
    {
    }

    inline MatrixType singularValues() const
    {
        return b_singularValues;
    }

    inline MatrixType matrixU(bool trans = false) const
    {
        if (trans)
        {
            return b_rightSingularVectors;
        }
        else
        {
            return b_leftSingularVectors;
        }
    }

    inline MatrixType matrixV(bool trans = false) const
    {
        if (trans)
        {
            return b_leftSingularVectors;
        }
        else
        {
            return b_rightSingularVectors;
        }
    }

    // G = D * Omega; H = D.transpose() * G;
    void computeUSV(int p, int windows = 0)
    {
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
        MatrixType B = R.transpose().colPivHouseholderQr().solve(H.transpose());
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
class RsvdOne
{
private:
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix mat;
    uint k, os, rand;
    bool trans;

    RsvdOpOnePass<MatrixType>* op;
    RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>* rsvd;

public:
    RsvdOne(ConstGenericMatrix& mat_, uint k_, uint os_ = 10, uint rand_ = 1) : mat(mat_), k(k_), os(os_), rand(rand_)
    {
        if (mat.rows() >= mat.cols())
            trans = false;
        else
            trans = true;
        op = new RsvdOpOnePass<MatrixType>(mat, k, os, rand);
        rsvd = new RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>(*op);
    }

    ~RsvdOne()
    {
        delete op;
        delete rsvd;
    }

    inline void compute(int p, int windows = 0)
    {
        rsvd->computeUSV(p, windows);
    }

    inline MatrixType matrixU() const
    {
        return rsvd->matrixU(trans);
    }

    inline MatrixType matrixV() const
    {
        return rsvd->matrixV(trans);
    }

    inline MatrixType singularValues() const
    {
        return rsvd->singularValues();
    }
};

#endif // PCAONE_RSVD_
