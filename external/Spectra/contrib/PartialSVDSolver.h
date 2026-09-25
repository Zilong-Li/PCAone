// Copyright (C) 2018-2022 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SPECTRA_PARTIAL_SVD_SOLVER_H
#define SPECTRA_PARTIAL_SVD_SOLVER_H

#include <Eigen/Core>
#include "../SymEigsSolver.h"

namespace Spectra {

// Abstract class for matrix operation
template <typename Scalar_>
class SVDMatOp
{
public:
    using Scalar = Scalar_;

private:
    using Index = Eigen::Index;

public:
    virtual Index rows() const = 0;
    virtual Index cols() const = 0;

    // y_out = A' * A * x_in or y_out = A * A' * x_in
    virtual void perform_op(const Scalar* x_in, Scalar* y_out) const = 0;

    virtual ~SVDMatOp() {}
};

// Operation of a tall matrix in SVD
// We compute the eigenvalues of A' * A
// MatrixType is either Eigen::Matrix<Scalar, ...> or Eigen::SparseMatrix<Scalar, ...>
template <typename Scalar, typename MatrixType>
class SVDTallMatOp : public SVDMatOp<Scalar>
{
private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix m_mat;
    const Index m_dim;
    mutable Vector m_cache;

public:
    // Constructor
    SVDTallMatOp(ConstGenericMatrix& mat) :
        m_mat(mat),
        m_dim((std::min)(mat.rows(), mat.cols())),
        m_cache(mat.rows())
    {}

    // These are the rows and columns of A' * A
    Index rows() const override { return m_dim; }
    Index cols() const override { return m_dim; }

    // y_out = A' * A * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const override
    {
        MapConstVec x(x_in, m_mat.cols());
        MapVec y(y_out, m_mat.cols());
        m_cache.noalias() = m_mat * x;
        y.noalias() = m_mat.transpose() * m_cache;
    }
};

// Operation of a wide matrix in SVD
// We compute the eigenvalues of A * A'
// MatrixType is either Eigen::Matrix<Scalar, ...> or Eigen::SparseMatrix<Scalar, ...>
template <typename Scalar, typename MatrixType>
class SVDWideMatOp : public SVDMatOp<Scalar>
{
private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix m_mat;
    const Index m_dim;
    mutable Vector m_cache;

public:
    // Constructor
    SVDWideMatOp(ConstGenericMatrix& mat) :
        m_mat(mat),
        m_dim((std::min)(mat.rows(), mat.cols())),
        m_cache(mat.cols())
    {}

    // These are the rows and columns of A * A'
    Index rows() const override { return m_dim; }
    Index cols() const override { return m_dim; }

    // y_out = A * A' * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const override
    {
        MapConstVec x(x_in, m_mat.rows());
        MapVec y(y_out, m_mat.rows());
        m_cache.noalias() = m_mat.transpose() * x;
        y.noalias() = m_mat * m_cache;
    }
};

// Partial SVD solver
// MatrixType is either Eigen::Matrix<Scalar, ...> or Eigen::SparseMatrix<Scalar, ...>
template <typename MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
class PartialSVDSolver
{
private:
    using Scalar = typename MatrixType::Scalar;
    using Index = Eigen::Index;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

    ConstGenericMatrix m_mat;
    const Index m_m;
    const Index m_n;
    SVDMatOp<Scalar>* m_op;
    SymEigsSolver<SVDMatOp<Scalar>>* m_eigs;
    Index m_nconv;
    Matrix m_evecs;

public:
    // Constructor
    PartialSVDSolver(ConstGenericMatrix& mat, Index ncomp, Index ncv) :
        m_mat(mat), m_m(mat.rows()), m_n(mat.cols()), m_evecs(0, 0)
    {
        // Determine the matrix type, tall or wide
        if (m_m > m_n)
        {
            m_op = new SVDTallMatOp<Scalar, MatrixType>(mat);
        }
        else
        {
            m_op = new SVDWideMatOp<Scalar, MatrixType>(mat);
        }

        // Solver object
        m_eigs = new SymEigsSolver<SVDMatOp<Scalar>>(*m_op, ncomp, ncv);
    }

    // Destructor
    virtual ~PartialSVDSolver()
    {
        delete m_eigs;
        delete m_op;
    }

    // Computation
    Index compute(Index maxit = 1000, Scalar tol = 1e-10)
    {
        // PCAone patch: matrix_U()/matrix_V() cache the eigenvectors in
        // m_evecs and only fetch them when the cache is empty. A solver whose
        // matrix is updated in place and re-computed -- which is exactly what
        // the EM-PCA loop in Arnoldi.cpp does -- then keeps handing back the
        // singular vectors of the FIRST decomposition forever. Invalidate the
        // cache here so every compute() is observed by the accessors.
        m_evecs.resize(0, 0);
        m_eigs->init();
        m_nconv = m_eigs->compute(SortRule::LargestAlge, maxit, tol);

        return m_nconv;
    }

    // The converged singular values
    Vector singular_values() const
    {
        Vector svals = m_eigs->eigenvalues().cwiseSqrt();

        return svals;
    }

    // The converged left singular vectors
    Matrix matrix_U(Index nu)
    {
        if (m_evecs.cols() < 1)
        {
            m_evecs = m_eigs->eigenvectors();
        }
        nu = (std::min)(nu, m_nconv);
        if (m_m <= m_n)
        {
            return m_evecs.leftCols(nu);
        }

        return m_mat * (m_evecs.leftCols(nu).array().rowwise() / m_eigs->eigenvalues().head(nu).transpose().array().sqrt()).matrix();
    }

    // The converged right singular vectors
    Matrix matrix_V(Index nv)
    {
        if (m_evecs.cols() < 1)
        {
            m_evecs = m_eigs->eigenvectors();
        }
        nv = (std::min)(nv, m_nconv);
        if (m_m > m_n)
        {
            return m_evecs.leftCols(nv);
        }

        // PCAone patch: the product in column panels of m_mat, one thread each.
        // A multithreaded m_mat' * W packed a kc x cols() block of m_mat' --
        // nearly a copy of m_mat when it has few rows -- and was slower.
        const Matrix W = (m_evecs.leftCols(nv).array().rowwise() / m_eigs->eigenvalues().head(nv).transpose().array().sqrt()).matrix();
        Matrix V(m_n, nv);
        const Index bs = (std::max)(Index(64), (Index(1) << 18) / (std::max)(Index(1), m_m));
        const Index nb = (m_n + bs - 1) / bs;
#pragma omp parallel for schedule(static)
        for (Index b = 0; b < nb; ++b)
        {
            const Index c = b * bs, w = (std::min)(bs, m_n - c);
            V.middleRows(c, w).noalias() = m_mat.middleCols(c, w).transpose() * W;
        }
        return V;
    }
};

}  // namespace Spectra

#endif  // SPECTRA_PARTIAL_SVD_SOLVER_H
