#ifndef __SVD_BLAS_H__
#define __SVD_BLAS_H__

#include <cblas.h>
#include <Eigen/Core>

using namespace Eigen;

typedef unsigned long long uint64;

// extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);


class SvdBlas
{
public:
    SvdBlas(const MatrixXd& mat_) : mat(mat_), m(mat.rows()), n(mat.cols()), dim(n), workm(m)
    {}

    ~SvdBlas() {}

    void perform_op(const double *x_in, double* y_out) {
        perform_blas_matprod(x_in, workm.data());
        perform_blas_mattprod(workm.data(), y_out);
    }

    // y_out = A * x_in
    void perform_blas_matprod(const double *x_in, double* y_out) {
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, mat.data(), m, x_in, 1, 0.0, y_out, 1);
    }
    // y_out = A' * x_in
    void perform_blas_mattprod(const double *x_in, double* y_out) {
        cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, mat.data(), m, x_in, 1, 0.0, y_out, 1);
    }

    // Number of rows and columns of the operator B'B, not of B itself
    inline unsigned int rows() const { return dim; }
    inline unsigned int cols() const { return dim; }

private:
    const MatrixXd& mat;
    const uint64 n;
    const uint64 m;
    const uint64 dim;
    VectorXd workm;
    VectorXd workn;
};

class SvdBlasWide
{
public:
    SvdBlasWide(MatrixXd& mat_) : mat(mat_), m(mat.rows()), n(mat.cols()), dim(m), workn(n)
    {}

    ~SvdBlasWide() {}

    void perform_op(const double *x_in, double* y_out) {
        Map<const VectorXd> x(x_in, dim);
        Map<VectorXd> y(y_out, dim);
        perform_blas_mattprod(x.data(), workn.data());
        perform_blas_matprod(workn.data(), y.data());
    }

    // y_out = A * x_in
    void perform_blas_matprod(const double *x_in, double* y_out) {
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, mat.data(), m, x_in, 1, 0.0, y_out, 1);
    }
    // y_out = A' * x_in
    void perform_blas_mattprod(const double *x_in, double* y_out) {
        cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, mat.data(), m, x_in, 1, 0.0, y_out, 1);
    }

    // Number of rows and columns of the operator BB', not of B itself
    inline unsigned int rows() const { return dim; }
    inline unsigned int cols() const { return dim; }

private:
    const MatrixXd& mat;
    const uint64 n;
    const uint64 m;
    const uint64 dim;
    VectorXd workm;
    VectorXd workn;
};

#endif
