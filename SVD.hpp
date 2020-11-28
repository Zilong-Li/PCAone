#ifndef __SVD_H__
#define __SVD_H__

#include <Eigen/Core>

using namespace Eigen;

typedef unsigned long long uint64;

class SvdOp
{
public:
    // this is for tall matrix, i.e mat.rows() > mat.cols()
    SvdOp(const MatrixXd& mat_): mat(mat_), n(mat_.rows()) {}

    ~SvdOp() {}

    inline unsigned int rows() const { return n; }
    inline unsigned int cols() const { return n; }

    // we calculate y_out = A' * A * x_in
    void perform_op(const double *x_in, double* y_out);

private:
    const MatrixXd& mat;
    const uint64 n;
};


#endif
