#include "SVD.hpp"

void SvdOp::perform_op(const double *x_in, double* y_out)
{
   Map<const VectorXd> x(x_in, n);
   Map<VectorXd> y(y_out, n);
   // initParallel();
   // setNbThreads(20);
   y.noalias() = mat * (mat.transpose() * x);
   // y.noalias() = mat.transpose() * (mat * x);
}