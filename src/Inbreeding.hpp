#ifndef PCAONE_INBREEDING_
#define PCAONE_INBREEDING_

#include "Cmd.hpp"
#include "Data.hpp"

void calc_inbreed_coef(Mat2D& Fd, const Mat1D& F, const Mat2D& PI, const Mat2D& GL);

void run_inbreeding_em(Mat1D& F, const Mat2D& PI, const Mat2D& GL, const Param& params);

void run_inbreeding(Data* data, const Param& params);

#endif  // PCAONE_INBREEDING_
