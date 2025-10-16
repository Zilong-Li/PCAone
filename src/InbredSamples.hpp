#ifndef PCAONE_INBREDSAMPLES_
#define PCAONE_INBREDSAMPLES_

#include "Cmd.hpp"
#include "Data.hpp"

// PI is N x M
// GL is (N x 2) x M
// D is Fnew-F
// F is updated as Fnew in-place
void inbreed_coef_sample(Mat1D& D, Mat1D& F, const Mat2D& PI, const Mat2D& GL, const int type,
                             const int size, const uint start);

// out of core version
void inbreed_coef_sample_ooc(Mat1D& D1, Mat1D& F, Data* data, Data* Pi, const Param& params);

void inbreed_coef_sample_em(int type, const Mat2D& GL, const Mat2D& PI, const Param& params);

void run_inbreed_coef_sample(Data* Pi, const Param& params);

void write_coef_per_sample(const std::string& fout, const std::string& fam, const Mat1D& coef);

#endif  // PCAONE_INBREDSAMPLES_
