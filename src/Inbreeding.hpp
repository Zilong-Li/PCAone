#ifndef PCAONE_INBREEDING_
#define PCAONE_INBREEDING_

#include "Cmd.hpp"
#include "Data.hpp"

// PI is N x M
// GL is (N x 2) x M
// D is Fnew-F
// F is updated as Fnew in-place
void calc_inbreed_coef(Mat1D& D, Mat1D& F, const Mat2D& PI, const Mat2D& GL, const int type, const int size,
                       const uint start);

// out of core version
void calc_inbreed_coef_outofcore(Mat1D& D1, Mat1D& F, Data* data, Data* Pi, const Param& params);

// log-likelihoods test per site
// F is the estimated inbreeding coef
// T is output
void calc_inbreed_site_lrt(Mat1D& T, const Mat1D& F, const Mat2D& PI, const Mat2D& GL, const int type,
                           const int size, const uint start);

void write_hwe_per_site(const std::string& fout, const std::string& fbim, const Mat1D& hwe, const Mat1D& lrt,
                        const Mat1D& coef);

void run_inbreeding_em(int type, const Mat2D& GL, const Mat2D& PI, const Param& params);

void run_inbreeding(Data* Pi, const Param& params);

#endif  // PCAONE_INBREEDING_
