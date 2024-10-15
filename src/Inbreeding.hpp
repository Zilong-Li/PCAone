#ifndef PCAONE_INBREEDING_
#define PCAONE_INBREEDING_

#include "Cmd.hpp"
#include "Data.hpp"

// PI is N x M
// GL is (N x 2) x M
// D is Fnew-F
// F is updated as Fnew in-place
void calc_inbreed_coef(Mat1D& D, Mat1D& F, const Mat2D& PI, const Mat2D& GL);

// log-likelihoods test per site
// F is the estimated inbreeding coef
// T is output
void calc_inbreed_site_lrt(Mat1D& T, const Mat1D& F, const Mat2D& PI, const Mat2D& GL);

void write_hwe_per_site(const std::string& fout, const std::string& fbim, const Mat1D& hwe, const Mat1D& lrt,
                        const Mat1D& coef);

void run_inbreeding_em(Mat1D& F, const Mat2D& PI, const Mat2D& GL, const Param& params);

void run_inbreeding(Data* data, const Param& params);

#endif  // PCAONE_INBREEDING_
