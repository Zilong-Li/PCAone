#ifndef PCAONE_PROJECTION_
#define PCAONE_PROJECTION_

#include "Cmd.hpp"
#include "Data.hpp"

Mat2D solve_bootstrap_projection_no_missing(const Mat2D& design, const Mat2D& G, const std::vector<uint>& counts);
Mat2D solve_bootstrap_projection_missing(const Mat2D& design,
                                         const ArrBool& C,
                                         const Mat2D& G,
                                         const std::vector<uint>& counts);
void write_projection_bootstrap_stats(
    const Mat2D& design, const ArrBool& C, const Mat2D& G, const Mat2D& U0, uint nreps, const Param& params);

void solve_projection_scores(const Mat2D& V, const ArrBool& C, const Mat2D& G, Mat2D& U);
void run_projection(Data* data, const Param& params);

#endif  // PCAONE_PROJECTION_
