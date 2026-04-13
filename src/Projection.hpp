#ifndef PCAONE_PROJECTION_
#define PCAONE_PROJECTION_

#include "Cmd.hpp"
#include "Data.hpp"

void solve_projection_scores(const Mat2D& V, const ArrBool& C, const Mat2D& G, Mat2D& U);
void run_projection(Data* data, const Param& params);

#endif  // PCAONE_PROJECTION_
