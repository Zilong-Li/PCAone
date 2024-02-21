#ifndef PCAONE_LD_
#define PCAONE_LD_

#include "Utils.hpp"

MyVector calc_sds(const MyMatrix & X);


void calc_ld_metrics(std::string fileout,
                     const MyMatrix & G,
                     const MyVector & F,
                     const Int1D & snp_pos,
                     const Int1D & chr_pos_end,
                     int ld_window_bp,
                     double r2_tol,
                     bool verbose);

void calc_ld_pairs(std::string fileout,
                   std::string filebim,
                   const MyMatrix & G,
                   const MyVector & F,
                   const Int1D & snp_pos,
                   const Int1D & chr_pos_end,
                   const std::vector<std::string> & chrs);

void calc_ld_clump(std::string fileout,
                   std::string fileassoc,
                   int clump_bp,
                   double clump_r2,
                   double clump_p1,
                   double clump_p2,
                   const MyMatrix & G,
                   const MyVector & F,
                   const Int1D & snp_pos,
                   const Int1D & chr_pos_end,
                   const std::vector<std::string> & chrs);

Int1D valid_assoc_file(const std::string & fileassoc);



#endif // PCAONE_LD_
