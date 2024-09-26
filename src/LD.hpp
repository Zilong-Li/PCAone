#ifndef PCAONE_LD_
#define PCAONE_LD_

#include "Cmd.hpp"
#include "Data.hpp"
#include "Utils.hpp"

MyArray calc_sds(const MyMatrix& X);

double calc_cor(const MyVector& x, const MyVector& y, const double df);

std::string get_snp_pos_bim(SNPld& snp, const std::string& filebim,
                            bool header = false, Int1D idx = Int1D{0, 3});

std::tuple<Int2D, Int2D> get_target_snp_idx(const SNPld& snp_t,
                                            const SNPld& snp);

void divide_pos_by_window(SNPld& snp, const int ld_window_bp);

void ld_prune_small(Data* data, const std::string& fileout,
                    const std::string& filebim, const SNPld& snp,
                    const double r2_tol);

void ld_prune_big(const std::string& fileout, const std::string& filebim,
                  const MyMatrix& G, const SNPld& snp, double r2_tol);

void ld_clump_single_pheno(std::string fileout, const std::string& head,
                           const int clump_bp, const double clump_r2,
                           const double clump_p1, const double clump_p2,
                           const MyMatrix& G, const Int2D& idx_per_chr,
                           const Int2D& bp_per_chr,
                           const std::vector<UMapIntPds>& pvals_per_chr);

Int1D valid_assoc_file(const std::string& fileassoc,
                       const std::string& colnames);

std::vector<UMapIntPds> map_index_snps(const std::string& fileassoc,
                                       const Int1D& colidx, double clump_p2);

std::vector<UMapIntString> map_assoc_file(const std::string& fileassoc,
                                          const Int1D& colidx);

void ld_r2_big(const MyMatrix& G, const SNPld& snp, const std::string& filebim,
               const std::string& fileout);

void ld_r2_small(Data* data, const SNPld& snp, const std::string& filebim,
                 const std::string& fileout);

void run_ld_stuff(const Param& params, Data* data);

#endif  // PCAONE_LD_
