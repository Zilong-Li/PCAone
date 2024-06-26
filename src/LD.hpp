#ifndef PCAONE_LD_
#define PCAONE_LD_

#include "Cmd.hpp"
#include "Data.hpp"
#include "Utils.hpp"

MyVector calc_sds(const MyMatrix& X);

std::string get_snp_pos_bim(const std::string& filebim, Int1D& pos,
                            Int1D& chr_pos_end, std::vector<std::string>& chrs,
                            bool header = false, Int1D idx = Int1D{0, 3});

std::tuple<Int2D, Int2D, std::string> get_target_snp_idx(
    const std::string& filebim, const Int1D& pos, const Int1D& chr_pos_end,
    const std::vector<std::string>& chrs, bool header = false,
    Int1D colidx = Int1D{0, 3});

void calc_ld_metrics(const std::string& fileout, const std::string& filebim,
                     const MyMatrix& G, const MyVector& F, const Int1D& snp_pos,
                     const Int1D& chr_pos_end, int ld_window_bp, double r2_tol,
                     bool verbose);

void calc_ld_clump(std::string fileout, std::string fileassoc,
                   std::string colnames, int clump_bp, double clump_r2,
                   double clump_p1, double clump_p2, const MyMatrix& G,
                   const Int1D& snp_pos, const Int1D& chr_pos_end,
                   const std::vector<std::string>& chrs);

Int1D valid_assoc_file(const std::string& fileassoc,
                       const std::string& colnames);

inline UMapIntInt vector2map(const Int1D& v) {
  UMapIntInt m;
  int i = 0;
  for (auto k : v) m.insert({k, i++});
  return m;
}

std::vector<UMapIntPds> map_index_snps(const std::string& fileassoc,
                                       const Int1D& colidx, double clump_p2);

std::vector<UMapIntString> map_assoc_file(const std::string& fileassoc,
                                          const Int1D& colidx);

void run_ld_pruning(const Param& params, const MyMatrix& G, const MyVector& F);

void run_ld_stuff(const Param& params, Data* data);

#endif  // PCAONE_LD_
