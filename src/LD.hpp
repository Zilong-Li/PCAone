#ifndef PCAONE_LD_
#define PCAONE_LD_

#include "Cmd.hpp"
#include "Data.hpp"
#include "Utils.hpp"

struct SNPld {
  std::vector<int> pos;          // pos of each SNP
  std::vector<int> end_pos;      // 0-based index for last snp pos
  std::vector<std::string> chr;  // chr sequences
  Double1D af;                   // allele frequency
  std::vector<int> ws;           //  the snp index, i.e the index for lead SNP
  std::vector<int> we;           // the number of SNPs (including lead SNP) in a window
};

// get the chr, bp, id
struct BIM {
  std::string data;
  operator std::string const &() const { return data; }
  friend std::istream& operator>>(std::istream& is, BIM& BIM) {
    std::getline(is, BIM.data);
    if (BIM.data.empty()) return is;
    auto token = split_string(BIM.data, " \t");
    BIM.data = token[0] + "\t" + token[3] + "\t" + token[1];
    return is;
  }
};

Arr1D calc_sds(const Mat2D& X);

double calc_cor(const Mat1D& x, const Mat1D& y, const double df);

double calc_cor_correct_n(const Mat1D& x, const Mat1D& y);

std::string get_snp_pos_bim(SNPld& snp, const std::string& filebim, bool header = false,
                            Int1D idx = Int1D{0, 3});

std::tuple<Int2D, Int2D> get_target_snp_idx(const SNPld& snp_t, const SNPld& snp);

void divide_pos_by_window(SNPld& snp, const int ld_window_bp);

void ld_prune_small(Data* data, const std::string& fileout, const std::string& filebim, const SNPld& snp,
                    const double r2_tol);

void ld_prune_big(const Mat2D& G, const SNPld& snp, double r2_tol, const std::string& fileout,
                  const std::string& filebim, const bool r2correction);

void ld_clump_single_pheno(const std::string& fileout, const std::string& head, const int clump_bp,
                           const double clump_r2, const double clump_p1, const double clump_p2,
                           const Mat2D& G, const Int2D& idx_per_chr, const Int2D& bp_per_chr,
                           const std::vector<UMapIntPds>& pvals_per_chr);

Int1D valid_assoc_file(const std::string& fileassoc, const std::string& colnames);

std::vector<UMapIntPds> map_index_snps(const std::string& fileassoc, const Int1D& colidx, double clump_p2);

std::vector<UMapIntString> map_assoc_file(const std::string& fileassoc, const Int1D& colidx);

void ld_r2_big(const Mat2D& G, const SNPld& snp, const std::string& filebim, const std::string& fileout);

void ld_r2_small(Data* data, const SNPld& snp, const std::string& filebim, const std::string& fileout);

void run_ld_stuff(Data* data, const Param& params);

#endif  // PCAONE_LD_
