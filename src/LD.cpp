/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/LD.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "LD.hpp"

#include <zlib.h>

#include <cstddef>

#include "Cmd.hpp"
#include "Data.hpp"
#include "Utils.hpp"

using namespace std;

namespace {
bool is_pvar_header(const std::string& line) { return !line.empty() && line[0] == '#'; }

String1D variant_tokens_for_bimish_line(const std::string& line, const std::string& path) {
  if (is_pvar_header(line)) return {};
  auto tokens = split_string(line, " \t");
  if ((int)tokens.size() >= 5 && path.size() >= 5 && path.substr(path.size() - 5) == ".pvar") {
    return split_string(pvar_line_to_bim_line(line, path), " \t");
  }
  return tokens;
}

// CHR, BP and SNP of each variant, as the .ld.gz columns want them
String1D variant_labels(const String1D& variants) {
  String1D labels;
  labels.reserve(variants.size());
  for (const auto& line : variants) {
    auto tokens = split_string(line, " \t");
    labels.push_back(tokens[0] + "\t" + tokens[3] + "\t" + tokens[1]);
  }
  return labels;
}

// append one variant; chromosomes are assumed to be contiguous
void push_snp_pos(SNPld& snp, const std::string& chr, int pos, std::string& chr_prev, int& i) {
  if (chr_prev.empty()) snp.chr.push_back(chr);
  // when starting a new chromosome
  if (!chr_prev.empty() && chr_prev != chr) {
    snp.end_pos.push_back(i - 1);
    snp.chr.push_back(chr);
  }
  chr_prev = chr;
  snp.pos.push_back(pos);
  i++;
}

// every out-of-core read in the LD code goes through here, so a block is never
// used before the PCs are removed from it
void read_ld_block(Data* data, uint b, const Mat2D& Q) {
  data->read_block_initial(data->start[b], data->stop[b], false);
  adjust_for_pcs(data->G, Q);
}
}  // namespace

// compute sample standard deviation
Arr1D calc_sds(const Mat2D& X) {
  const double df = 1.0 / (X.rows() - 1);  // N-1
  return (X.array().square().colwise().sum() * df).sqrt();
}

// Inverse sample standard deviation per column, with zero-variance columns set
// to 0 instead of infinity.
//
// A variant that is monomorphic, or whose genotypes lie entirely in the span of
// the PCs being adjusted for, has no residual variance. Taking 1/sd gives inf,
// and the correlation inf * inf * 0 is NaN, which then propagates into every
// pair involving that variant and into the .ld.gz output. Reporting 0 instead
// says what is actually true -- a variant with no variance carries no LD
// information -- and keeps the output numeric.
Arr1D calc_inv_sds(const Mat2D& X, Eigen::Index& nzero) {
  const Arr1D sds = calc_sds(X);
  Arr1D inv(sds.size());
  nzero = 0;
  for (Eigen::Index i = 0; i < sds.size(); ++i) {
    if (sds(i) > VAR_TOL) {
      inv(i) = 1.0 / sds(i);
    } else {
      inv(i) = 0.0;
      ++nzero;
    }
  }
  return inv;
}

// shared warning so every LD entry point reports this the same way
static void warn_zero_variance(Eigen::Index nzero, Eigen::Index ntotal) {
  if (nzero > 0)
    cao.warn(nzero, " of ", ntotal,
             " variants have no variance after ancestry adjustment (monomorphic, or fully explained by the "
             "PCs). their R2 is reported as 0 rather than NaN. consider --maf to remove them.");
}

// compute the peason correlation coefficient
// should check x.size()==y.size()
double calc_cor(const Mat1D& x, const Mat1D& y, const double df) {
  int N = x.size();
  double numerator = x.dot(y);
  // Calculate variances with N-1 degrees of freedom
  double var_x = x.dot(x) / (N - 1);
  double var_y = y.dot(y) / (N - 1);
  // const double epsilon = std::numeric_limits<double>::epsilon();
  // const double min_variance = epsilon * epsilon;
  // if (var_x < min_variance || var_y < min_variance) {
  //   if (var_x < min_variance && var_y < min_variance) {
  //     // Both variances are essentially zero
  //     // Perfect correlation if both vectors are constant
  //     return 1.0;
  //   } else {
  //     // No correlation if one variance is essentially zero, the other is not
  //     return 0.0;
  //   }
  // }

  // no variance, no LD information: 0 rather than NaN, as calc_inv_sds() does
  if (!(std::sqrt(var_x) > VAR_TOL) || !(std::sqrt(var_y) > VAR_TOL)) return 0.0;
  double sd_x = 1.0 / std::sqrt(var_x);
  double sd_y = 1.0 / std::sqrt(var_y);
  return numerator * sd_x * sd_y * df;
}

std::string get_snp_pos_bim(SNPld& snp, const std::string& filebim, bool header, Int1D idx) {
  std::ifstream fin(filebim);
  if (!fin.is_open()) cao.error("can not open " + filebim);
  std::string ret, line{""}, chr_prev;
  int i = 0;
  if (header) getline(fin, line);
  ret = line;
  while (getline(fin, line)) {
    if (line.empty() || is_pvar_header(line)) continue;
    auto tokens = variant_tokens_for_bimish_line(line, filebim);
    push_snp_pos(snp, tokens[idx[0]], std::stoi(tokens[idx[1]]), chr_prev, i);
  }
  snp.end_pos.push_back(i - 1);  // add the last SNP
  return ret;
}

void get_snp_pos_variants(SNPld& snp, const String1D& variants) {
  std::string chr_prev;
  int i = 0;
  for (const auto& line : variants) {
    auto tokens = split_string(line, " \t");
    push_snp_pos(snp, tokens[0], std::stoi(tokens[3]), chr_prev, i);
  }
  snp.end_pos.push_back(i - 1);  // add the last SNP
}

String1D read_ld_variants(const Data* data, const Param& params) {
  const bool is_pvar = params.file_t == FileType::PGEN;
  const std::string path = params.filein + (is_pvar ? ".pvar" : ".bim");
  std::ifstream fin(path);
  if (!fin.is_open()) cao.error("can not open " + path);
  const Int1D& keep = data->keepSNPs;  // empty unless --maf removed sites
  String1D variants;
  variants.reserve(data->nsnps);
  std::string line;
  for (int idx = 0, kept = 0; getline(fin, line);) {
    if (line.empty() || is_pvar_header(line)) continue;
    if (keep.empty() || (kept < (int)keep.size() && keep[kept] == idx)) {
      variants.push_back(is_pvar ? pvar_line_to_bim_line(line, path) : line);
      if ((int)split_string(variants.back(), " \t").size() < 6)
        cao.error("the input variant file is not valid!\n => " + path);
      kept++;
    }
    idx++;
  }
  if (variants.size() != data->nsnps)
    cao.error(path, "lists", variants.size(), "variants but the genotype matrix has", data->nsnps);
  return variants;
}

Mat2D read_ld_pcs(const Param& params, uint nsamples) {
  if (params.ld_stats == 1) {
    if (!params.fileU.empty()) cao.warn("--ld-stats 1 is the standard LD. the PCs in", params.fileU, "are not used");
    cao.print(tick.date(), "compute the standard LD, without ancestry adjustment");
    return Mat2D();
  }
  const Mat2D U = read_usv(params.fileU);
  if (U.rows() != nsamples)
    cao.error(params.fileU, "has", U.rows(), "rows but the genotypes have", nsamples,
              "samples. the PCs must come from the same samples, in the same order");
  // .eigvecs is text with 6 significant digits, so its columns are orthonormal
  // only to ~1e-6. Q spans the same PCs and is orthonormal to machine
  // precision, which makes I - QQ' an exact projector.
  Eigen::HouseholderQR<Mat2D> qr(U);
  Mat2D Q = qr.householderQ() * Mat2D::Identity(U.rows(), U.cols());
  cao.print(tick.date(), "compute the ancestry adjusted LD, removing", U.cols(), "PCs in", params.fileU);
  return Q;
}

// G holds centred genotypes (missing calls imputed to the site mean), so the
// residuals of regressing each site on the PCs are (I - QQ')G. The PCs come from
// centred data and are orthogonal to the intercept, so the residuals stay
// centred; re-centring only removes rounding. Per-site scaling commutes with
// the projection, so whether the PCA standardised the sites does not matter.
void adjust_for_pcs(Mat2D& G, const Mat2D& Q) {
  if (Q.size() == 0) return;
  G.noalias() -= Q * (Q.transpose() * G);
  G.rowwise() -= G.colwise().mean();
}

// given a list of snps, find its index per chr in the original pos
// assume chromosomes are continuous
// TODO: check duplicated POS
std::tuple<Int2D, Int2D> get_target_snp_idx(const SNPld& snp_t, const SNPld& snp) {
  cao.print(tick.date(), "try to match target SNPs to the SNPs in LD matrix");
  Int1D idx, bp;
  UMapIntInt mpos;
  int c, s, e, p, i;
  Int1D ord;
  for (c = 0; c < (int)snp.chr.size(); c++) {
    i = 0;
    while (snp_t.chr[i] != snp.chr[c]) {
      i++;
      if (i == (int)snp_t.chr.size()) break;
    }
    ord.push_back(i);
  }
  if (!std::is_sorted(ord.begin(), ord.end()))
    cao.error(
        "the association file may be not sorted, hence not matching the bim"
        "file ");

  Int2D idx_per_chr(snp_t.chr.size());
  Int2D bp_per_chr(snp_t.chr.size());
  for (int tc = 0; tc < (int)snp_t.chr.size(); tc++) {
    for (c = 0; c < (int)snp.chr.size(); c++)
      if (snp.chr[c] == snp_t.chr[tc]) break;
    if (c == (int)snp.chr.size()) continue;  // a chromosome the genotypes do not have
    e = snp.end_pos[c];
    s = c > 0 ? snp.end_pos[c - 1] : 0;
    for (i = s; i <= e; i++) mpos[snp.pos[i]] = i;
    e = snp_t.end_pos[tc];
    s = tc > 0 ? snp_t.end_pos[tc - 1] : 0;
    for (i = s; i <= e; i++) {
      p = snp_t.pos[i];
      if (mpos.count(p)) {
        bp.push_back(p);
        idx.push_back(mpos[p]);
      }
    }
    idx_per_chr[tc] = idx;
    bp_per_chr[tc] = bp;
    bp.clear();
    idx.clear();
    mpos.clear();
  }
  return std::make_tuple(idx_per_chr, bp_per_chr);
}

/// return ws, we
void divide_pos_by_window(SNPld& snp, const int ld_window_bp) {
  int nsnp = snp.pos.size();
  int j{0}, c{0}, nsites;
  for (int i = 0; i < nsnp; i++) {
    if (snp.pos[i] == snp.pos[snp.end_pos[c]]) {
      c++;
      continue;
    }
    for (j = i; j <= snp.end_pos[c]; j++)
      if (snp.pos[j] - snp.pos[i] > ld_window_bp) break;
    nsites = j - i;
    snp.ws.push_back(i);       // start pos in the window
    snp.we.push_back(nsites);  // the number of sites
  }
}

void write_pruned_snp_ids(const String1D& variants, const std::string& fileout, const ArrBool& keep) {
  cao.print(tick.date(), keep.count(), " sites will be kept");
  std::ofstream ofs_out(fileout + ".ld.prune.out");
  std::ofstream ofs_in(fileout + ".ld.prune.in");
  int i = 0;
  for (const auto& line : variants) {
    const auto fields = split_string(line, " \t");
    if (keep(i))
      ofs_in << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << fields[4] << "\t"
             << fields[5] << std::endl;
    else
      ofs_out << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << fields[4] << "\t"
              << fields[5] << std::endl;
    i++;
  }
}

// of each pair above the cutoff, the site with the lower MAF is removed. F is
// read as the blocks are: every site compared has been read by then.
void ld_prune_small(
    Data* data, const Mat2D& Q, const String1D& variants, const SNPld& snp, double r2_tol, const std::string& fileout) {
  cao.print(tick.date(), "LD pruning, the site with the lower MAF of each pair is removed");
  const Mat1D& F = data->F;
  ArrBool keep = ArrBool::Constant(data->nsnps, true);
  const double df = 1.0 / (data->nsamples - 1);  // N-1
  uint b = 0, start = 0, end = 0;
  data->check_file_offset_first_var();
  read_ld_block(data, b, Q);
  Mat2D G = data->G;  // make a copy. we need two G anyway
  b++;
  read_ld_block(data, b, Q);
  for (size_t w = 0; w < snp.ws.size(); w++) {
    uint i = snp.ws[w];
    if (!keep(i)) continue;
    start = i;  // ensure start >= data->start[b-1]
    if (start < data->start[b - 1]) cao.error("BUG: ld_prune_small ");
    end = snp.we[w] + i - 1;
    if (end > data->stop[b]) {  // renew G and data->G
      G = data->G;              // copy
      b++;
      read_ld_block(data, b, Q);
    }
#pragma omp parallel for
    for (int j = 1; j < snp.we[w]; j++) {
      uint k = i + j;
      if (!keep(k)) continue;
      double r = 0;
      if (i >= data->start[b - 1] && k <= data->stop[b - 1]) {
        r = calc_cor(G.col(i - data->start[b - 1]), G.col(k - data->start[b - 1]), df);
      } else if (i >= data->start[b] && k <= data->stop[b]) {
        r = calc_cor(data->G.col(i - data->start[b]), data->G.col(k - data->start[b]), df);
      } else {
        r = calc_cor(G.col(i - data->start[b - 1]), data->G.col(k - data->start[b]), df);
      }
      if (r * r > r2_tol) {
        const uint o = MAF(F(k)) > MAF(F(i)) ? i : k;
        keep(o) = false;
      }
    }
  }
  write_pruned_snp_ids(variants, fileout, keep);
}

void ld_prune_big(const Mat2D& G,
                  const Mat1D& F,
                  const String1D& variants,
                  const SNPld& snp,
                  double r2_tol,
                  const std::string& fileout) {
  if ((long int)snp.pos.size() != G.cols()) cao.error("The number of variants is not matching the LD matrix");
  cao.print(tick.date(), "LD pruning, the site with the lower MAF of each pair is removed");
  Eigen::Index nzero = 0;
  Arr1D sds = calc_inv_sds(G, nzero);
  warn_zero_variance(nzero, G.cols());
  ArrBool keep = ArrBool::Constant(G.cols(), true);
  const double df = 1.0 / (G.rows() - 1);  // N-1
  for (int w = 0; w < (int)snp.ws.size(); w++) {
    int i = snp.ws[w];
    if (!keep(i)) continue;
#pragma omp parallel for
    for (int j = 1; j < snp.we[w]; j++) {
      int k = i + j;
      if (!keep(k)) continue;
      double r = G.col(i).dot(G.col(k)) * (sds(i) * sds(k) * df);
      if (r * r > r2_tol) {
        const int o = MAF(F(k)) > MAF(F(i)) ? i : k;
        keep(o) = false;
      }
    }
  }
  write_pruned_snp_ids(variants, fileout, keep);
}

Int1D valid_assoc_file(const std::string& fileassoc, const std::string& colnames) {
  std::ifstream fin(fileassoc);
  if (!fin.is_open()) cao.error("can not open " + fileassoc);
  std::string line, sep{"\t"}, sep2{","};
  getline(fin, line);
  const auto fields = split_string(line, sep);
  std::vector<std::string> fields_users{"CHR", "BP", "P"};
  if (!colnames.empty()) fields_users = split_string(colnames, sep2);
  Int1D idx(3, -1);
  int j = 0;
  for (auto col : fields) {
    if (col == fields_users[0]) {
      idx[0] = j;
    } else if (col == fields_users[1]) {
      idx[1] = j;
    } else if (col == fields_users[2]) {
      idx[2] = j;
    }
    j++;
  }
  j = 0;
  for (auto i : idx) {
    if (i < 0) cao.error("the assoc-like file has no " + fields_users[j] + " column");
    j++;
  }
  return idx;
}

std::vector<UMapIntPds> map_index_snps(const std::string& fileassoc, const Int1D& colidx, double clump_p2) {
  std::ifstream fin(fileassoc);
  if (!fin.is_open()) cao.error("can not open " + fileassoc);
  std::string line, chr_cur, chr_prev, sep{"\t"};
  getline(fin, line);
  vector<UMapIntPds> vm;
  UMapIntPds m;
  int bp;
  double pval;
  while (getline(fin, line)) {
    auto tokens = split_string(line, sep);
    chr_cur = tokens[colidx[0]];
    bp = std::stoi(tokens[colidx[1]]);
    pval = std::stod(tokens[colidx[2]]);
    if (pval <= clump_p2) m.insert({bp, {pval, line}});
    if (!chr_prev.empty() && chr_prev != chr_cur) {
      vm.push_back(m);
      m.clear();
    }
    chr_prev = chr_cur;
  }
  vm.push_back(m);  // add the last chr
  return vm;
}

void ld_clump_single_pheno(const std::string& fileout,
                           const std::string& head,
                           const int clump_bp,
                           const double clump_r2,
                           const double clump_p1,
                           const double clump_p2,
                           const Mat2D& G,
                           const Int2D& idx_per_chr,
                           const Int2D& bp_per_chr,
                           const std::vector<UMapIntPds>& pvals_per_chr) {
  // sort by pvalues and get new idx
  Eigen::Index nzero = 0;
  const Arr1D sds = calc_inv_sds(G, nzero);
  warn_zero_variance(nzero, G.cols());
  const double df = 1.0 / (G.rows() - 1);  // N-1
  std::ofstream ofs(fileout);
  ofs << head + "\tSP2" << std::endl;
  for (int c = 0; c < (int)bp_per_chr.size(); c++) {
    const auto idx = idx_per_chr[c];
    const auto bp = bp_per_chr[c];
    const auto mbp = vector2map(bp);
    // greedy clumping algorithm
    auto mpp = pvals_per_chr[c];  // key: pos, val: pval
    Double1D pp;
    Int1D ps;
    for (auto it = mpp.begin(); it != mpp.end(); it++) {
      if (it->second.first <= clump_p1) {
        ps.push_back(it->first);
        pp.push_back(it->second.first);
      }
    }
    int p, p2, j, k;
    for (auto i : sortidx(pp)) {  // snps sorted by p value
      p = ps[i];
      if (mpp.count(p) == 0) continue;  // if snps with pval < clump_p1 are already clumped
      Int1D clumped;
      bool backward = true;
      if (mbp.count(p) == 0) continue;
      j = mbp.at(p);  // j:current
      k = j;          // k:forward or backward
      while (true) {
        if (backward) {
          --k;
          if (k < 0 || (bp[k] < p - clump_bp)) {
            backward = false;
            k = j;
            continue;
          }
        } else {
          ++k;
          if (k >= (int)bp.size() || (bp[k] > p + clump_bp)) break;
        }
        p2 = bp[k];
        if (mpp.count(p2) == 0) continue;
        double r = G.col(idx[j]).dot(G.col(idx[k])) * (sds(idx[j]) * sds(idx[k]) * df);
        if (r * r >= clump_r2) {
          clumped.push_back(p2);
          mpp.erase(p2);
        }
      }
      // what we do with clumped SNPs. sort them by pval?
      ofs << pvals_per_chr[c].at(p).second << "\t";
      if (clumped.empty()) {
        ofs << "NONE";
      } else {
        Double1D opp;
        for (auto op : clumped) opp.push_back(pvals_per_chr[c].at(op).first);
        k = 0;
        for (auto oi : sortidx(opp)) {
          if (k == (int)opp.size() - 1)
            ofs << clumped[oi];
          else
            ofs << clumped[oi] << ",";
          k++;
        }
      }
      ofs << std::endl;
    }
    // end current chr
  }
}

void ld_r2_small(Data* data, const Mat2D& Q, const String1D& variants, const SNPld& snp, const std::string& fileout) {
  const String1D bims = variant_labels(variants);
  const double df = 1.0 / (data->nsamples - 1);  // N-1
  int b = 0;
  data->check_file_offset_first_var();
  read_ld_block(data, b, Q);
  Mat2D G = data->G;  // make a copy. we need two G anyway
  b++;
  read_ld_block(data, b, Q);

  gzFile gzfp = gzopen(fileout.c_str(), "wb");
  std::string line{"CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2\n"};
  gzwrite(gzfp, line.c_str(), line.size());

  for (size_t w = 0; w < snp.ws.size(); w++) {
    uint i = snp.ws[w];
    if (snp.we[w] + i - 1 > data->stop[b]) {  // renew G and data->G
      G = data->G;                            // copy
      b++;
      read_ld_block(data, b, Q);
    }

    if (data->params.verbose) cao.print(tick.date(), "process window", w, ", the index SNP is", i);

    for (int j = 1; j < snp.we[w]; j++) {
      uint k = i + j;
      double r = 0;
      if (i >= data->start[b - 1] && k <= data->stop[b - 1]) {
        r = calc_cor(G.col(i - data->start[b - 1]), G.col(k - data->start[b - 1]), df);
      } else if (i >= data->start[b] && k <= data->stop[b]) {
        r = calc_cor(data->G.col(i - data->start[b]), data->G.col(k - data->start[b]), df);

      } else {
        r = calc_cor(G.col(i - data->start[b - 1]), data->G.col(k - data->start[b]), df);
      }
      //// FIXME: this would be slow
      line = bims[i] + "\t" + bims[k] + "\t" + std::to_string(r * r) + "\n";

      if (gzwrite(gzfp, line.c_str(), line.size()) != static_cast<int>(line.size()))
        cao.error("failed to write data to ld.gz file");
    }
  }
  gzclose(gzfp);
}

void ld_r2_big(const Mat2D& G, const String1D& variants, const SNPld& snp, const std::string& fileout) {
  const String1D bims = variant_labels(variants);
  Eigen::Index nzero = 0;
  Arr1D sds = calc_inv_sds(G, nzero);
  warn_zero_variance(nzero, G.cols());
  const double df = 1.0 / (G.rows() - 1);  // N-1
  gzFile gzfp = gzopen(fileout.c_str(), "wb");
  std::string line{"CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2\n"};
  gzwrite(gzfp, line.c_str(), line.size());

  for (int w = 0; w < (int)snp.ws.size(); w++) {
    int i = snp.ws[w];

    for (int j = 1; j < snp.we[w]; j++) {
      int k = i + j;
      double r = G.col(i).dot(G.col(k)) * (sds(i) * sds(k) * df);
      //// FIXME: this would be slow
      line = bims[i] + "\t" + bims[k] + "\t" + std::to_string(r * r) + "\n";
      if (gzwrite(gzfp, line.c_str(), line.size()) != static_cast<int>(line.size()))
        cao.error("failed to write data to ld.gz file");
    }
  }
  gzclose(gzfp);
}

// The LD statistics are correlations between the columns of (I - QQ')G, where G
// is the centred genotypes and Q the PCs of --USV (none for --ld-stats 1). The
// residuals are formed as the genotypes are read, the whole matrix in-core or
// block by block with -m, so no residual matrix is ever written to disk.
void run_ld_stuff(Data* data, const Param& params) {
  cao.print(tick.date(), "run LD stuff");
  // read the PCs before the genotypes, so a mismatched .eigvecs fails fast
  const Mat2D Q = read_ld_pcs(params, data->nsamples);
  data->prepare();
  const String1D variants = read_ld_variants(data, params);
  SNPld snp;  // SNPs information for LD prunning
  get_snp_pos_variants(snp, variants);
  if (!params.out_of_core) adjust_for_pcs(data->G, Q);

  if (params.clump.empty()) {
    divide_pos_by_window(snp, params.ld_bp);

    if (params.out_of_core) {
      if (params.print_r2) {
        ld_r2_small(data, Q, variants, snp, params.fileout + ".ld.gz");
      } else {
        ld_prune_small(data, Q, variants, snp, params.ld_r2, params.fileout);
      }
    } else {
      if (params.print_r2) {
        ld_r2_big(data->G, variants, snp, params.fileout + ".ld.gz");
      } else {
        ld_prune_big(data->G, data->F, variants, snp, params.ld_r2, params.fileout);
      }
    }

  } else {
    const auto assocfiles = split_string(params.clump, ",");
    for (size_t i = 0; i < assocfiles.size(); i++) {
      cao.print(tick.date(), "LD-based clumping for associated file:", assocfiles[i]);
      auto colidx = valid_assoc_file(assocfiles[i], params.assoc_colnames);
      SNPld snp_t;
      std::string head = get_snp_pos_bim(snp_t, assocfiles[i], true, colidx);
      const auto pvals_per_chr = map_index_snps(assocfiles[i], colidx, params.clump_p2);
      Int2D idx_per_chr, bp_per_chr;
      std::tie(idx_per_chr, bp_per_chr) = get_target_snp_idx(snp_t, snp);
      if (params.out_of_core) {
        Mat2D G(data->nsamples, snp_t.pos.size());
        int b = 0;
        data->check_file_offset_first_var();
        read_ld_block(data, b, Q);
        for (int sidx = 0, c = 0; c < (int)idx_per_chr.size(); c++) {
          int i = 0;
          for (auto icol : idx_per_chr[c]) {
            // the target sites are sorted, but may skip a whole block
            while (icol > (int)data->stop[b]) {
              b++;
              read_ld_block(data, b, Q);
            }
            G.col(sidx) = data->G.col(icol - data->start[b]);
            idx_per_chr[c][i] = sidx;
            sidx++;
            i++;
          }
        }

        ld_clump_single_pheno(params.fileout + ".p" + std::to_string(i) + ".clump", head, params.clump_bp,
                              params.clump_r2, params.clump_p1, params.clump_p2, G, idx_per_chr, bp_per_chr,
                              pvals_per_chr);

      } else {
        ld_clump_single_pheno(params.fileout + ".p" + std::to_string(i) + ".clump", head, params.clump_bp,
                              params.clump_r2, params.clump_p1, params.clump_p2, data->G, idx_per_chr, bp_per_chr,
                              pvals_per_chr);
      }
    }
  }
}
