/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/LD.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "LD.hpp"

#include "Cmd.hpp"
#include "Data.hpp"
#include "Utils.hpp"

using namespace std;

MyArray calc_sds(const MyMatrix& X) {
  // compute degree of freedom
  const double df = 1.0 / (X.rows() - 1);  // N-1
  return (X.array().square().colwise().sum() * df).sqrt();
}

// df = 1.0 / (N - 1)
// should check x.size()==y.size()
double calc_cor(const MyVector& x, const MyVector& y, const double df) {
  int N = x.size();
  double numerator = x.dot(y);
  double sd_x = 1.0 / std::sqrt(x.dot(x) / (N - 1));
  double sd_y = 1.0 / std::sqrt(y.dot(y) / (N - 1));
  return numerator * sd_x * sd_y * df;
}

std::string get_snp_pos_bim(SNPld& snp, const std::string& filebim, bool header,
                            Int1D idx) {
  std::ifstream fin(filebim);
  if (!fin.is_open()) cao.error("can not open " + filebim);
  std::string ret, line{""}, chr_cur, chr_prev, sep{" \t"};
  int i = 0;
  if (header) getline(fin, line);
  ret = line;
  while (getline(fin, line)) {
    auto tokens = split_string(line, sep);
    if ((int)tokens.size() == 7 && !header)
      snp.af.push_back(std::stod(tokens[6]));
    chr_cur = tokens[idx[0]];
    if (chr_prev.empty()) snp.chr.push_back(chr_cur);
    if (!chr_prev.empty() && chr_prev != chr_cur) {
      snp.end_pos.push_back(i - 1);
      snp.chr.push_back(chr_cur);
    }
    chr_prev = chr_cur;
    snp.pos.push_back(std::stoi(tokens[idx[1]]));
    i++;
  }
  snp.end_pos.push_back(i - 1);  // add the last SNP
  return ret;
}

// given a list of snps, find its index per chr in the original pos
// assume chromosomes are continuous
// TODO: check duplicated POS
std::tuple<Int2D, Int2D> get_target_snp_idx(const SNPld& snp_t,
                                            const SNPld& snp) {
  Int1D idx, bp;
  UMapIntInt mpos;
  int c, s, e, p, i;
  Int2D ret(snp_t.chr.size());
  Int2D ret2(snp_t.chr.size());
  for (int tc = 0; tc < (int)snp_t.chr.size(); tc++) {
    for (c = 0; c < (int)snp.chr.size(); c++)
      if (snp.chr[c] == snp_t.chr[tc]) break;
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
    ret[tc] = idx;
    ret2[tc] = bp;
    bp.clear();
    idx.clear();
    mpos.clear();
  }
  return std::make_tuple(ret, ret2);
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

void ld_prune_small(Data* data, const std::string& fileout,
                    const std::string& filebim, const SNPld& snp,
                    const double r2_tol) {
  const bool pick_random_one = snp.af.size() > 0 ? false : true;
  cao.print(tick.date(), "LD pruning, pick_random_one =", pick_random_one);
  ArrayXb keep = ArrayXb::Constant(data->nsnps, true);
  const double df = 1.0 / (data->nsamples - 1);  // N-1
  int b = 0, start = 0, end = 0;
  data->check_file_offset_first_var();
  data->read_block_initial(data->start[b], data->stop[b], false);
  MyMatrix G = data->G;  // make a copy. we need two G anyway
  b++;
  data->read_block_initial(data->start[b], data->stop[b], false);
  for (int w = 0; w < (int)snp.ws.size(); w++) {
    int i = snp.ws[w];
    if (!keep(i)) continue;
    start = i;  // ensure start >= data->start[b-1]
    if (start < data->start[b - 1]) cao.error("BUG: ld_prune_small ");
    end = snp.we[w] + i - 1;
    if (end > data->stop[b]) {  // renew G and data->G
      G = data->G;              // copy
      b++;
      data->read_block_initial(data->start[b], data->stop[b], false);
    }
#pragma omp parallel for
    for (int j = 1; j < snp.we[w]; j++) {
      int k = i + j;
      if (!keep(k)) continue;
      double r = 0;
      if (i >= data->start[b - 1] && k <= data->stop[b - 1]) {
        r = calc_cor(G.col(i - data->start[b - 1]),
                     G.col(k - data->start[b - 1]), df);
      } else if (i >= data->start[b] && k <= data->stop[b]) {
        r = calc_cor(data->G.col(i - data->start[b]),
                     data->G.col(k - data->start[b]), df);

      } else {
        r = calc_cor(G.col(i - data->start[b - 1]),
                     data->G.col(k - data->start[b]), df);
      }
      if (r * r > r2_tol) {
        int o = k;  // or i
        if (!pick_random_one) o = MAF(snp.af[k]) > MAF(snp.af[i]) ? i : k;
        keep(o) = false;
      }
    }
  }
  std::ifstream fin(filebim);
  if (!fin.is_open()) cao.error("can not open " + filebim);
  std::ofstream ofs_out(fileout + ".ld.prune.out");
  std::ofstream ofs_in(fileout + ".ld.prune.in");
  std::string line;
  int i = 0;
  while (getline(fin, line)) {
    if (keep(i))
      ofs_in << line << std::endl;
    else
      ofs_out << line << std::endl;
    i++;
  }
}

void ld_prune_big(const std::string& fileout, const std::string& filebim,
                  const MyMatrix& G, const SNPld& snp, double r2_tol) {
  if ((long int)snp.pos.size() != G.cols())
    cao.error("The number of variants is not matching the LD matrix");
  // TODO: maybe add an option in CLI
  const bool pick_random_one = snp.af.size() > 0 ? false : true;
  cao.print(tick.date(), "LD pruning, pick_random_one =", pick_random_one);
  MyArray sds = 1.0 / calc_sds(G);
  ArrayXb keep = ArrayXb::Constant(G.cols(), true);
  const double df = 1.0 / (G.rows() - 1);  // N-1
  for (int w = 0; w < (int)snp.ws.size(); w++) {
    int i = snp.ws[w];
    if (!keep(i)) continue;
#pragma omp parallel for
    for (int j = 1; j < snp.we[w]; j++) {
      int k = i + j;
      if (!keep(k)) continue;
      // double r = calc_cor2(G.col(i), G.col(k), df);
      double r = G.col(i).dot(G.col(k)) * (sds(i) * sds(k) * df);
      // (G.col(i).array() * G.col(k).array() * (sds(i) * sds(k))).sum() * df;
      if (r * r > r2_tol) {
        int o = k;  // or i
        if (!pick_random_one) o = MAF(snp.af[k]) > MAF(snp.af[i]) ? i : k;
        keep(o) = false;
      }
    }
  }
  std::ifstream fin(filebim);
  if (!fin.is_open()) cao.error("can not open " + filebim);
  std::ofstream ofs_out(fileout + ".ld.prune.out");
  std::ofstream ofs_in(fileout + ".ld.prune.in");
  std::string line;
  int i = 0;
  while (getline(fin, line)) {
    if (keep(i))
      ofs_in << line << std::endl;
    else
      ofs_out << line << std::endl;
    i++;
  }
}

Int1D valid_assoc_file(const std::string& fileassoc,
                       const std::string& colnames) {
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
    if (i < 0)
      cao.error("the assoc-like file has no " + fields_users[j] + " column");
    j++;
  }
  return idx;
}

std::vector<UMapIntPds> map_index_snps(const std::string& fileassoc,
                                       const Int1D& colidx, double clump_p2) {
  std::ifstream fin(fileassoc);
  if (!fin.is_open()) cao.error("can not open " + fileassoc);
  std::string line, chr_cur, chr_prev, sep{"\t"};
  getline(fin, line);
  vector<UMapIntPds> vm;
  UMapIntPds m;
  int c = 0, bp;
  double pval;
  while (getline(fin, line)) {
    auto tokens = split_string(line, sep);
    chr_cur = tokens[colidx[0]];
    bp = std::stoi(tokens[colidx[1]]);
    pval = std::stod(tokens[colidx[2]]);
    if (pval <= clump_p2) m.insert({bp, {pval, line}});
    if (!chr_prev.empty() && chr_prev != chr_cur) {
      c++;
      vm.push_back(m);
      m.clear();
    }
    chr_prev = chr_cur;
  }
  vm.push_back(m);  // add the last chr
  return vm;
}

void ld_clump_big(std::string fileout, std::string fileassoc,
                  std::string colnames, int clump_bp, double clump_r2,
                  double clump_p1, double clump_p2, const SNPld& snp,
                  const MyMatrix& G) {
  if ((long int)snp.pos.size() != G.cols())
    cao.error("The number of variants is not matching the LD matrix");
  cao.print(tick.date(), "LD-based clumping: p1 =", clump_p1,
            ", p2 =", clump_p2, ", r2 =", clump_r2, ", bp =", clump_bp,
            ", associated file:", fileassoc);
  // 0: chr, 1: pos, 2: pvalue
  auto colidx = valid_assoc_file(fileassoc, colnames);
  SNPld snp_t;
  std::string head = get_snp_pos_bim(snp_t, fileassoc, true, colidx);

  Int2D idx_per_chr, bp_per_chr;
  std::tie(idx_per_chr, bp_per_chr) = get_target_snp_idx(snp_t, snp);
  const auto pvals_per_chr = map_index_snps(fileassoc, colidx, clump_p2);
  // sort by pvalues and get new idx
  const MyArray sds = 1.0 / calc_sds(G);
  const double df = 1.0 / (G.rows() - 1);  // N-1
  std::ofstream ofs(fileout);
  ofs << head + "\tSP2" << std::endl;
  for (int c = 0; c < (int)idx_per_chr.size(); c++) {
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
      if (mpp.count(p) == 0)
        continue;  // if snps with pval < clump_p1 are already clumped with
      // previous snps
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
        double r =
            G.col(idx[j]).dot(G.col(idx[k])) * (sds(idx[j]) * sds(idx[k]) * df);
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
          if (k == opp.size() - 1)
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

void run_ld_stuff(const Param& params, Data* data) {
  cao.print(tick.date(), "run ld stuff");
  SNPld snp;  // SNPs information for LD prunning
  get_snp_pos_bim(snp, params.filebim);
  if (params.clump.empty()) {
    divide_pos_by_window(snp, params.ld_bp);
    if (params.out_of_core) {
      ld_prune_small(data, params.fileout, params.filebim, snp, params.ld_r2);
    } else {
      ld_prune_big(params.fileout, params.filebim, data->G, snp, params.ld_r2);
    }
  } else {
    const auto assocfiles = split_string(params.clump, ",");
    for (size_t i = 0; i < assocfiles.size(); i++) {
      // TODO: subset G for each pheno and compared with it
      if (params.out_of_core) {
        auto colidx = valid_assoc_file(assocfiles[i], params.assoc_colnames);
        Int2D idx_per_chr, bp_per_chr;
        SNPld snp_t;
        std::string head = get_snp_pos_bim(snp_t, assocfiles[i], true, colidx);
        std::tie(idx_per_chr, bp_per_chr) = get_target_snp_idx(snp_t, snp);
        const auto pvals_per_chr =
            map_index_snps(assocfiles[i], colidx, params.clump_p2);
        MyMatrix G(data->nsamples, snp_t.pos.size());
        int b = 0, sidx = 0;
        data->check_file_offset_first_var();
        data->read_block_initial(data->start[b], data->stop[b], false);
        Int1D cumlen(idx_per_chr.size());
        int cl = 0;
        for (int c = 0; c < (int)idx_per_chr.size(); c++) {
          for (auto icol : idx_per_chr[c]) {
            if (icol >= data->start[b] && icol <= data->stop[b]) {
              G.col(sidx) = data->G.col(icol - data->start[b]);
            } else {
              b++;
              data->read_block_initial(data->start[b], data->stop[b], false);
              G.col(sidx) = data->G.col(icol - data->start[b]);
            }
            sidx++;
          }
          cl += idx_per_chr[c].size();
          cumlen[c] = cl;
        }
        const MyArray sds = 1.0 / calc_sds(G);
        const double df = 1.0 / (G.rows() - 1);  // N-1
        std::string fileout =
            params.fileout + ".p" + std::to_string(i) + ".clump";
        std::ofstream ofs(fileout);
        ofs << head + "\tSP2" << std::endl;
        for (int c = 0; c < (int)idx_per_chr.size(); c++) {
          const auto bp = bp_per_chr[c];
          const auto mbp = vector2map(bp);
          // greedy clumping algorithm
          auto mpp = pvals_per_chr[c];  // key: pos, val: pval
          Double1D pp;
          Int1D ps;
          for (auto it = mpp.begin(); it != mpp.end(); it++) {
            if (it->second.first <= params.clump_p1) {
              ps.push_back(it->first);
              pp.push_back(it->second.first);
            }
          }
          int p, p2, j, k;
          for (auto i : sortidx(pp)) {  // snps sorted by p value
            p = ps[i];
            if (mpp.count(p) == 0)
              continue;  // if snps with pval < clump_p1 are already clumped
                         // with
            // previous snps
            Int1D clumped;
            bool backward = true;
            if (mbp.count(p) == 0) continue;
            j = mbp.at(p);  // j:current
            k = j;          // k:forward or backward
            while (true) {
              if (backward) {
                --k;
                if (k < 0 || (bp[k] < p - params.clump_bp)) {
                  backward = false;
                  k = j;
                  continue;
                }
              } else {
                ++k;
                if (k >= (int)bp.size() || (bp[k] > p + params.clump_bp)) break;
              }
              p2 = bp[k];
              if (mpp.count(p2) == 0) continue;
              int cj = c > 0 ? cumlen[c - 1] + j - 1 : j;
              int ck = c > 0 ? cumlen[c - 1] + k - 1 : k;
              double r = G.col(cj).dot(G.col(ck)) * (sds(cj) * sds(ck) * df);
              if (r * r >= params.clump_r2) {
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
              for (auto op : clumped)
                opp.push_back(pvals_per_chr[c].at(op).first);
              k = 0;
              for (auto oi : sortidx(opp)) {
                if (k == opp.size() - 1)
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
      } else {
        ld_clump_big(params.fileout + ".p" + std::to_string(i) + ".clump",
                     assocfiles[i], params.assoc_colnames, params.clump_bp,
                     params.clump_r2, params.clump_p1, params.clump_p2, snp,
                     data->G);
      }
    }
  }
}
