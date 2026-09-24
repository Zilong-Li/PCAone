/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/LD.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "LD.hpp"

#include <omp.h>
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

// scale each column to unit norm, so the correlation of two sites is the dot
// product of their columns. A site with no variance -- monomorphic, or fully
// explained by the PCs -- becomes a zero column: its R2 is 0, not NaN.
void normalize_ld_columns(Mat2D& G, Eigen::Index& nzero) {
  const double n1 = (double)G.rows() - 1.0;
  Eigen::Index nz = 0;
#pragma omp parallel for reduction(+ : nz)
  for (Eigen::Index j = 0; j < G.cols(); ++j) {
    const double ss = G.col(j).squaredNorm();
    if (std::sqrt(ss / n1) > VAR_TOL) {
      G.col(j) /= std::sqrt(ss);
    } else {
      G.col(j).setZero();
      ++nz;
    }
  }
  nzero += nz;
}

// one gzip member (RFC 1952). Concatenated members are a valid .gz file, as
// zcat, gzip -d, R and Python read it, so every thread compresses its own part
// of the .ld.gz output.
void gzip_member(const std::string& in, std::string& out) {
  z_stream s{};
  if (deflateInit2(&s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK)
    cao.error("failed to initialize zlib");
  out.resize(deflateBound(&s, in.size()));
  s.next_in = (Bytef*)in.data();
  s.avail_in = in.size();
  s.next_out = (Bytef*)&out[0];
  s.avail_out = out.size();
  if (deflate(&s, Z_FINISH) != Z_STREAM_END) cao.error("failed to compress the .ld.gz output");
  out.resize(s.total_out);
  deflateEnd(&s);
}

}  // namespace

// shared warning so every LD entry point reports this the same way
static void warn_zero_variance(Eigen::Index nzero, Eigen::Index ntotal) {
  if (nzero > 0)
    cao.warn(nzero, " of ", ntotal,
             " variants have no variance after ancestry adjustment (monomorphic, or fully explained by the "
             "PCs). their R2 is reported as 0 rather than NaN. consider --maf to remove them.");
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
// Why this equals the old G - USV': docs/ld-ancestry-adjusted.md
void adjust_for_pcs(Mat2D& G, const Mat2D& Q) {
  if (Q.size() == 0) return;
  G.noalias() -= Q * (Q.transpose() * G);
  G.rowwise() -= G.colwise().mean();
}

LDColumns::LDColumns(Data* data_, const Mat2D& Q_)
    : data(data_), Q(Q_), ooc(data_->params.out_of_core) {
  if (ooc) {
    rewind();
  } else {
    adjust_for_pcs(data->G, Q);
    normalize_ld_columns(data->G, nzero);
  }
}

void LDColumns::load(uint blk) {
  data->read_block_initial(data->start[blk], data->stop[blk], false);
  adjust_for_pcs(data->G, Q);
  normalize_ld_columns(data->G, nzero);
}

void LDColumns::rewind() {
  if (!ooc) return;
  nzero = 0;
  b = 0;
  data->check_file_offset_first_var();
  load(0);
}

uint LDColumns::max_reach(uint lo) const {
  if (!ooc) return data->nsnps - 1;
  const uint bl = std::min<uint>(lo / data->blocksize + 1, data->nblocks - 1);
  return data->stop[bl];
}

void LDColumns::need(uint lo, uint hi) {
  if (!ooc) return;
  while (hi > data->stop[b]) {
    prev.swap(data->G);  // keep block b; its buffer is reused for block b + 1
    load(++b);
  }
  if (lo < (b > 0 ? data->start[b - 1] : 0))
    cao.error("an LD window spans more than two blocks of -m. please increase -m or decrease --ld-bp");
}

// Batches of sites whose columns are multiplied at once: C = X_rows' X_cols is
// one matrix product (BLAS-3) instead of one dot product per pair, which reuses
// each column R times while it is in cache. Each batch buffer holds this many
// columns: about 64 MB, and with -m at most an eighth of a block, so the three
// buffers add under half a block to the two blocks that -m budgets for.
uint LDColumns::batch_cols() const {
  uint n = (1u << 23) / std::max<uint>(1, data->nsamples);
  if (ooc) n = std::min<uint>(n, data->blocksize / 8);
  return std::max<uint>(16, n);
}

Eigen::Ref<const Mat1D> LDColumns::col(uint k) const {
  if (!ooc) return data->G.col(k);
  if (k >= data->start[b]) return data->G.col(k - data->start[b]);
  return prev.col(k - data->start[b - 1]);
}

Eigen::Ref<const Mat2D> LDColumns::span(uint lo, uint hi, Mat2D& buf) const {
  const Eigen::Index n = hi - lo + 1;
  if (!ooc) return data->G.middleCols(lo, n);
  if (lo >= data->start[b]) return data->G.middleCols(lo - data->start[b], n);
  if (hi < data->start[b]) return prev.middleCols(lo - data->start[b - 1], n);
  // the two blocks in memory are separate matrices: copy the span across them
  const Eigen::Index n1 = data->start[b] - lo;
  buf.resize(data->nsamples, n);
  buf.leftCols(n1) = prev.middleCols(lo - data->start[b - 1], n1);
  buf.rightCols(n - n1) = data->G.leftCols(n - n1);
  return buf;
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

// Greedy pruning, in window order: of each pair above the cutoff, the site with
// the lower MAF is removed. The windows are taken in batches of R whose lead
// site is still kept, and the correlations of a batch with the kept sites in
// its span come from one matrix product. Processing the batch in order makes
// the same decisions as one window at a time. Leads removed within their own
// batch waste their row, so R grows while most leads survive and shrinks when
// they do not. F is filled as blocks are read (-m); every site used is read.
void ld_prune(LDColumns& X, const Mat1D& F, const String1D& variants, const SNPld& snp, double r2_tol,
              const std::string& fileout) {
  cao.print(tick.date(), "LD pruning, the site with the lower MAF of each pair is removed");
  const uint N = X.nsamples();
  const uint tile = X.batch_cols(), Rmax = std::min<uint>(1024, tile), Rmin = std::min<uint>(16, Rmax);
  ArrBool keep = ArrBool::Constant(variants.size(), true);
  Mat2D C, Xr, buf, Cg;
  std::vector<size_t> wins;
  std::vector<uint> kept;
  uint R = std::min<uint>(64, Rmax);
  const size_t nw = snp.ws.size();
  for (size_t w = 0; w < nw;) {
    // the next R windows whose lead site is still kept, all within reach of -m
    wins.clear();
    uint lo = 0, hi = 0, reach = 0;
    size_t w1 = w;
    for (; w1 < nw && wins.size() < R; ++w1) {
      const uint i = snp.ws[w1], e = i + snp.we[w1] - 1;
      if (!keep(i)) continue;
      if (wins.empty()) {
        lo = i;
        hi = e;
        reach = X.max_reach(i);
      }
      if (e > reach) {
        if (wins.empty())
          cao.error("an LD window spans more than two blocks of -m. please increase -m or decrease --ld-bp");
        break;
      }
      wins.push_back(w1);
      hi = std::max(hi, e);
    }
    w = w1;
    if (wins.empty()) continue;
    X.need(lo, hi);
    Xr.resize(N, wins.size());
    for (size_t t = 0; t < wins.size(); ++t) Xr.col(t) = X.col(snp.ws[wins[t]]);
    // C(t, k - lo) for the kept sites k in lo..hi. A tile that is mostly kept is
    // multiplied as it is; a sparse one is gathered first.
    C.resize(wins.size(), hi - lo + 1);
    for (uint c0 = lo; c0 <= hi; c0 += tile) {
      const uint c1 = std::min(hi, c0 + tile - 1);
      kept.clear();
      for (uint k = c0; k <= c1; ++k)
        if (keep(k)) kept.push_back(k);
      if (kept.empty()) continue;
      if (2 * kept.size() >= c1 - c0 + 1) {
        C.middleCols(c0 - lo, c1 - c0 + 1).noalias() = Xr.transpose() * X.span(c0, c1, buf);
      } else {
        buf.resize(N, kept.size());
        for (size_t t = 0; t < kept.size(); ++t) buf.col(t) = X.col(kept[t]);
        Cg.noalias() = Xr.transpose() * buf;
        for (size_t t = 0; t < kept.size(); ++t) C.col(kept[t] - lo) = Cg.col(t);
      }
    }
    size_t used = 0;
    for (size_t t = 0; t < wins.size(); ++t) {
      const uint i = snp.ws[wins[t]];
      if (!keep(i)) continue;  // removed by an earlier window of this batch
      ++used;
      for (int j = 1; j < snp.we[wins[t]]; ++j) {
        const uint k = i + j;
        if (!keep(k)) continue;
        const double r = C(t, k - lo);
        if (r * r > r2_tol) keep(MAF(F(k)) > MAF(F(i)) ? i : k) = false;
      }
    }
    if (used == wins.size())
      R = std::min(2 * R, Rmax);
    else if (2 * used < wins.size())
      R = std::max(R / 2, Rmin);
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
  // the columns of G have unit norm (LDColumns), so r is their dot product
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
        const double r = G.col(idx[j]).dot(G.col(idx[k]));
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

// All pairs within each window, written in window order. A batch of R
// consecutive windows is one matrix product over the span of their sites; the
// lines are then formatted and compressed by all threads, each into its own
// gzip member, and written in order.
void ld_r2(LDColumns& X, const String1D& variants, const SNPld& snp, const std::string& fileout, uint verbose) {
  const String1D bims = variant_labels(variants);
  const uint tile = X.batch_cols(), R = std::min<uint>(1024, tile);
  FILE* fp = fopen(fileout.c_str(), "wb");
  if (!fp) cao.error("can not open " + fileout);
  const int nthreads = omp_get_max_threads();
  std::vector<std::string> text(nthreads), gz(nthreads);
  gzip_member("CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2\n", gz[0]);
  fwrite(gz[0].data(), 1, gz[0].size(), fp);
  Mat2D C, bufr, bufc;
  const size_t nw = snp.ws.size();
  for (size_t w = 0; w < nw;) {
    const uint lo = snp.ws[w], reach = X.max_reach(lo);
    uint hi = lo;
    size_t w1 = w;
    for (; w1 < nw && w1 - w < R; ++w1) {
      const uint e = snp.ws[w1] + snp.we[w1] - 1;
      if (e > reach) break;
      hi = std::max(hi, e);
    }
    if (w1 == w) cao.error("an LD window spans more than two blocks of -m. please increase -m or decrease --ld-bp");
    if (verbose > 1) cao.print(tick.date(), "process windows", w, "to", w1 - 1);
    X.need(lo, hi);
    const uint rlast = snp.ws[w1 - 1];
    const auto Xr = X.span(lo, rlast, bufr);
    C.resize(rlast - lo + 1, hi - lo + 1);
    for (uint c0 = lo; c0 <= hi; c0 += tile) {
      const uint c1 = std::min(hi, c0 + tile - 1);
      C.middleCols(c0 - lo, c1 - c0 + 1).noalias() = Xr.transpose() * X.span(c0, c1, bufc);
    }
    const size_t nb = w1 - w;
    for (auto& g : gz) g.clear();
#pragma omp parallel num_threads(nthreads)
    {
      const int t = omp_get_thread_num(), T = omp_get_num_threads();
      std::string& s = text[t];
      s.clear();
      for (size_t x = w + nb * t / T; x < w + nb * (t + 1) / T; ++x) {
        const uint i = snp.ws[x];
        for (int j = 1; j < snp.we[x]; ++j) {
          const uint k = i + j;
          const double r = C(i - lo, k - lo);
          s += bims[i];
          s += '\t';
          s += bims[k];
          s += '\t';
          s += std::to_string(r * r);
          s += '\n';
        }
      }
      if (!s.empty()) gzip_member(s, gz[t]);
    }
    for (int t = 0; t < nthreads; ++t)
      if (!gz[t].empty() && fwrite(gz[t].data(), 1, gz[t].size(), fp) != gz[t].size())
        cao.error("failed to write data to ld.gz file");
    w = w1;
  }
  if (fclose(fp) != 0) cao.error("failed to write data to ld.gz file");
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
  LDColumns X(data, Q);  // in-core: every site, adjusted and normalized, now

  if (params.clump.empty()) {
    divide_pos_by_window(snp, params.ld_bp);
    if (params.print_r2) {
      ld_r2(X, variants, snp, params.fileout + ".ld.gz", params.verbose);
    } else {
      ld_prune(X, data->F, variants, snp, params.ld_r2, params.fileout);
    }
    warn_zero_variance(X.zero_variance(), data->nsnps);
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
      const std::string out = params.fileout + ".p" + std::to_string(i) + ".clump";
      if (params.out_of_core) {
        // gather the target sites; they are sorted, but may skip whole blocks
        if (i > 0) X.rewind();
        Mat2D G(data->nsamples, snp_t.pos.size());
        for (int sidx = 0, c = 0; c < (int)idx_per_chr.size(); c++) {
          for (auto& icol : idx_per_chr[c]) {
            X.need(icol, icol);
            G.col(sidx) = X.col(icol);
            icol = sidx++;
          }
        }
        ld_clump_single_pheno(out, head, params.clump_bp, params.clump_r2, params.clump_p1, params.clump_p2, G,
                              idx_per_chr, bp_per_chr, pvals_per_chr);
      } else {
        ld_clump_single_pheno(out, head, params.clump_bp, params.clump_r2, params.clump_p1, params.clump_p2, data->G,
                              idx_per_chr, bp_per_chr, pvals_per_chr);
      }
    }
    warn_zero_variance(X.zero_variance(), data->nsnps);
  }
}
