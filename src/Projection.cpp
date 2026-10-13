/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Projection.cpp
 * @author      Zilong Li
 * Copyright (C) 2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Projection.hpp"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <memory>

#include "Cmd.hpp"
#include "Common.hpp"
#include "Data.hpp"
#include "LinAlg.hpp"
#include "PackedScore.hpp"
#include "TextIO.hpp"
#include "Utils.hpp"
#include "Variants.hpp"

Mat2D solve_bootstrap_projection_no_missing(const Mat2D& design, const Mat2D& G, const std::vector<uint>& counts) {
  const int M = design.rows();
  const int K = design.cols();
  const int N = G.rows();
  Mat2D ata = Mat2D::Zero(K, K);
  Mat2D atg = Mat2D::Zero(K, N);
  int distinct = 0;
  for (int j = 0; j < M; ++j) {
    const uint w = counts[j];
    if (w == 0) continue;
    ++distinct;
    ata.noalias() += (double)w * design.row(j).transpose() * design.row(j);
    atg.noalias() += (double)w * design.row(j).transpose() * G.col(j).transpose();
  }
  // fewer distinct sites than PCs: no unique solution, see below
  if (distinct < K) return Mat2D::Constant(N, K, std::numeric_limits<double>::quiet_NaN());
  return ata.ldlt().solve(atg).transpose();
}

Mat2D solve_bootstrap_projection_missing(const Mat2D& design,
                                         const ArrBool& C,
                                         const Mat2D& G,
                                         const std::vector<uint>& counts) {
  const int M = design.rows();
  const int K = design.cols();
  const int N = G.rows();
  Mat2D U(N, K);
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    Mat2D ata = Mat2D::Zero(K, K);
    Mat1D atg = Mat1D::Zero(K);
    // the called sites drawn, each counted once: a site drawn twice adds no
    // rank, and counting the draws let a replicate of 2 sites drawn twice each
    // through for K = 4, with a singular system and arbitrary scores
    int observed = 0;
    for (int j = 0; j < M; ++j) {
      const uint w = counts[j];
      if (w == 0 || C((Eigen::Index)j * N + i)) continue;
      ata.noalias() += (double)w * design.row(j).transpose() * design.row(j);
      atg.noalias() += (double)w * design.row(j).transpose() * G(i, j);
      ++observed;
    }
    if (observed < K) {
      U.row(i).setConstant(std::numeric_limits<double>::quiet_NaN());
    } else {
      U.row(i) = ata.ldlt().solve(atg).transpose();
    }
  }
  return U;
}

namespace {

std::vector<std::pair<int, int>> pc_pairs(int K) {
  std::vector<std::pair<int, int>> pairs;
  pairs.reserve(static_cast<size_t>(K) * (K - 1) / 2);
  for (int a = 0; a < K; ++a) {
    for (int b = a + 1; b < K; ++b) pairs.emplace_back(a, b);
  }
  return pairs;
}

// "<i + 1>\t<FID>\t<IID>": the first columns of the row of sample i in the bootstrap files
void put_sample(TextBuf& buf, Eigen::Index i, const SampleIds& ids) {
  char* p = buf.reserve(22);
  p = put_uint(i + 1, p);
  *p++ = '\t';
  buf.commit(p);
  buf.append(ids.fid[i]);
  buf.append("\t", 1);
  buf.append(ids.iid[i]);
}

// "<i + 1>\t<FID>\t<IID>\tPC<a>" then the values, each after a tab
void put_row(TextBuf& buf,
             Eigen::Index i,
             const SampleIds& ids,
             std::initializer_list<int> pcs,
             std::initializer_list<double> vals) {
  put_sample(buf, i, ids);
  char* p = buf.reserve(pcs.size() * 24 + vals.size() * (TEXT_NUM_MAX + 1) + 1);
  for (int pc : pcs) std::memcpy(p, "\tPC", 3), p = put_uint(pc + 1, p + 3);
  for (double v : vals) *p++ = '\t', p = put_g6(v, p);
  *p++ = '\n';
  buf.commit(p);
}

const char* const BOOT_TSV_HEADER = "sample\tFID\tIID\tpc\tbaseline\tmean\tbootstrap_se\tmin\tmax\trmsd\n";
const char* const BOOT_COV_HEADER =
    "sample\tFID\tIID\tpc_x\tpc_y\tbaseline_x\tbaseline_y\tmean_x\tmean_y\tvar_x\tvar_y\tcov\tcorr\n";
const char* const BOOT_RAW_HEADER = "replicate\tsample\tFID\tIID";

// the replicate r (0-based) of the samples first, first + 1, ... as rows of
// .proj.bootstrap.eigvecs
void write_bootstrap_replicate(
    std::ostream& raw, uint r, const Eigen::Ref<const Mat2D>& Ub, Eigen::Index first, const SampleIds& ids) {
  const int K = Ub.cols();
  write_lines(raw, Ub.rows(), 48 + 12 * K, [&](Eigen::Index i, TextBuf& buf) {
    char* p = buf.reserve(22);
    p = put_uint(r + 1, p);
    *p++ = '\t';
    buf.commit(p);
    put_sample(buf, first + i, ids);
    p = buf.reserve(K * (TEXT_NUM_MAX + 1) + 1);
    for (int k = 0; k < K; ++k) *p++ = '\t', p = put_g6(Ub(i, k), p);
    *p++ = '\n';
    buf.commit(p);
  });
}

// The per-sample summaries of the bootstrap replicates of n samples, which
// write() appends to .proj.bootstrap.tsv and .proj.bootstrap.cov.tsv. The
// out-of-core projection keeps one per group of samples.
class BootstrapSummary {
 public:
  BootstrapSummary(Eigen::Index n, int K)
      : pairs(pc_pairs(K)),
        sum(Mat2D::Zero(n, K)),
        sumsq(Mat2D::Zero(n, K)),
        diff_sumsq(Mat2D::Zero(n, K)),
        minv(Mat2D::Constant(n, K, std::numeric_limits<double>::infinity())),
        maxv(Mat2D::Constant(n, K, -std::numeric_limits<double>::infinity())),
        cross_sum(Mat2D::Zero(n, pairs.size())) {}

  void add(const Eigen::Ref<const Mat2D>& Ub, const Eigen::Ref<const Mat2D>& U0) {
    // no sign alignment: the design V*S is fixed, so Ub is an ordinary
    // least-squares estimate with no arbitrary sign. Flipping it towards U0
    // would fold the replicates onto the baseline's side and bias the mean,
    // variance and covariance.
    sum += Ub;
    sumsq += Ub.array().square().matrix();
    diff_sumsq += (Ub - U0).array().square().matrix();
    minv = minv.cwiseMin(Ub);
    maxv = maxv.cwiseMax(Ub);
    for (int p = 0; p < (int)pairs.size(); ++p) {
      cross_sum.col(p).array() += Ub.col(pairs[p].first).array() * Ub.col(pairs[p].second).array();
    }
  }

  // the rows of these samples, the first of them being sample `first` (0-based)
  void write(std::ostream& out,
             std::ostream& covout,
             const Eigen::Ref<const Mat2D>& U0,
             uint nreps,
             Eigen::Index first,
             const SampleIds& ids) const {
    const Eigen::Index N = sum.rows();
    const int K = sum.cols();
    const double R = (double)nreps;
    Mat2D mean = sum / R;
    Mat2D variance = (sumsq.array() - (sum.array().square() / R)).max(0.0) / std::max(1.0, R - 1.0);
    Mat2D sd = variance.array().sqrt().matrix();
    Mat2D rmsd = (diff_sumsq.array() / R).sqrt().matrix();

    write_lines(out, N * K, 96, [&](Eigen::Index l, TextBuf& buf) {
      const Eigen::Index i = l / K, k = l % K;
      put_row(buf, first + i, ids, {(int)k}, {U0(i, k), mean(i, k), sd(i, k), minv(i, k), maxv(i, k), rmsd(i, k)});
    });
    const Eigen::Index np = (Eigen::Index)pairs.size();
    write_lines(covout, N * np, 160, [&](Eigen::Index l, TextBuf& buf) {
      const Eigen::Index i = l / np;
      const int a = pairs[l % np].first, b = pairs[l % np].second;
      const double cov = (cross_sum(i, l % np) - R * mean(i, a) * mean(i, b)) / std::max(1.0, R - 1.0);
      const double denom = std::sqrt(variance(i, a) * variance(i, b));
      const double corr = denom > 0.0 ? cov / denom : std::numeric_limits<double>::quiet_NaN();
      put_row(buf, first + i, ids, {a, b},
              {U0(i, a), U0(i, b), mean(i, a), mean(i, b), variance(i, a), variance(i, b), cov, corr});
    });
  }

 private:
  std::vector<std::pair<int, int>> pairs;
  Mat2D sum, sumsq, diff_sumsq, minv, maxv, cross_sum;
};

// the files of the bootstrap: .proj.bootstrap.tsv and .cov.tsv, each with its header
struct BootstrapFiles {
  std::string ftsv, fcov;
  std::ofstream out, covout;
  explicit BootstrapFiles(const Param& params)
      : ftsv(params.fileout + ".proj.bootstrap.tsv"),
        fcov(params.fileout + ".proj.bootstrap.cov.tsv") {
    out.open(ftsv, std::ios::binary);
    if (!out.is_open()) cao.error("can not open " + ftsv);
    covout.open(fcov, std::ios::binary);
    if (!covout.is_open()) cao.error("can not open " + fcov);
    out << BOOT_TSV_HEADER;
    covout << BOOT_COV_HEADER;
  }
  void close() {
    out.close();
    if (!out) cao.error("error writing " + ftsv);
    covout.close();
    if (!covout) cao.error("error writing " + fcov);
  }
};

// the IDs of the n target samples, or NA for each if they could not be read
SampleIds ids_or_na(SampleIds ids, Eigen::Index n) {
  if (ids.error.empty() && (Eigen::Index)ids.size() == n) return ids;
  cao.warn(
      "the FID and IID of the bootstrap files are NA:",
      ids.error.empty() ? std::to_string(ids.size()) + " sample IDs for " + std::to_string(n) + " samples" : ids.error);
  ids.fid.assign(n, "NA");
  ids.iid.assign(n, "NA");
  ids.error.clear();
  return ids;
}

void open_bootstrap_raw(std::ofstream& raw, const std::string& path, int K) {
  raw.open(path, std::ios::binary);
  if (!raw.is_open()) cao.error("can not open " + path);
  raw << BOOT_RAW_HEADER;
  for (int k = 0; k < K; ++k) raw << "\tPC" << k + 1;
  raw << '\n';
}

}  // namespace

void write_projection_bootstrap_stats(const Mat2D& design,
                                      const ArrBool& C,
                                      const Mat2D& G,
                                      const Mat2D& U0,
                                      uint nreps,
                                      const Param& params,
                                      const SampleIds& ids) {
  const int M = design.rows();
  const int N = G.rows();
  const int K = design.cols();
  if (nreps < 2) cao.error("--project-bootstrap must be at least 2");
  if (M < K) cao.error("projection bootstrap requires at least as many matched SNPs as PCs");

  cao.print(tick.date(), "run projection SNP bootstrap with", nreps, "replicates");
  const bool has_missing = C.size() && C.count() > 0;
  PortableRng rng(params.seed);  // uniform_int_distribution differs between standard libraries

  BootstrapSummary summary(N, K);
  std::ofstream raw;
  if (params.project_bootstrap_save) open_bootstrap_raw(raw, params.fileout + ".proj.bootstrap.eigvecs", K);

  std::vector<uint> counts(M);
  for (uint r = 0; r < nreps; ++r) {
    std::fill(counts.begin(), counts.end(), 0);
    for (int draw = 0; draw < M; ++draw) counts[rng.below(M)]++;

    Mat2D Ub = has_missing ? solve_bootstrap_projection_missing(design, C, G, counts)
                           : solve_bootstrap_projection_no_missing(design, G, counts);
    summary.add(Ub, U0);
    if (raw.is_open()) write_bootstrap_replicate(raw, r, Ub, 0, ids);
  }
  if (raw.is_open()) {
    raw.close();
    if (!raw) cao.error("error writing " + params.fileout + ".proj.bootstrap.eigvecs");
  }

  BootstrapFiles files(params);
  summary.write(files.out, files.covout, U0, nreps, 0, ids);
  files.close();
  cao.print(tick.date(), "projection bootstrap diagnostics saved to", params.fileout + ".proj.bootstrap.tsv");
}

void solve_projection_scores(const Mat2D& V, const ArrBool& C, const Mat2D& G, Mat2D& U) {
  if (U.rows() == 0 || U.cols() == 0) return;
  double p_miss = C.size() ? (double)C.count() / (double)C.size() : 0.0;
  if (p_miss == 0.0) {
    Eigen::ColPivHouseholderQR<Mat2D> qr(V);
#pragma omp parallel for
    for (uint i = 0; i < (uint)U.rows(); i++) {
      U.row(i) = qr.solve(G.row(i).transpose());
    }
  } else {
    // with fewer called sites than PCs the least-squares system has no unique
    // solution; the pivoted QR returned 0 for every PC (a sample with no calls
    // at all was projected to the origin), indistinguishable from a real score
    uint64 nunder = 0;
#pragma omp parallel for reduction(+ : nunder)
    for (uint i = 0; i < (uint)U.rows(); i++) {
      Int1D idx;
      for (int j = 0; j < V.rows(); j++) {
        if (!C((Eigen::Index)j * U.rows() + i)) idx.push_back(j);
      }
      if ((Eigen::Index)idx.size() < V.cols()) {
        U.row(i).setConstant(std::numeric_limits<double>::quiet_NaN());
        ++nunder;
        continue;
      }
      U.row(i) = V(idx, Eigen::all).colPivHouseholderQr().solve(G(i, idx).transpose());
    }
    if (nunder > 0)
      cao.warn(std::to_string(nunder) + " sample(s) have fewer called sites than PCs and are written as NA");
  }
}

// Per-site factor a_j for --project 3, where X_ref = a_j * (g/2 - f) is what the
// reference PCA decomposed in terms of PCAone's 0..1 centred coding. The EM
// regresses a_j * (E[g]/2 - f) on V*S, so U lands on the reference's scale, and
// maps the reconstruction back with pi = f + (U*S*V')_j / a_j. PCAngsd's
// Algorithm 1 is the a_j = 2 case: pi = f + recon/2 on centred dosages.
//
//   standardized (scale=-9)      sqrt(ploidy)/sd   as standardize_E_ref()
//   centred dosages (gscale=2)   2                 pcangsd: dosage - 2f
//   centred 0..1                 1                 -D/--ld, --scale 0
static Mat1D gl_site_scale(const Mat1D& F, bool standardize, const UsvTransform& t) {
  Mat1D a = Mat1D::Constant(F.size(), (double)t.gscale);
  if (!standardize) return a;
  const double rploidy = sqrt((double)t.ploidy);
  for (Eigen::Index j = 0; j < F.size(); ++j) {
    const double sd = sqrt(F(j) * (1.0 - F(j)));
    if (sd > VAR_TOL) a(j) = rploidy / sd;
  }
  return a;
}

namespace {

using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// the K x K symmetric matrix packed by rows of its upper triangle, (a, b) with a <= b
void unpack_sym(const double* p, int K, Mat2D& A) {
  for (int a = 0; a < K; ++a)
    for (int b = a; b < K; ++b) A(a, b) = A(b, a) = *p++;
}

// least squares from its normal equations A z = b, with A packed; u = z / S
template <typename Out>
void solve_packed(const double* a, const double* b, const Mat1D& S, Mat2D& A, Eigen::LDLT<Mat2D>& ldlt, Out&& u) {
  const int K = S.size();
  unpack_sym(a, K, A);
  ldlt.compute(A);
  u = ldlt.solve(Eigen::Map<const Mat1D>(b, K));
  u.array() /= S.array();
}

// out += X * Y, for a block X of the target (samples x sites) and a thin Y, in
// row panels of about 8 MB of X, one thread each. Eigen's GEMM (and so
// mul_X_Y) packs up to half the rows of such an X at once when Y is thin: a
// copy of half the block, outside -m. The panels depend on X.cols() only, so
// the sums do not depend on the number of threads.
void add_X_Y_panels(const Eigen::Ref<const Mat2D>& X, const Eigen::Ref<const Mat2D>& Y, Eigen::Ref<Mat2D> out) {
  const Eigen::Index n = X.rows();
  const Eigen::Index h =
      std::max<Eigen::Index>(48, (Eigen::Index(1) << 20) / std::max<Eigen::Index>(1, X.cols()) / 48 * 48);
  const Eigen::Index np = (n + h - 1) / h;
#pragma omp parallel for schedule(static)
  for (Eigen::Index p = 0; p < np; ++p) {
    const Eigen::Index r = p * h, w = std::min(h, n - r);
    out.middleRows(r, w).noalias() += X.middleRows(r, w) * Y;
  }
}


// the upper triangle of the K x K row-major sq, by rows, added to out
inline void add_packed(const double* sq, int K, double* out) {
  for (int a = 0, q = 0; a < K; ++a)
    for (int c = a; c < K; ++c) out[q++] += sq[a * K + c];
}

// One pass of the out-of-core projection over the 2-bit calls (packed_calls())
// of the samples [g0, g0 + gn): sums(i - g0, :) += the calls of sample i,
// centred and scaled as the reference was, times the design Vd. That is the
// product of the decoded blocks with Vd. For --project 2 (A0 not empty), A_i
// and the missing calls as from the decoded blocks: A_i is the block's sum of
// v_j v_j' less the sites the sample misses, or the sum over its called sites
// when it misses most of the block. Returns the missing calls of the group.
uint64 project_packed_pass(Data* data,
                           const Mat2D& Vd,
                           bool standardize,
                           const UsvTransform& usv,
                           Eigen::Index g0,
                           Eigen::Index gn,
                           RowMat& sums,
                           RowMat& A0,
                           std::vector<uint64>& nmiss) {
  const int K = Vd.cols();
  const bool lsq = A0.size() > 0;
  const Eigen::Index g1 = g0 + gn, kp = sums.cols();
  const PackedCoding coding = data->packed_coding();
  const double rploidy = std::sqrt((double)usv.ploidy);
  const uint64 bytes = ((uint64)data->nsamples + 3) / 4;
  const int nthreads = std::max(1, omp_get_max_threads());
  // chunks of whole words of the calls (32 samples), up to 4096 samples so
  // that their sums and calls stay in the cache. A chunk tabulates each group
  // of 4 sites (packed_product_add), so with fewer chunks than threads the
  // sites are split in parts as well, each part of a chunk with its own sums.
  const Eigen::Index csz =
      std::min<Eigen::Index>(4096, std::max<Eigen::Index>(1024, gn / (4 * nthreads))) / 32 * 32;
  const Eigen::Index nchunk = (gn + csz - 1) / csz;
  const int nparts = (int)std::max<Eigen::Index>(1, (nthreads + nchunk - 1) / nchunk);
  std::vector<RowMat> part_sums(nparts - 1, RowMat::Zero(gn, kp));
  std::vector<double> value;
  RowMat Vr;            // the rows of V of the block, padded to kp
  Mat1D Asum;           // their sum of v_j v_j', packed
  uint64 nmiss_group = 0;
  using clock = std::chrono::steady_clock;
  double t_read = 0, t_work = 0, t_miss = 0, t_prod = 0;
  auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };

  for (uint bi = 0; bi < data->nblocks; ++bi) {
    const uint64 start = data->start[bi];
    const Eigen::Index cols = data->stop[bi] - start + 1;
    auto t0 = clock::now();
    const PCAone::PackedBlock blk{data->read_block_packed(start, data->stop[bi]), bytes, cols, coding.missing};
    t_read += since(t0);
    t0 = clock::now();
    // the value of each code at each site: the frequency of the code less F,
    // scaled as standardize_block_ref(), and 0 for a missing call
    value.resize(4 * cols);
    for (Eigen::Index j = 0; j < cols; ++j) {
      const double f = data->F(start + j), sd = std::sqrt(f * (1.0 - f));
      const double scale = standardize && sd > VAR_TOL ? rploidy / sd : 1.0;
      for (int c = 0; c < 4; ++c) value[4 * j + c] = c == coding.missing ? 0.0 : scale * (coding.value[c] - f);
    }
    const auto Vb = Vd.middleRows(start, cols);
    if (lsq) {
      Vr = RowMat::Zero(cols, kp);  // padded as add_gram() takes it
      Vr.leftCols(K) = Vb;
      const Mat2D VtV = Vb.transpose() * Vb;
      Asum = Mat1D::Zero(K * (K + 1) / 2);
      add_packed(Mat2D(VtV.transpose()).data(), K, Asum.data());  // row-major from column-major
    }
    const Eigen::Index ngroups = (cols + 3) / 4;

    uint64 nmiss_block = 0;
#pragma omp parallel for collapse(2) schedule(dynamic) reduction(+ : nmiss_block, t_miss, t_prod)
    for (Eigen::Index ch = 0; ch < nchunk; ++ch) {
      for (int p = 0; p < nparts; ++p) {
        double tm = omp_get_wtime();
        const Eigen::Index i0 = g0 + ch * csz, i1 = std::min(g1, i0 + csz);
        const uint64 w0 = (uint64)i0 / 32, w1 = ((uint64)i1 + 31) / 32;
        // the missing calls of the chunk, by its first part, site by site
        // through the calls, which are then in the cache for the sums
        if (p == 0 && !lsq) {
          for (Eigen::Index j = 0; j < cols; ++j)
            for (uint64 w = w0; w < w1; ++w)
              nmiss_block += __builtin_popcountll(PCAone::missing_word(blk, j, w, i1));
        } else if (p == 0) {
          // the missing calls (sample, site) in one pass, then their sites by
          // sample, in site order (a counting sort)
          static thread_local std::vector<uint64> miss;  // sample << 32 | site
          static thread_local std::vector<uint32_t> sites;
          miss.clear();
          for (Eigen::Index j = 0; j < cols; ++j)
            for (uint64 w = w0; w < w1; ++w)
              for (uint64 m = PCAone::missing_word(blk, j, w, i1); m; m &= m - 1)
                miss.push_back((w * 32 + (__builtin_ctzll(m) >> 1) - i0) << 32 | (uint64)j);
          std::vector<uint32_t> mb(i1 - i0 + 1, 0);
          for (uint64 e : miss) ++mb[(e >> 32) + 1];
          for (Eigen::Index i = 0; i < i1 - i0; ++i) mb[i + 1] += mb[i];  // the first of each sample
          sites.resize(miss.size());
          {
            std::vector<uint32_t> at(mb.begin(), mb.end() - 1);
            for (uint64 e : miss) sites[at[e >> 32]++] = (uint32_t)e;
          }
          std::vector<double> sq((size_t)kp * kp);
          std::vector<uint32_t> called;
          for (Eigen::Index i = 0; i < i1 - i0; ++i) {
            const uint32_t n = mb[i + 1] - mb[i];
            const Eigen::Index li = i0 + i - g0;
            nmiss[li] += n, nmiss_block += n;
            std::fill(sq.begin(), sq.end(), 0.0);
            if (2 * (Eigen::Index)n <= cols) {
              for (Eigen::Index q = 0; q < Asum.size(); ++q) A0(li, q) += Asum(q);
              if (n == 0) continue;
              PCAone::add_gram(Vr.data(), &sites[mb[i]], n, kp, -1.0, sq.data());  // less the missing sites
            } else {
              // most of the block missing: the sum over the called sites
              called.clear();
              for (Eigen::Index j = 0; j < cols; ++j)
                if (PCAone::packed_code(blk, j, i0 + i) != coding.missing) called.push_back((uint32_t)j);
              PCAone::add_gram(Vr.data(), called.data(), called.size(), kp, 1.0, sq.data());
            }
            for (int a = 0, q = 0; a < K; ++a)
              for (int c = a; c < K; ++c) A0(li, q++) += sq[a * kp + c];
          }
        }
        const Eigen::Index j0 = 4 * (ngroups * p / nparts), j1 = std::min(cols, 4 * (ngroups * (p + 1) / nparts));
        RowMat& out = p == 0 ? sums : part_sums[p - 1];
        t_miss += omp_get_wtime() - tm;
        tm = omp_get_wtime();
        PCAone::packed_product_add(blk, j0, j1, value.data(), Vb, i0, i1, &out(i0 - g0, 0), kp);
        t_prod += omp_get_wtime() - tm;
      }
    }
    nmiss_group += nmiss_block;
    t_work += since(t0);
  }
  for (const auto& ps : part_sums) sums += ps;
  if (data->params.verbose > 1)
    cao.print(tick.date(), "out-of-core projection: waited", t_read, "s for the calls, worked", t_work,
              "s; thread seconds: missing calls", t_miss, ", sums", t_prod);
  return nmiss_group;
}

}  // namespace

/**
 * Out-of-core --project 1 and 2 (-m), on PLINK and PGEN input: one pass over
 * the matched sites, block by block, holding sums per sample instead of the
 * N x M genotypes, so a target of any size is projected in one run.
 *
 * --project 1 is the product u = g V S^-1, summed over the blocks.
 *
 * --project 2 solves the least squares (V S) u = g over the called sites of
 * each sample as V z = g, u = S^-1 z, from its normal equations
 *     A_i = sum over the called sites j of v_j v_j',  b_i = sum over j of g_ij v_j
 * (v_j the row of V of site j), and z_i = A_i^-1 b_i. b is a product of the
 * block with V. A_i, packed to K(K+1)/2 values, is the sum over all the sites
 * of the block less those the sample misses, or the sum over the called ones
 * when it misses most of the block, so the work is in the fewer of the two. In
 * exact arithmetic this is the in-core solution, which takes the pivoted QR of
 * the called rows of V S; the two agree to rounding.
 *
 * --project-bootstrap does the same with the sites weighted by the number of
 * times each replicate drew them, all replicates in the same pass. The draws
 * are those of the in-core bootstrap, from --seed.
 *
 * When the sums of all samples do not fit in -m with a block, the samples are
 * taken in groups, one pass over the genotypes each.
 */
static Mat2D run_projection_blockwise(Data* data,
                                      const Param& params,
                                      const Mat2D& V,
                                      const Mat1D& S,
                                      const UsvTransform& usv,
                                      bool standardize,
                                      const SampleIds& ids) {
  const Eigen::Index N = data->nsamples, M = data->nsnps;
  const int K = V.cols(), KK = K * (K + 1) / 2;
  const bool lsq = params.project == 2;
  const uint R = lsq ? params.project_bootstrap : 0;
  if (V.rows() != M) cao.error("the loadings have", V.rows(), "rows for", M, "matched sites");
  if (R > 0 && M < K) cao.error("projection bootstrap requires at least as many matched SNPs as PCs");

  // the 2-bit calls of the block straight into the sums (PackedScore.hpp),
  // instead of decoding the block to doubles and multiplying: PLINK, and PGEN
  // hard calls. The bootstrap's sums per replicate still take the product.
  const bool packed = R == 0 && data->packed_calls();
  const int kp = (K + 3) / 4 * 4;  // the sums of the 2-bit calls, padded as their tables
  const uint64 bytes = ((uint64)N + 3) / 4;

  // --- memory, in doubles of -m: what is held whatever the blocks are, per
  // sample of a group, and per site of a block
  const double D = 134217728.0;  // doubles per GB, as Data::prepare()
  const double budget = params.memory * D;
  const double npairs = K * (K - 1) / 2.0;
  // U, V and its copy, F with the lookup table of the reader, the draws, and
  // the panels packed by each thread in add_X_Y_panels() (the tables of the
  // 2-bit calls, a quarter of that); --project 1 sums the 2-bit calls apart from U
  const int nthreads = std::max(1, omp_get_max_threads());
  const double fixed = (double)N * K + 2.0 * M * K + 5.0 * M + (R ? (double)M * R / 8.0 : 0.0) +
                       nthreads * (double)(packed ? 1 << 18 : 1 << 20) + (packed && !lsq ? (double)N * kp : 0.0);
  // the count of missing calls, A and b, then for each replicate A, b (then u) and
  // the number of called sites drawn, and the summaries of BootstrapSummary
  const double per_sample =
      1.0 + (lsq ? KK + (packed ? kp : K) : 0) + (R ? R * (KK + K + 1.0) + 5.0 * K + npairs : 0);
  // the block (N rows) and the raw records read ahead, W, V * draws; or the
  // 2-bit calls (twice, as read ahead), W and the values of the 4 codes
  auto per_site = [&](double ng) {
    if (packed) return N / 16.0 + (lsq ? kp : 0) + 4.0;
    return 1.0625 * N + ng / 8.0 + (lsq ? KK : 0) + (double)R * K + K;
  };
  const double avail = budget - fixed;
  double ng = (double)N;  // --project 1 keeps no sums besides U: all samples in one pass
  if (lsq && per_sample * N > avail / 2) ng = std::floor(std::max(0.0, avail / 2) / per_sample);
  if (packed && ng < N && ng >= 32) ng = std::floor(ng / 32) * 32;  // groups start at a word of the calls
  double b = std::floor((avail - per_sample * ng) / per_site(ng));
  // blocks of up to 128 MB of 2-bit calls, so that reading the next one
  // overlaps this one; for --project 2 up to twice the A_i of the group, which
  // each block reads and writes once
  if (packed) b = std::min(b, std::max(64.0, std::floor(std::max(134217728.0, lsq ? 16.0 * KK * ng : 0.0) / bytes)));
  const double ng_min = std::min<double>(N, 256), b_min = std::min<double>(M, 64);
  if (ng < ng_min || b < b_min) {
    ng = std::max(ng, ng_min);
    b = std::max(b, b_min);
    cao.warn("-m is too small for the out-of-core projection, which needs about",
             (fixed + per_sample * ng + per_site(ng) * b) / D, "GB with groups of", (uint64)ng, "samples and",
             (uint64)b, "sites per block");
  }
  const Eigen::Index gmax = (Eigen::Index)std::min<double>(ng, N);
  data->set_blocks((uint)std::min<double>(b, M));
  const Eigen::Index npass = (N + gmax - 1) / gmax;
  cao.print(tick.date(), "out-of-core projection by -m/--memory: blocksize =", data->blocksize,
            ", nblocks =", data->nblocks, ", samples per pass =", gmax, ", passes =", npass,
            packed ? ", from the 2-bit calls" : "");

  // the design: --project 1 multiplies by V S^-1, --project 2 solves V z = g
  const Mat2D Vd = lsq ? V : Mat2D(V * S.cwiseInverse().asDiagonal());
  // the draws of the bootstrap, as write_projection_bootstrap_stats() makes
  // them, site-major: counts[j * R + r] is the weight of site j in replicate r.
  // A byte each: M draws over M sites put 256 on one site with probability
  // about M / 256!, i.e. never, but it is checked.
  std::vector<uint8_t> counts;
  if (R > 0) {
    cao.print(tick.date(), "run projection SNP bootstrap with", R, "replicates");
    counts.assign((size_t)M * R, 0);
    PortableRng rng(params.seed);
    for (uint r = 0; r < R; ++r)
      for (Eigen::Index draw = 0; draw < M; ++draw)
        if (++counts[rng.below(M) * R + r] == 0) cao.error("a site was drawn 256 times in one bootstrap replicate");
  }

  data->mark_missing = true;  // a missing call is read as NaN
  Mat2D U = Mat2D::Zero(N, K);
  RowMat Usum, bsum;  // the sums of the 2-bit calls: --project 1 of all samples, b of the group
  if (packed && !lsq) Usum = RowMat::Zero(N, kp);
  uint64 nmiss_all = 0, nunder = 0;

  std::unique_ptr<BootstrapFiles> bfiles;
  std::ofstream raw;
  std::vector<std::string> raw_parts;  // with more than one pass, the replicates of each, merged at the end
  if (R > 0) bfiles = std::make_unique<BootstrapFiles>(params);
  const std::string fraw = params.fileout + ".proj.bootstrap.eigvecs";

  for (Eigen::Index pass = 0; pass < npass; ++pass) {
    const Eigen::Index g0 = pass * gmax, gn = std::min(gmax, N - g0);
    if (npass > 1) cao.print(tick.date(), "projection pass", pass + 1, "of", npass, ": samples", g0 + 1, "to", g0 + gn);
    // the sums of the group
    RowMat A0, Ab, od;
    Mat2D b0, bb;
    if (lsq) {
      A0 = RowMat::Zero(gn, KK);
      if (packed)
        bsum = RowMat::Zero(gn, kp);
      else
        b0 = Mat2D::Zero(gn, K);
    }
    if (R > 0) {
      Ab = RowMat::Zero(gn, (Eigen::Index)R * KK);
      bb = Mat2D::Zero(gn, (Eigen::Index)R * K);
      od = RowMat::Zero(gn, R);
    }
    std::vector<uint64> nmiss(gn, 0);
    std::vector<uint32_t> mb(gn);  // missing calls of the sample in the block
    std::vector<uint8_t> sub(gn);  // A_i from the block's sum less the missing sites
    std::vector<uint8_t> miss;     // the block's missing calls of the group, column-major
    RowMat W, WcSum;               // v_j v_j' packed, per site; their sums per replicate
    Eigen::Matrix<double, 1, Eigen::Dynamic> Wsum, dsum;
    Mat2D Vc;  // V weighted by the draws, [replicate][PC] columns
    // chunks of samples, each worked by one thread: they own the rows of A
    const Eigen::Index nchunk = std::min<Eigen::Index>(gn, 8 * nthreads), csz = (gn + nchunk - 1) / nchunk;

    data->check_file_offset_first_var();
    if (packed) nmiss_all += project_packed_pass(data, Vd, standardize, usv, g0, gn, lsq ? bsum : Usum, A0, nmiss);
    for (uint bi = 0; bi < (packed ? 0u : data->nblocks); ++bi) {
      const uint64 start = data->start[bi];
      const Eigen::Index cols = data->stop[bi] - start + 1;
      data->read_block_initial(start, data->stop[bi], false);
      if (standardize) data->standardize_block_ref(usv, start, cols);  // NaN stays NaN
      auto Gg = data->G.block(g0, 0, gn, cols);
      const auto Vb = Vd.middleRows(start, cols);
      miss.resize((size_t)gn * cols);
      if (lsq) {
        W.resize(cols, KK);
        for (Eigen::Index j = 0; j < cols; ++j)
          for (int a = 0, q = 0; a < K; ++a)
            for (int c = a; c < K; ++c) W(j, q++) = Vb(j, a) * Vb(j, c);
        Wsum = W.colwise().sum();
      }
      const uint8_t* cb = R ? counts.data() + start * R : nullptr;
      if (R > 0) {
        Vc.resize(cols, (Eigen::Index)R * K);
        for (uint r = 0; r < R; ++r)
          for (int k = 0; k < K; ++k)
            for (Eigen::Index j = 0; j < cols; ++j) Vc(j, r * K + k) = cb[j * R + r] * Vb(j, k);
        WcSum = RowMat::Zero(R, KK);
        dsum = Eigen::Matrix<double, 1, Eigen::Dynamic>::Zero(R);
        for (Eigen::Index j = 0; j < cols; ++j)
          for (uint r = 0; r < R; ++r)
            if (cb[j * R + r]) {
              WcSum.row(r) += cb[j * R + r] * W.row(j);
              dsum(r) += 1;
            }
      }

      uint64 nmiss_block = 0;
#pragma omp parallel for schedule(static) reduction(+ : nmiss_block)
      for (Eigen::Index ch = 0; ch < nchunk; ++ch) {
        const Eigen::Index i0 = ch * csz, i1 = std::min(gn, i0 + csz);
        if (i0 >= i1) continue;
        // the missing calls: NaN to 0, the mean, and marked
        for (Eigen::Index j = 0; j < cols; ++j) {
          double* g = &Gg(0, j);
          uint8_t* m = &miss[(size_t)j * gn];
          for (Eigen::Index i = i0; i < i1; ++i) {
            m[i] = std::isnan(g[i]);
            if (m[i]) g[i] = 0.0, ++mb[i];
          }
        }
        if (lsq) {
          for (Eigen::Index i = i0; i < i1; ++i) {
            sub[i] = 2 * (Eigen::Index)mb[i] <= cols;
            if (!sub[i]) continue;
            A0.row(i) += Wsum;
            if (R > 0) {
              Ab.row(i) += Eigen::Map<const Eigen::Matrix<double, 1, Eigen::Dynamic>>(WcSum.data(), WcSum.size());
              od.row(i) += dsum;
            }
          }
          // subtract the missing sites, or add the called ones
          for (Eigen::Index j = 0; j < cols; ++j) {
            const uint8_t* m = &miss[(size_t)j * gn];
            for (Eigen::Index i = i0; i < i1; ++i) {
              if (m[i] != sub[i]) continue;
              const double sign = sub[i] ? -1.0 : 1.0;
              A0.row(i) += sign * W.row(j);
              for (uint r = 0; r < R; ++r) {
                const double c = cb[j * R + r];
                if (!c) continue;
                Ab.row(i).segment(r * KK, KK) += (sign * c) * W.row(j);
                od(i, r) += sign;
              }
            }
          }
        }
        for (Eigen::Index i = i0; i < i1; ++i) nmiss[i] += mb[i], nmiss_block += mb[i], mb[i] = 0;
      }
      nmiss_all += nmiss_block;

      auto Ug = U.middleRows(g0, gn);
      add_X_Y_panels(Gg, Vb, lsq ? Eigen::Ref<Mat2D>(b0) : Eigen::Ref<Mat2D>(Ug));
      if (R > 0) add_X_Y_panels(Gg, Vc, bb);
    }

    if (!lsq) continue;
    // solve each sample of the group, then each of its replicates
    uint64 under = 0;
#pragma omp parallel reduction(+ : under)
    {
      Mat2D A(K, K);
      Eigen::LDLT<Mat2D> ldlt(K);
      Mat1D z(K);
#pragma omp for schedule(static)
      for (Eigen::Index i = 0; i < gn; ++i) {
        // with fewer called sites than PCs there is no unique solution
        if ((uint64)M - nmiss[i] < (uint64)K) {
          U.row(g0 + i).setConstant(std::numeric_limits<double>::quiet_NaN());
          ++under;
        } else {
          const Mat1D bi = packed ? Mat1D(Eigen::Map<const Mat1D>(&bsum(i, 0), K)) : Mat1D(b0.row(i).transpose());
          solve_packed(&A0(i, 0), bi.data(), S, A, ldlt, z);
          U.row(g0 + i) = z.transpose();
        }
        for (uint r = 0; r < R; ++r) {
          auto ub = bb.row(i).segment(r * K, K);
          if (od(i, r) < K) {  // as solve_bootstrap_projection_missing()
            ub.setConstant(std::numeric_limits<double>::quiet_NaN());
          } else {
            const Mat1D rhs = ub.transpose();
            solve_packed(&Ab(i, r * KK), rhs.data(), S, A, ldlt, z);
            ub = z.transpose();
          }
        }
      }
    }
    nunder += under;

    if (R == 0) continue;
    const auto U0 = U.middleRows(g0, gn);
    BootstrapSummary summary(gn, K);
    for (uint r = 0; r < R; ++r) summary.add(bb.middleCols(r * K, K), U0);
    summary.write(bfiles->out, bfiles->covout, U0, R, g0, ids);
    if (params.project_bootstrap_save) {
      std::string path = fraw;
      if (npass > 1) {
        path = fraw + ".part" + std::to_string(pass);
        raw_parts.push_back(path);
        raw.open(path, std::ios::binary);
        if (!raw.is_open()) cao.error("can not open " + path);
      } else {
        open_bootstrap_raw(raw, path, K);
      }
      for (uint r = 0; r < R; ++r) write_bootstrap_replicate(raw, r, bb.middleCols(r * K, K), g0, ids);
      raw.close();
      if (!raw) cao.error("error writing " + path);
    }
  }

  if (packed && !lsq) U = Usum.leftCols(K);
  const double p_miss = (double)nmiss_all / ((double)N * M);
  cao.print(tick.date(), "the proportion of missingness  is", p_miss);
  if (!lsq && p_miss > 0) cao.warn("there are missing genotypes. recommend using --project 2 or 3.");
  if (lsq && p_miss == 0.0) cao.warn("there is no missing genotypes");
  if (nunder > 0)
    cao.warn(std::to_string(nunder) + " sample(s) have fewer called sites than PCs and are written as NA");

  if (R > 0) {
    bfiles->close();
    if (!raw_parts.empty()) {
      // the replicates in the order of the in-core run: replicate by replicate,
      // each with all the samples, so the parts of the passes are interleaved
      open_bootstrap_raw(raw, fraw, K);
      std::vector<std::ifstream> parts(raw_parts.size());
      for (size_t p = 0; p < parts.size(); ++p) {
        parts[p].open(raw_parts[p], std::ios::binary);
        if (!parts[p].is_open()) cao.error("can not open " + raw_parts[p]);
      }
      std::string line;
      for (uint r = 0; r < R; ++r)
        for (size_t p = 0; p < parts.size(); ++p) {
          const Eigen::Index gn = std::min(gmax, N - (Eigen::Index)p * gmax);
          for (Eigen::Index i = 0; i < gn; ++i) {
            if (!std::getline(parts[p], line)) cao.error("can not read " + raw_parts[p]);
            raw << line << '\n';
          }
        }
      raw.close();
      if (!raw) cao.error("error writing " + fraw);
      for (auto& part : parts) part.close();
      for (const auto& path : raw_parts) std::remove(path.c_str());
    }
    cao.print(tick.date(), "projection bootstrap diagnostics saved to", params.fileout + ".proj.bootstrap.tsv");
  }
  return U;
}

/**
 * options:
 * 1: simple, assume no missingness
 * 2: like smartPCA, solving g=Vx, can take missing genotypes
 * 3: iterative GL-aware projection (EM): alternates between updating individual allele frequencies
 *    from current PC scores and re-solving for PC scores with updated expected genotypes (BEAGLE only)
   // out-of-core (-m): --project 1 and 2 only, see run_projection_blockwise()
 */
void run_projection(Data* data, const Param& params) {
  const BimMatch match = match_variants_to_mbim(params);
  if (!match.identical) {
    data->keepSNPs = match.bim_indices;
    data->keepRefSNPs = match.mbim_indices;
    for (int k = 0; k < (int)match.flip.size(); ++k) {
      if (match.flip[k]) data->flipSNPs.push_back(k);
    }
    cao.warn("SNP info is not fully identical between input and reference .mbim");
    cao.print(tick.date(), "projection will use", match.bim_indices.size(), " overlapped sites. there are",
              data->flipSNPs.size(), " flipped alleles");
    if (!data->flipSNPs.empty())
      cao.warn(data->flipSNPs.size(), " SNPs have flipped ref/alt alleles and will be corrected");
  }
  cao.print(tick.date(), "run projection");
  // the FID and IID of the target (.fam, .psam or the BEAGLE header), for
  // .eigvecs2 and the first columns of the bootstrap files
  const SampleIds ids = read_sample_ids(params);
  const SampleIds boot_ids = params.project_bootstrap > 0 ? ids_or_na(ids, data->nsamples) : SampleIds{};
  cao.print(tick.date(), "start parsing V:", params.fileV, ", S:", params.fileS);
  uint nsamples, nsnps;
  Mat1D S;
  // read the reference's transform before touching G: the target has to be put
  // on the same scale as the matrix that produced V and S, which is not
  // necessarily what this run's -C/--scale says. --project 3 also accepts a
  // dosage-scale (pcangsd) reference, and applies the transform per site
  // inside its EM (gl_site_scale) rather than to G up front.
  UsvTransform usv;
  read_sigvals(params.fileS, nsamples, nsnps, S, &usv);
  if (!params.out_of_core && params.project != 3) {
    // a double and a missing flag (ArrBool) per genotype
    const double gb = 9.0 * data->nsamples * (match.identical ? data->nsnps : match.bim_indices.size()) / 1073741824.0;
    if (gb >= 1.0)
      cao.print(tick.date(), "projection in-core: the genotypes take", gb,
                "GB. -m/--memory streams them in blocks instead");
  }
  data->prepare();  // read AF and resize F to matched size; in-core, read G
  const bool standardize = data->resolve_ref_scaling(usv, params.fileS, params.project == 3);
  if (params.project != 3 && standardize && !params.out_of_core) data->standardize_E_ref(usv);
  // target number of PCs for getting individual allele frequency
  const int K = ref_pcs(params, S.size(), params.fileS);  // -k: the leading PCs, else all
  cao.print(tick.date(), "projection: using", K, "of the", S.size(), "PCs in the reference");
  // V has K columns: with -k below the reference's PC count, V * S.asDiagonal()
  // did not conform and corrupted the heap
  S = S.head(K).eval();
  Mat2D V = read_eigvecs(params.fileV, nsnps, K);
  if (!match.identical) {
    Mat2D V_overlap(match.mbim_indices.size(), K);
    for (int i = 0; i < (int)match.mbim_indices.size(); ++i) {
      V_overlap.row(i) = V.row(match.mbim_indices[i]);
      if (match.flip[i]) V_overlap.row(i) = -V_overlap.row(i);
    }
    V = V_overlap;
  }
  TextFormat fmt;
  fmt.na = "NA";  // a sample without a score
  if (params.out_of_core) {
    // -m: --project 1 and 2 on PLINK or PGEN (Cmd.cpp)
    const Mat2D U = run_projection_blockwise(data, params, V, S, usv, standardize, boot_ids);
    write_matrix(params.fileout + ".eigvecs", U, fmt);
    write_eigvecs2(params.fileout + ".eigvecs2", U, ids, fmt.na);
    return;
  }
  Mat2D U(data->nsamples, K);
  double p_miss = data->C.size() ? (double)data->C.count() / (double)data->C.size() : data->p_miss;

  if (params.project == 1) {
    if (p_miss > 0) cao.warn("there are missing genotypes. recommend using --project 2 or 3.");
    // get 1 / Singular = sqrt(Eigen * M)
    V = V * (S.array().inverse().matrix().asDiagonal());
    // G V = U D
    U = data->G * V;
  } else if (params.project == 2) {
    V = V * S.asDiagonal();
    if (p_miss == 0.0) cao.warn("there is no missing genotypes");
    solve_projection_scores(V, data->C, data->G, U);
    if (params.project_bootstrap > 0)
      write_projection_bootstrap_stats(V, data->C, data->G, U, params.project_bootstrap, params, boot_ids);
  } else if (params.project == 3) {
    // project == 3: iterative GL-aware projection (EM)
    // E-step: update expected G using individual allele frequencies PI
    // M-step: solve U with new G
    if (params.file_t != FileType::BEAGLE) cao.error("--project 3 requires BEAGLE genotype likelihood input");
    const bool filter = !data->keepSNPs.empty();

    // FileBeagle::read_all() leaves G = E[g]/2 - f; put it on the scale of the
    // matrix the reference decomposed before the first solve, as modes 1 and 2 do
    const Mat1D a = gl_site_scale(data->F, standardize, usv);
    data->G.array().rowwise() *= a.transpose().array();

    V = V * S.asDiagonal();  // VS
    solve_projection_scores(V, data->C, data->G, U);

    cao.print(tick.date(), "run EM to update expected G and solve U iteratively");
    for (uint iter = 0; iter < params.maxiter; ++iter) {
      Mat2D Uprev = U;

      // E-step: update G using individual allele frequencies, from U * V' formed
      // a panel of sites at a time instead of a scalar loop per genotype
      for_each_product_column(U, V.transpose(), 0, data->nsnps, [&](Eigen::Index jj, const auto& recon) {
        const uint j = (uint)jj;
        uint s = filter ? data->keepSNPs[j] : (uint)j;
        for (uint i = 0; i < data->nsamples; ++i) {
          double pt = recon(i);
          // pt is on the reference's scale, the same as G; undo a_j to get pi
          pt = fmin(fmax(pt / a(j) + data->F(j), 1e-4), 1.0 - 1e-4);
          const double p0 = data->P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
          const double p1 = data->P(2 * i + 1, s) * 2.0 * pt * (1.0 - pt);
          const double p2 = (1.0 - data->P(2 * i + 0, s) - data->P(2 * i + 1, s)) * pt * pt;
          const double psum = p0 + p1 + p2;
          if (!std::isfinite(psum) || psum <= 0.0) {
            data->C[(uint64)j * data->nsamples + i] = 1;
            data->G(i, j) = 0.0;
            continue;
          }
          data->C[(uint64)j * data->nsamples + i] = 0;
          data->G(i, j) = a(j) * ((p1 + 2.0 * p2) / (2.0 * psum) - data->F(j));
        }
      });

      // M-step: solve for U using expected G
      solve_projection_scores(V, data->C, data->G, U);
      // // flip signs
      // for (int k = 0; k < K; ++k) {
      //   if (U.col(k).dot(Uprev.col(k)) < 0.0) U.col(k) *= -1.0;
      // }
      double denom = Uprev.norm();
      if (denom < 1e-12) denom = 1.0;
      double diff = (U - Uprev).norm() / denom;
      cao.print(tick.date(), "GL projection iter", iter + 1, ", diff =", diff);
      if (diff < params.tolem) break;
    }
  } else {
    cao.error("unsupported --project mode: " + std::to_string(params.project));
  }

  write_matrix(params.fileout + ".eigvecs", U, fmt);
  // the scores again, with the FID and IID of each sample, as a PCA writes them
  write_eigvecs2(params.fileout + ".eigvecs2", U, ids, fmt.na);
}
