/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/EvalAdmix.cpp
 * @author      Anders Albrechtsen
 * Copyright (C) 2026. Use of this code is governed by the LICENSE file.
 *
 * Correlation of residuals for evaluating the fit of a PCA/admixture model
 * (evalAdmix, Garcia-Erill & Albrechtsen 2020; projection estimator of
 * van Waaij et al. 2023, Genetics 225:iyad157).
 *
 * The statistic is  corres = bhat - chat  where, with P the projection onto
 * [PC_1..PC_k, 1] and R = G(I-P) the residuals,
 *
 *   bhat = cov2cor( Rtilde' Rtilde ),  Rtilde = column-centred R
 *   chat = cov2cor( (I-P) Dhat (I-P) ),  Dhat = diag(mean heterozygosity)
 *
 * Rtilde'Rtilde can be written without ever forming R:
 *
 *   Rtilde'Rtilde = (I-P) [ G'G - M gbar gbar' ] (I-P)
 *
 * so one streaming pass accumulating the N x N Gram matrix G'G, the per-sample
 * mean genotype gbar and the per-sample mean heterozygosity is sufficient.
 * Memory is one N x N matrix, independent of the number of sites, and since P
 * has rank k+1 the rest costs O(N^2 k): (I-P) is never formed.
 *
 * Missing genotypes are imputed to the site mean, which is what keeps the
 * projection -- taken across individuals within a site -- well defined. The
 * identity above then holds exactly for the imputed matrix. Each missing call
 * is then replaced by its fit from the PCs, and since it carries no
 * relatedness, each pair is rescaled by the sites both samples are genotyped
 * at, counted in the same pass (a second N x N matrix, only when genotypes are
 * missing).
 *
 * Note PCAone codes genotypes on the 0..1 scale (BED2GENO = {1, NA, 0.5, 0}),
 * i.e. x = g/2. Both bhat and chat are correlation matrices, so the constant
 * factor between the x and g scales cancels and no rescaling is needed.
 *
 * The PC scores are read with Utils::read_usv(). -k below the reference's PC
 * count is a regression test for its old row-major/column-major bug, which
 * transposed any file with more than one column.
 *
 * --evaladmix-kin computes the same statistic for biobank-scale samples, where
 * the N x N matrix does not fit in memory, and writes only the pairs whose
 * kinship reaches a cutoff, plus a maximal unrelated set; see
 * run_evaladmix_pairs() below.
 ******************************************************************************/
#include "EvalAdmix.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <new>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Utils.hpp"

// "%.6f" of x, rounded exactly as printf rounds it (to nearest, ties to even),
// about ten times faster than an ostream. nan is written as "nan".
static inline char* put_fixed6(double x, char* p) {
  if (std::isnan(x)) {
    std::memcpy(p, "nan", 3);
    return p + 3;
  }
  const double ax = std::fabs(x);
  if (!(ax < 1e9)) return p + std::snprintf(p, 16, "%.6g", x);  // not reached: the entries are in [-1, 1]
  if (std::signbit(x)) *p++ = '-';
  const double y = ax * 1e6;                // 1e6 is exact, so ax * 1e6 = y + e exactly
  const double e = std::fma(ax, 1e6, -y);
  double n = std::floor(y);
  const double s = ((y - n) - 0.5) + e;     // y - n - 0.5 is exact; the sign of s is too
  if (s > 0 || (s == 0 && std::fmod(n, 2.0) != 0)) n += 1;
  uint64 q = (uint64)n, ip = q / 1000000, fp = q % 1000000;
  char tmp[20];
  int t = 0;
  do tmp[t++] = char('0' + ip % 10); while (ip /= 10);
  while (t) *p++ = tmp[--t];
  *p++ = '.';
  for (int k = 5; k >= 0; --k, fp /= 10) p[k] = char('0' + fp % 10);
  return p + 6;
}

// M is symmetric, so row i is written from column i, which is contiguous. Rows
// are formatted in parallel, a block of about 64 MB at a time, and written in order.
static void write_matrix(const std::string& fn, const Mat2D& M, const std::vector<std::string>& ids,
                         double scale = 1.0) {
  std::ofstream ofs(fn, std::ios::binary);
  if (!ofs.is_open()) cao.error("can not open file for writing: " + fn);
  if (!ids.empty()) {
    for (size_t i = 0; i < ids.size(); ++i) ofs << (i ? "\t" : "") << ids[i];
    ofs << "\n";
  }
  const Eigen::Index n = M.rows();
  const Eigen::Index width = 16 * n + 1;  // "-0.123456\t" is 10 bytes
  const Eigen::Index R = std::max<Eigen::Index>(1, std::min<Eigen::Index>(n, (Eigen::Index(1) << 26) / width));
  std::vector<std::vector<char>> buf(R, std::vector<char>(width));
  std::vector<Eigen::Index> len(R);
  for (Eigen::Index r0 = 0; r0 < n; r0 += R) {
    const Eigen::Index h = std::min(R, n - r0);
#pragma omp parallel for schedule(dynamic, 1)
    for (Eigen::Index t = 0; t < h; ++t) {
      char* p = buf[t].data();
      const double* c = M.col(r0 + t).data();
      for (Eigen::Index j = 0; j < n; ++j) {
        if (j) *p++ = '\t';
        p = put_fixed6(c[j] * scale, p);
      }
      *p++ = '\n';
      len[t] = p - buf[t].data();
    }
    for (Eigen::Index t = 0; t < h; ++t) ofs.write(buf[t].data(), len[t]);
  }
  if (!ofs) cao.error("error writing " + fn);
}

// bytes of one N x N double matrix, in GiB
static double nn_gib(Eigen::Index N) { return (double)N * (double)N * 8.0 / 1073741824.0; }

// sample IDs for the header of the output matrices: the IID column of the
// .fam, or of the .psam, whose #-header names its columns (without one it is
// laid out like a .fam). Empty when the count does not match.
static std::vector<std::string> read_sample_ids(const Param& params, Eigen::Index N) {
  const bool pgen = params.file_t == FileType::PGEN;
  std::ifstream ifs(params.filein + (pgen ? ".psam" : ".fam"));
  std::vector<std::string> ids;
  std::string line, tok;
  int col = 1;  // FID IID ...
  while (std::getline(ifs, line)) {
    if (line.empty() || (line[0] == '#' && line.size() > 1 && line[1] == '#')) continue;
    std::istringstream iss(line[0] == '#' ? line.substr(1) : line);
    if (line[0] == '#') {  // #FID IID ... or #IID ...
      for (int c = 0; iss >> tok; ++c)
        if (tok == "IID") col = c;
      continue;
    }
    for (int c = 0; c <= col && iss >> tok; ++c)
      if (c == col) ids.push_back(tok);
  }
  if ((Eigen::Index)ids.size() != N) ids.clear();
  return ids;
}

// Missing calls reach run_evaladmix() imputed to the site mean, i.e. as an
// exact 0.0 in the centred genotypes G. An observed call is exactly 0.0 only
// when it equals the site frequency f, which needs f in {0, 0.5, 1}; at every
// other site a 0.0 is a missing call. Over those sites, this
//
//   * counts ninf, the number of them, and nfull, how many have no missing call;
//   * adds o o' to Nobs for the rest, o the observed-call indicator (lower
//     triangle, allocated at the first missing call). nfull + Nobs(i, j) is
//     then the number of sites where i and j are both genotyped;
//   * replaces each missing call by its fit from the PCs, (Q Q' g)_i, clamped to
//     the genotype range. The site mean leaves the call a residual f - pi_i, the
//     ancestry deviation, which pairs share when they miss the same sites, as in
//     a genotyping batch. The fit leaves next to none.
//
// Complete data never allocates Nobs, changes nothing and costs one scan of G.
//
// impute_missing() does the counting and the imputation. Each chunk of sites
// with a missing call is shown to on_chunk(O, w) first, O holding the
// observed-call indicators of those sites in its first w columns: handle_missing()
// adds O O' to Nobs there, and --evaladmix-kin lists the missing calls instead.
template <class OnChunk>
static void impute_missing(Mat2D& G, const double* f, const Mat2D& Q, double& nfull, uint64& ninf,
                           OnChunk&& on_chunk) {
  const Eigen::Index N = G.rows(), m = G.cols();
  std::vector<char> kind(m);  // 0: missingness not visible, 1: complete, 2: has a missing call
#pragma omp parallel for schedule(static)
  for (Eigen::Index i = 0; i < m; ++i) {
    if (f[i] == 0.0 || f[i] == 0.5 || f[i] == 1.0)
      kind[i] = 0;
    else
      kind[i] = (G.col(i).array() == 0.0).any() ? 2 : 1;
  }
  std::vector<Eigen::Index> cols;
  for (Eigen::Index i = 0; i < m; ++i) {
    if (kind[i] == 0) continue;
    ++ninf;
    if (kind[i] == 1)
      nfull += 1.0;
    else
      cols.push_back(i);
  }
  if (cols.empty()) return;
  // blocks of at most ~64 MB
  const Eigen::Index ch = std::max<Eigen::Index>(64, (Eigen::Index(1) << 23) / std::max<Eigen::Index>(1, N));
  const Eigen::Index w0 = std::min<Eigen::Index>(ch, (Eigen::Index)cols.size());
  Mat2D O(N, w0), X(N, w0);
  for (size_t c0 = 0; c0 < cols.size(); c0 += ch) {
    const Eigen::Index w = std::min<Eigen::Index>(ch, (Eigen::Index)(cols.size() - c0));
#pragma omp parallel for schedule(static)
    for (Eigen::Index c = 0; c < w; ++c) {
      X.col(c) = G.col(cols[c0 + c]);
      O.col(c) = (X.col(c).array() != 0.0).cast<double>();
    }
    on_chunk(O, w);
    const Mat2D fit = Q * (Q.transpose() * X.leftCols(w));  // centred, since P1 = 1
#pragma omp parallel for schedule(static)
    for (Eigen::Index c = 0; c < w; ++c) {
      const Eigen::Index s = cols[c0 + c];
      for (Eigen::Index i = 0; i < N; ++i)
        if (O(i, c) == 0.0) G(i, s) = std::min(1.0 - f[s], std::max(-f[s], fit(i, c)));
    }
  }
}

static void handle_missing(Mat2D& G, const double* f, const Mat2D& Q, Mat2D& Nobs, double& nfull, uint64& ninf) {
  const Eigen::Index N = G.rows();
  impute_missing(G, f, Q, nfull, ninf, [&](const Mat2D& O, Eigen::Index w) {
    if (Nobs.size() == 0) Nobs = Mat2D::Zero(N, N);
    syrk_lower_add(Nobs, O.leftCols(w));
  });
}

// ---- shared by the dense and the pairwise output ------------------------

// FID and IID of each sample, from the .fam or the .psam (whose #-header names
// its columns; without one it is laid out like a .fam). Without an FID column
// has_fid is false. Empty when the count does not match.
struct SampleIds {
  std::vector<std::string> fid, iid;
  bool has_fid = true;
};

static SampleIds read_fid_iid(const Param& params, Eigen::Index N) {
  const bool pgen = params.file_t == FileType::PGEN;
  std::ifstream ifs(params.filein + (pgen ? ".psam" : ".fam"));
  SampleIds ids;
  std::string line, tok;
  int cf = 0, ci = 1;  // FID IID ...
  while (std::getline(ifs, line)) {
    if (line.empty() || (line[0] == '#' && line.size() > 1 && line[1] == '#')) continue;
    std::istringstream iss(line[0] == '#' ? line.substr(1) : line);
    if (line[0] == '#') {  // #FID IID ... or #IID ...
      cf = -1;
      for (int c = 0; iss >> tok; ++c) {
        if (tok == "IID") ci = c;
        if (tok == "FID") cf = c;
      }
      ids.has_fid = cf >= 0;
      continue;
    }
    std::string f, id;
    int c = 0;
    for (; c <= std::max(ci, cf) && iss >> tok; ++c) {
      if (c == ci) id = tok;
      if (c == cf) f = tok;
    }
    if (c <= std::max(ci, cf)) break;  // a short line: give up, as for a count mismatch
    ids.fid.push_back(f);
    ids.iid.push_back(id);
  }
  if ((Eigen::Index)ids.iid.size() != N) ids.fid.clear(), ids.iid.clear();
  return ids;
}

// The PC scores have to be of these samples in this order, which their number
// alone cannot tell. A PCA of PLINK or PGEN input writes the FID and IID of
// each row to .eigvecs2 beside .eigvecs, so the IIDs are compared with the
// .fam/.psam. Without that file only the count is checked, with a warning.
static void check_sample_order(const Param& params, Eigen::Index N) {
  const std::string& fu = params.fileU;
  const std::string ext = ".eigvecs", fids = params.filein + (params.file_t == FileType::PGEN ? ".psam" : ".fam");
  const bool named = fu.size() > ext.size() && fu.compare(fu.size() - ext.size(), ext.size(), ext) == 0;
  const std::string f2 = fu + "2";
  std::ifstream ifs;
  if (named) ifs.open(f2);
  if (!ifs.is_open()) {
    cao.warn("evalAdmix: there is no", named ? f2 : "<prefix>.eigvecs2",
             "with the sample IDs of the PC scores, so only their number is checked. the scores must be of the "
             "samples of", fids, "in the same order");
    return;
  }
  const SampleIds ids = read_fid_iid(params, N);
  if (ids.iid.empty()) {
    cao.warn("evalAdmix: could not read", N, "sample IDs from", fids, "to compare with", f2);
    return;
  }
  std::string line, fid, iid;
  Eigen::Index i = 0;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;  // #FID IID PC1 ...
    std::istringstream iss(line);
    if (!(iss >> fid >> iid)) cao.error("evalAdmix: cannot read the FID and IID on line", i + 2, "of", f2);
    if (i < N && iid != ids.iid[i])
      cao.error("evalAdmix: sample", i + 1, "is", ids.iid[i], "in", fids, "but", iid, "in", f2,
                ". the PC scores must be of the same samples in the same order");
    ++i;
  }
  if (i != N) cao.error("evalAdmix:", f2, "has", i, "samples but", fids, "has", N);
  cao.print(tick.date(), "evalAdmix: the samples of the PC scores match", fids, ", in order");
}

// Q, an orthonormal basis of span[PC_1..PC_k, 1], from the reference scores
// (-P/--USV or --read-U), which must hold these samples in genotype-file order
// (check_sample_order()).
// The projection is P = Q Q'. A rank-revealing QR drops collinear columns as
// the pseudo-inverse would. P is never formed: with r = rank <= k+1, every
// product with (I-P) below is a rank-r correction of an N x N matrix, O(N^2 r)
// instead of the O(N^3) of multiplying by a dense N x N (I-P).
static Mat2D projection_basis(const Param& params, Eigen::Index N) {
  const std::string& fpcs = params.fileU;
  if (!std::ifstream(fpcs).good()) cao.error("evalAdmix: cannot read PC scores from", fpcs);
  Mat2D U = read_usv(fpcs);  // N x kmax
  if (U.cols() == 0 || !U.allFinite()) cao.error("evalAdmix: reference PC scores must be nonempty and finite");
  cao.print(tick.date(), "evalAdmix: read", U.rows(), "x", U.cols(), "PC scores from", fpcs);
  if (U.rows() != N)
    cao.error("evalAdmix:", fpcs, "has", U.rows(), "rows but the genotype file has", N, "samples");
  check_sample_order(params, N);
  const Eigen::Index k = ref_pcs(params, U.cols(), fpcs);  // -k: the leading PCs, else all
  cao.print(tick.date(), "evalAdmix: using", k, "PC(s) + intercept =", k + 1, "dimensions");
  Mat2D V(N, k + 1);
  V.leftCols(k) = U.leftCols(k);
  V.col(k).setOnes();
  Eigen::ColPivHouseholderQR<Mat2D> qr(V);
  const Eigen::Index r = qr.rank();
  if (r < k + 1)
    cao.warn("evalAdmix: the", k, "PC(s) and the intercept span only", r,
             "dimensions; the collinear direction(s) are dropped");
  return qr.householderQ() * Mat2D::Identity(N, r);
}

// "zero" up to the rounding of (x - f) + f
static void warn_nohet(const Mat1D& d) {
  const Eigen::Index nohet = (d.array() <= 1e-9 * d.mean()).count();
  if (nohet > 0)
    cao.warn(nohet,
             "sample(s) have no heterozygous genotype. the evalAdmix null model takes each sample's variance from "
             "its heterozygosity, which assumes diploid genotypes; for these samples it is zero and their "
             "statistic is not meaningful");
}

// The per-sample pieces of step 4 of run_evaladmix(), from T = S Q, the
// diagonal of S and Dhat: off the diagonal,
//
//   corres_ij = sb_i sb_j S_ij - (L R')_ij
struct Factors {
  Mat1D sb, sc;  // 1/sqrt of the diagonals of Bcov and Ccov
  Mat2D L, R;    // N x 4r
};

static Factors make_factors(const Mat2D& Q, const Mat2D& T, const Mat1D& Sdiag, const Mat1D& d) {
  const Eigen::Index N = Q.rows(), r = Q.cols();
  const Mat2D QtT = Q.transpose() * T;
  Mat2D W = T;
  W.noalias() -= 0.5 * Q * QtT;
  const Mat2D dQ = d.asDiagonal() * Q;
  const Mat2D H = Q.transpose() * dQ;
  const Mat2D QH = Q * H;
  // cov2cor: 1/sqrt of each diagonal, floored as before
  Factors F;
  F.sb.resize(N);
  F.sc.resize(N);
  for (Eigen::Index i = 0; i < N; ++i) {
    const double bii = Sdiag(i) - 2.0 * Q.row(i).dot(W.row(i));
    const double cii = d(i) - 2.0 * dQ.row(i).dot(Q.row(i)) + Q.row(i).dot(QH.row(i));
    F.sb(i) = 1.0 / std::sqrt(std::max(bii, 1e-300));
    F.sc(i) = 1.0 / std::sqrt(std::max(cii, 1e-300));
  }
  F.L.resize(N, 4 * r);
  F.R.resize(N, 4 * r);
  F.L << F.sb.asDiagonal() * Q, F.sb.asDiagonal() * W, F.sc.asDiagonal() * Q, F.sc.asDiagonal() * dQ;
  F.R << F.sb.asDiagonal() * W, F.sb.asDiagonal() * Q, F.sc.asDiagonal() * (QH - dQ), -(F.sc.asDiagonal() * Q);
  return F;
}

static void run_evaladmix_pairs(Data* data, const Param& params);

void run_evaladmix(Data* data, const Param& params) {
  if (params.evaladmix_pairs) return run_evaladmix_pairs(data, params);
  const Eigen::Index N = data->nsamples;
  const uint M = data->nsnps;

  // ---- 0. size -----------------------------------------------------------
  // One dense N x N matrix at peak, and two dense N x N text files on the way
  // out. Both grow with the square of the sample count and neither depends on
  // the number of sites, so a cohort whose PCA runs comfortably out-of-core can
  // still be far out of reach here. Say so before spending a pass over the
  // genotypes rather than dying in the allocator afterwards.
  const double ram = nn_gib(N);  // twice that if genotypes are missing, see handle_missing()
  const double perfile = (double)N * (double)N * 9.0 / 1073741824.0;  // ~9 bytes per printed value
  cao.print(tick.date(), "evalAdmix:", N, "samples needs about", ram,
            "GB of RAM (twice that with missing genotypes), and writes two files of about", perfile, "GB each");
  if (ram > 8.0)
    cao.warn("evalAdmix needs about", ram, "GB of RAM for the", N, "x", N,
             "matrix. this does not go down with --memory, which only bounds the genotype blocks. if that is more "
             "than this machine has, use --evaladmix-kin <cutoff>, which computes the matrix in stripes within -m "
             "and writes only the pairs whose kinship reaches the cutoff.");

  // ---- 1. principal component scores, 2. projection onto [PCs, intercept] --
  // see projection_basis()
  const Mat2D Q = projection_basis(params, N);

  // ---- 3. one streaming pass: Gram matrix, mean genotype, heterozygosity --
  //
  // Both branches below see *centred* genotypes (x - f) with missing calls
  // imputed to the site mean (centred value 0) by the readers, and then to
  // their fit from the PCs by handle_missing(). That is deliberate:
  //
  //   * A and b need no correction. Per-site centring subtracts c_s * 1 from
  //     column s, and every resulting term carries a factor (I-P)1 = 0 because
  //     the projection contains the intercept. So (I-P)[A - M gbar gbar'](I-P)
  //     is identical for centred and raw genotypes.
  //   * Only d is affected, so the centring is undone column by column to get
  //     the genotype back onto the 0..1 scale.
  //
  // Asking the readers for raw genotypes instead (params.center = false) is
  // NOT an option: read_all() then leaves missing calls as BED_MISSING_VALUE
  // (-9), which poisons G*G' and makes d negative. See Cmd.cpp.
  //
  // Imputation keeps the projection (I-P), which mixes individuals within a
  // site, well defined when a genotype is absent. Imputed to its fit, a missing
  // call has next to no residual, so it adds nothing to the covariance of i and
  // j, which then runs over the n_ij sites where both are genotyped, while each
  // variance runs over n_i or n_j: the correlation shrinks by
  // n_ij / sqrt(n_i n_j) -- by the missing fraction, if calls are missing at
  // random. handle_missing() accumulates n_ij in the same pass and step 4
  // divides it out, which is what pairwise-complete sites, as evalAdmix uses,
  // would give. It is pairwise, not per sample, so that it also holds when
  // missingness is shared, e.g. by genotyping batch.
  Mat2D A;
  try {
    A = Mat2D::Zero(N, N);  // G'G, lower triangle
  } catch (const std::bad_alloc&) {
    cao.error("evalAdmix: out of memory allocating the", N, "x", N, "Gram matrix of", nn_gib(N), "GB");
  }
  Mat2D Nobs;          // pairs genotyped at the same sites; see handle_missing()
  double nfull = 0.0;   // sites where every sample is genotyped
  uint64 ninf = 0;      // sites where missingness is visible
  Mat1D b = Mat1D::Zero(N);  // sum_s g_s
  Mat1D d = Mat1D::Zero(N);  // sum_s g_s .* (1 - g_s)   (0..1 scale)
  tick.clock();
  if (!params.out_of_core) {
    // G is nsamples x nsnps, centred, not standardized.
    handle_missing(data->G, data->F.data(), Q, Nobs, nfull, ninf);
    const Mat2D& G = data->G;
    syrk_lower_add(A, G);  // lower triangle, half the flops of G * G'
    b = G.rowwise().sum();
    for (Eigen::Index i = 0; i < G.cols(); ++i) {
      const double f = data->F(i);
      d.array() += (G.col(i).array() + f) * (1.0 - (G.col(i).array() + f));
    }
    data->G.resize(0, 0);  // N x M doubles, not needed past this point
  } else {
    // Out-of-core. read_block_initial() estimates F on the fly, so F and the
    // lookup table must be allocated here before the first block.
    data->F = Mat1D::Zero(data->nsnps);
    data->centered_geno_lookup = Arr2D::Zero(4, data->nsnps);
    data->check_file_offset_first_var();
    for (uint bi = 0; bi < data->nblocks; ++bi) {
      data->read_block_initial(data->start[bi], data->stop[bi], false);
      handle_missing(data->G, data->F.data() + data->start[bi], Q, Nobs, nfull, ninf);
      const Mat2D& G = data->G;
      syrk_lower_add(A, G);
      b += G.rowwise().sum();
      for (Eigen::Index i = 0; i < G.cols(); ++i) {
        const double f = data->F(data->start[bi] + i);
        d.array() += (G.col(i).array() + f) * (1.0 - (G.col(i).array() + f));
      }
    }
    data->G.resize(0, 0);
  }
  b /= (double)M;
  d /= (double)M;
  cao.print(tick.date(), "evalAdmix: accumulated summary statistics over", M, "sites in", tick.reltime(),
            "seconds");
  warn_nohet(d);

  // ---- 4. corres = bhat - chat -------------------------------------------
  //
  //   S    = A - M gbar gbar'                   (Rtilde'Rtilde before (I-P))
  //   Bcov = (I-P) S (I-P)    = S - L R'        L = [Q, W], R = [W, Q]
  //   Ccov = (I-P) Dhat (I-P) = Dhat + Lc Rc'   Lc = [Q, dQ], Rc = [QH - dQ, -Q]
  //
  // with T = S Q, W = T - Q (Q'T)/2, dQ = Dhat Q and H = Q' Dhat Q, all N x r.
  // bhat and chat are Bcov and Ccov scaled to unit diagonal, whose diagonals
  // follow from the same factors in O(N r). So, off the diagonal,
  //
  //   corres = diag(sb) S diag(sb) - [sb.L, sc.Lc] [sb.R, sc.Rc]'
  //
  // with sb, sc the inverse square roots of the two diagonals: one scaling of A
  // and one N x N x 4r product, done in place. Only A is ever N x N.
  tick.clock();
  A.selfadjointView<Eigen::Lower>().rankUpdate(b, -(double)M);  // S, lower triangle
  mirror_lower(A);
  Factors F;
  {
    Mat2D T;
    T.noalias() = A * Q;
    F = make_factors(Q, T, A.diagonal(), d);
  }
#pragma omp parallel for schedule(static)
  for (Eigen::Index j = 0; j < N; ++j) A.col(j).array() *= F.sb.array() * F.sb(j);
  A.noalias() -= F.L * F.R.transpose();
  // undo the attenuation by missing calls (see step 3): multiply by
  // sqrt(n_i n_j) / n_ij. a pair never genotyped at the same site has no
  // estimate and is written as nan
  uint64 nopair = 0;
  if (Nobs.size() > 0) {
    mirror_lower(Nobs);
    const Mat1D ni = Nobs.diagonal().array() + nfull;
    cao.print(tick.date(), "evalAdmix: missing genotypes were imputed from the PCs; samples are genotyped at",
              ni.minCoeff() / ninf, "to", ni.maxCoeff() / ninf,
              "of the sites, and each pair is rescaled by its sites in common");
#pragma omp parallel for schedule(static) reduction(+ : nopair)
    for (Eigen::Index j = 0; j < N; ++j) {
      for (Eigen::Index i = 0; i < N; ++i) {
        const double nij = Nobs(i, j) + nfull;
        if (nij > 0) {
          A(i, j) *= std::sqrt(ni(i) * ni(j)) / nij;
        } else {
          A(i, j) = std::numeric_limits<double>::quiet_NaN();
          nopair += i != j;
        }
      }
    }
    Nobs.resize(0, 0);
  }
  // a correlation cannot leave [-1, 1]; neither may the difference of two.
  // scalar, so that nan passes through as nan
#pragma omp parallel for schedule(static)
  for (Eigen::Index j = 0; j < N; ++j)
    for (Eigen::Index i = 0; i < N; ++i)
      if (!std::isnan(A(i, j))) A(i, j) = std::min(1.0, std::max(-1.0, A(i, j)));
  A.diagonal().setZero();
  mirror_lower(A);  // exactly symmetric, as write_matrix() relies on
  if (nopair > 0)
    cao.warn(nopair / 2, "pair(s) of samples are never genotyped at the same site; their entries are nan");
  cao.print(tick.date(), "evalAdmix: correlation of residuals computed in", tick.reltime(), "seconds");

  // ---- 5. output ---------------------------------------------------------
  const std::vector<std::string> ids = read_sample_ids(params, N);
  if (ids.empty())
    cao.warn("evalAdmix: could not read", N, "sample IDs from",
             params.filein + (params.file_t == FileType::PGEN ? ".psam" : ".fam") + "; the output has no header line");
  write_matrix(params.fileout + ".corres", A, ids);
  // the correlation of residuals estimates 2*phi, so halve it for kinship. done
  // at write time rather than in another N x N copy; the diagonal is already 0.
  write_matrix(params.fileout + ".kinship", A, ids, 0.5);
  cao.print(tick.date(), "evalAdmix: correlation of residuals saved to", params.fileout + ".corres");
  cao.print(tick.date(), "evalAdmix: kinship (corres/2) saved to", params.fileout + ".kinship");
}

// ===========================================================================
// --evaladmix-kin: the pairs above a kinship cutoff, for biobank-scale data
// ===========================================================================
//
// The dense statistic above holds the N x N Gram matrix: 480 GB at N = 245,000.
// But in corres_ij = sb_i sb_j S_ij - (L R')_ij only S_ij is a pair quantity.
// Everything else is per sample: gbar, Dhat, the diagonal of S and T = S Q,
// which is G (G' Q) - M gbar (gbar' Q) and never needs S itself. So
//
//   pass 1     accumulates those, an O(N M r) pass, and each sample's missing
//              calls, and make_factors() turns them into sb, sc, L and R;
//   pass 2..   compute S one stripe of samples at a time: S(j, i) for every i
//              in the stripe and every j > i, the product of the stripe's rows
//              of G with the rows from the stripe down. A stripe is sized to
//              fit in -m; each is completed to corres in place and only the
//              pairs whose kinship (corres / 2) reaches the cutoff are written.
//
// Entry by entry the arithmetic is that of run_evaladmix(), so the two agree to
// rounding, and the total work is its one O(N^2 M) Gram product, cut into
// triangular stripes. Each stripe costs one more read of the genotypes: from
// RAM in-core, from the file with -m.
//
// With missing genotypes the rescaling needs the sites each pair has in common,
//
//   n_ij = ninf - m_i - m_j + mm_ij,
//
// with m_i the missing calls of sample i and mm_ij the sites where i and j are
// both missing. At biobank call rates nearly every site has a missing call but
// only a few of the samples miss it, so mm is counted pair by pair from the
// lists of missing calls (add_both_missing()) rather than as the dense O O'
// product of handle_missing(), which would double the cost of the Gram product.

using MatF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;

// The missing calls of a block of sites, one ascending list of samples per site
// that has any. Only the sites where missingness is visible; see handle_missing().
struct MissingCalls {
  std::vector<uint64> ptr{0};  // site s lists idx[ptr[s] .. ptr[s+1])
  std::vector<uint32_t> idx;
  size_t nsites() const { return ptr.size() - 1; }
  void clear() {
    ptr.assign(1, 0);
    idx.clear();
  }
  // the samples missing in each of the first w columns of O, observed-call indicators
  void add(const Mat2D& O, Eigen::Index w) {
    const size_t s0 = nsites();
    ptr.resize(s0 + w + 1);
#pragma omp parallel for schedule(static)
    for (Eigen::Index c = 0; c < w; ++c) ptr[s0 + c + 1] = (O.col(c).array() == 0.0).count();
    for (Eigen::Index c = 0; c < w; ++c) ptr[s0 + c + 1] += ptr[s0 + c];
    idx.resize(ptr.back());
#pragma omp parallel for schedule(static)
    for (Eigen::Index c = 0; c < w; ++c) {
      uint64 p = ptr[s0 + c];
      for (Eigen::Index i = 0; i < O.rows(); ++i)
        if (O(i, c) == 0.0) idx[p++] = (uint32_t)i;
    }
  }
};

// MM(j - i0, i - i0) += the number of sites of ms where sample i of the stripe
// [i0, i0 + nb) and sample j > i are both missing. Entries with j <= i are left
// unspecified. A site missing in few samples of the stripe is scattered pair by
// pair, each sample of the stripe adding to its own column; a site missing in
// so many that this would cost more than a dense product goes into a float
// product of missing indicators. Counts are exact in float up to 2^24 sites.
static void add_both_missing(MatF& MM, const MissingCalls& ms, Eigen::Index i0, Eigen::Index nb) {
  const Eigen::Index Np = MM.rows();
  const uint32_t lo = (uint32_t)i0, hi = (uint32_t)(i0 + nb);
  const size_t S = ms.nsites();
  const uint32_t* all = ms.idx.data();
  std::vector<uint64> first(S), mid(S);  // positions of the first sample >= i0, and >= i0 + nb
  std::vector<size_t> heavy;
  std::vector<uint64> cnt(nb + 1, 0);
  for (size_t s = 0; s < S; ++s) {
    const uint32_t* e = all + ms.ptr[s + 1];
    const uint32_t* x = std::lower_bound(all + ms.ptr[s], e, lo);
    const uint32_t* y = std::lower_bound(x, e, hi);
    first[s] = x - all;
    mid[s] = y - all;
    const double nI = (double)(y - x), nJ = (double)(e - x);
    if (nI == 0) continue;
    if (nI * nJ * 32.0 > (double)nb * (double)Np) {
      heavy.push_back(s);
      mid[s] = first[s];  // nothing to scatter
    } else {
      for (const uint32_t* t = x; t < y; ++t) ++cnt[*t - lo + 1];
    }
  }
  // for each sample of the stripe, where it sits in the lists of its sites
  for (Eigen::Index c = 0; c < nb; ++c) cnt[c + 1] += cnt[c];
  std::vector<uint64> pos(cnt.back()), end(cnt.back());
  {
    std::vector<uint64> fill(cnt.begin(), cnt.end() - 1);
    for (size_t s = 0; s < S; ++s)
      for (uint64 p = first[s]; p < mid[s]; ++p) {
        const uint64 q = fill[all[p] - lo]++;
        pos[q] = p;
        end[q] = ms.ptr[s + 1];
      }
  }
#pragma omp parallel for schedule(dynamic, 16)
  for (Eigen::Index c = 0; c < nb; ++c) {
    float* col = MM.col(c).data();
    for (uint64 q = cnt[c]; q < cnt[c + 1]; ++q)
      for (uint64 p = pos[q] + 1; p < end[q]; ++p) col[all[p] - lo] += 1.0f;  // the samples after i: j > i
  }
  for (size_t h0 = 0; h0 < heavy.size(); h0 += 256) {
    const Eigen::Index g = (Eigen::Index)std::min<size_t>(256, heavy.size() - h0);
    MatF XJ = MatF::Zero(Np, g), XI = MatF::Zero(nb, g);
#pragma omp parallel for schedule(static)
    for (Eigen::Index t = 0; t < g; ++t) {
      const size_t s = heavy[h0 + t];
      for (uint64 p = first[s]; p < ms.ptr[s + 1]; ++p) {
        const uint32_t j = all[p];
        XJ(j - lo, t) = 1.0f;
        if (j < hi) XI(j - lo, t) = 1.0f;
      }
    }
    MM.noalias() += XJ * XI.transpose();
  }
}

// unsigned integer to decimal
static inline char* put_uint(uint64 v, char* p) {
  char tmp[20];
  int t = 0;
  do tmp[t++] = char('0' + v % 10);
  while (v /= 10);
  while (t) *p++ = tmp[--t];
  return p;
}

// A maximal set of samples without a pair among the edges, by the greedy rule
// of Hail's maximal_independent_set and plink2's --king-cutoff: drop the sample
// with the most relatives left, until none has any; then take back each dropped
// sample none of whose relatives was kept. Ties drop the sample with more
// missing calls, then the later one in the file.
static std::vector<char> unrelated_set(Eigen::Index N, const std::vector<std::pair<uint32_t, uint32_t>>& edges,
                                       const std::vector<uint32_t>& miss) {
  std::vector<uint64> ptr(N + 1, 0);
  for (const auto& e : edges) ++ptr[e.first + 1], ++ptr[e.second + 1];
  for (Eigen::Index i = 0; i < N; ++i) ptr[i + 1] += ptr[i];
  std::vector<uint32_t> adj(ptr[N]), deg(N);
  {
    std::vector<uint64> fill(ptr.begin(), ptr.end() - 1);
    for (const auto& e : edges) adj[fill[e.first]++] = e.second, adj[fill[e.second]++] = e.first;
  }
  using Key = std::tuple<uint32_t, uint32_t, uint32_t>;  // relatives left, missing calls, sample
  std::priority_queue<Key> heap;
  for (Eigen::Index i = 0; i < N; ++i) {
    deg[i] = (uint32_t)(ptr[i + 1] - ptr[i]);
    if (deg[i] > 0) heap.emplace(deg[i], miss[i], (uint32_t)i);
  }
  std::vector<char> keep(N, 1);
  std::vector<uint32_t> dropped;
  while (!heap.empty()) {
    const auto [dg, mi, v] = heap.top();
    heap.pop();
    if (!keep[v] || dg != deg[v]) continue;  // stale
    keep[v] = 0;
    dropped.push_back(v);
    for (uint64 p = ptr[v]; p < ptr[v + 1]; ++p) {
      const uint32_t u = adj[p];
      if (keep[u] && --deg[u] > 0) heap.emplace(deg[u], miss[u], u);
    }
  }
  for (auto it = dropped.rbegin(); it != dropped.rend(); ++it) {
    bool free = true;
    for (uint64 p = ptr[*it]; p < ptr[*it + 1] && free; ++p) free = !keep[adj[p]];
    if (free) keep[*it] = 1;
  }
  return keep;
}

static void run_evaladmix_pairs(Data* data, const Param& params) {
  const Eigen::Index N = data->nsamples;
  const uint M = data->nsnps;
  const double kmin = params.evaladmix_kin, kunrel = params.evaladmix_unrel;
  const bool want_unrel = kunrel > 0;  // unrelated pairs scatter about 0
  cao.print(tick.date(), "evalAdmix: writing the pairs with kinship >=", kmin, "of", N, "samples");
  std::remove((params.fileout + ".unrelated").c_str());  // not to be mistaken for this run's
  const Mat2D Q = projection_basis(params, N);
  const Eigen::Index r = Q.cols();

  // ---- pass 1: the per-sample sums ----------------------------------------
  // the same centred, imputed genotypes as run_evaladmix() step 3
  Mat1D a = Mat1D::Zero(N);     // sum_s g_s^2, the diagonal of A = G'G
  Mat1D b = Mat1D::Zero(N);     // sum_s g_s
  Mat1D d = Mat1D::Zero(N);     // sum_s g_s .* (1 - g_s)   (0..1 scale)
  Mat2D T = Mat2D::Zero(N, r);  // A Q
  std::vector<uint32_t> miss(N, 0);  // m_i: missing calls where missingness is visible
  double nfull = 0.0;  // sites where every sample is genotyped
  uint64 ninf = 0;     // sites where missingness is visible
  MissingCalls ms;     // in-core, of all the sites and kept for the stripes; else of one block
  auto record = [&ms](const Mat2D& O, Eigen::Index w) { ms.add(O, w); };
  auto accumulate = [&](Mat2D& G, const double* f) {
    const size_t first = ms.idx.size();
    impute_missing(G, f, Q, nfull, ninf, record);
    for (size_t t = first; t < ms.idx.size(); ++t) ++miss[ms.idx[t]];
    const Eigen::Index w = G.cols();
#pragma omp parallel for schedule(static)
    for (Eigen::Index i0 = 0; i0 < N; i0 += 1024) {
      const Eigen::Index h = std::min<Eigen::Index>(1024, N - i0);
      for (Eigen::Index s = 0; s < w; ++s) {
        const auto g = G.col(s).segment(i0, h).array();
        a.segment(i0, h).array() += g.square();
        b.segment(i0, h).array() += g;
        d.segment(i0, h).array() += (g + f[s]) * (1.0 - (g + f[s]));
      }
    }
    const Mat2D Z = G.transpose() * Q;  // w x r
    T.noalias() += G * Z;
  };
  tick.clock();
  if (!params.out_of_core) {
    accumulate(data->G, data->F.data());
  } else {
    // read_block_initial() estimates F on the fly, see run_evaladmix()
    data->F = Mat1D::Zero(data->nsnps);
    data->centered_geno_lookup = Arr2D::Zero(4, data->nsnps);
    data->check_file_offset_first_var();
    for (uint bi = 0; bi < data->nblocks; ++bi) {
      data->read_block_initial(data->start[bi], data->stop[bi], false);
      ms.clear();
      accumulate(data->G, data->F.data() + data->start[bi]);
    }
    ms.clear();
  }
  b /= (double)M;
  d /= (double)M;
  cao.print(tick.date(), "evalAdmix: accumulated per-sample sums over", M, "sites in", tick.reltime(), "seconds");
  warn_nohet(d);
  // S = A - M gbar gbar': its product with Q, and its diagonal
  T.noalias() -= ((double)M * b) * (b.transpose() * Q);
  const Mat1D Sd = a.array() - (double)M * b.array().square();
  const Factors F = make_factors(Q, T, Sd, d);
  T.resize(0, 0);
  a.resize(0);

  const uint64 nmsites = ninf - (uint64)nfull;  // sites with a missing call
  const bool missing = nmsites > 0;
  Mat1D ni;  // n_i, the sites sample i is genotyped at
  if (missing) {
    if (nmsites >= (uint64(1) << 24))
      cao.error("evalAdmix:", nmsites,
                "sites have missing genotypes, more than the 16,777,216 the pair counts of --evaladmix-kin hold "
                "exactly. LD-prune or filter the sites first");
    ni.resize(N);
    for (Eigen::Index i = 0; i < N; ++i) ni(i) = (double)ninf - miss[i];
    cao.print(tick.date(), "evalAdmix: missing genotypes were imputed from the PCs; samples are genotyped at",
              ni.minCoeff() / ninf, "to", ni.maxCoeff() / ninf,
              "of the sites, and each pair is rescaled by its sites in common");
  }

  // ---- the stripes ----------------------------------------------------------
  // a stripe of nb samples from i0 holds (N - i0) x nb pairs: S in double, and
  // the pair counts in float if genotypes are missing. In-core the genotypes
  // are in RAM already and the stripes take up to 2 GB; with -m they take what
  // is left of it by the genotype block, the per-sample factors and, with
  // missing genotypes, the three chunks impute_missing() works in, a block's
  // lists of missing calls and the indicator panels of add_both_missing().
  const double gib = 1073741824.0, elem = missing ? 12.0 : 8.0;
  double budget = 2.0 * gib, fixed = 0.0;
  if (params.out_of_core) {
    const double bs = data->blocksize, chunk = std::min(bs, std::max(64.0, std::floor(8388608.0 / N)));
    double pmiss = 0.0;  // the share of calls missing
    for (Eigen::Index i = 0; i < N; ++i) pmiss += miss[i];
    pmiss /= std::max(1.0, (double)N * ninf);
    fixed = 8.0 * N * bs + 8.0 * N * (12.0 * r + 8.0) +
            (missing ? 24.0 * N * chunk + 20.0 * N * bs * pmiss + 4.0 * 256 * N : 0.0);
    budget = params.memory * gib - fixed;
  }
  std::vector<Eigen::Index> s0{0};
  bool floor64 = false;
  for (Eigen::Index i0 = 0; i0 < N;) {
    const Eigen::Index Np = N - i0, fit = (Eigen::Index)std::max(0.0, std::floor(budget / (elem * Np)));
    const Eigen::Index nb = std::min(Np, std::max(std::min<Eigen::Index>(64, Np), fit));
    floor64 |= fit < nb;
    i0 += nb;
    s0.push_back(i0);
  }
  const size_t nstripes = s0.size() - 1;
  const double maxgb = elem * (double)N * (double)(s0[1] - s0[0]) / gib;
  cao.print(tick.date(), "evalAdmix:", (uint64)N * (N - 1) / 2, "pairs in", nstripes, "stripe(s) of up to", maxgb,
            "GB, each one more pass over the genotypes", params.out_of_core ? "in the file" : "in RAM");
  if (floor64 && params.out_of_core)
    cao.warn("-m", params.memory, "GB leaves room for stripes of only 64 samples, which makes", nstripes,
             "passes over the genotypes, and needs about", (fixed + elem * 64.0 * N) / gib,
             "GB at the least. raise -m to make fewer passes");

  SampleIds ids = read_fid_iid(params, N);
  if (ids.iid.empty()) {
    cao.warn("evalAdmix: could not read", N, "sample IDs from",
             params.filein + (params.file_t == FileType::PGEN ? ".psam" : ".fam") +
                 "; the samples are numbered 1..N instead");
    ids.has_fid = false;
    for (Eigen::Index i = 0; i < N; ++i) ids.iid.push_back(std::to_string(i + 1));
  }
  std::vector<std::string> pre(N);  // "FID\tIID\t" or "IID\t"
  for (Eigen::Index i = 0; i < N; ++i) pre[i] = (ids.has_fid ? ids.fid[i] + "\t" : "") + ids.iid[i] + "\t";
  const std::string fkin = params.fileout + ".kin0";
  std::ofstream out(fkin, std::ios::binary);
  if (!out.is_open()) cao.error("can not open file for writing: " + fkin);
  out << (ids.has_fid ? "#FID1\tIID1\tFID2\tIID2" : "#IID1\tIID2") << "\tNSNP\tKINSHIP\n";

  // pairs by KING degree bin (Manichaikul et al. 2010), and a 4th degree
  const double edges[5] = {0.354, 0.177, 0.0884, 0.0442, 0.0221};
  uint64 nbin0 = 0, nbin1 = 0, nbin2 = 0, nbin3 = 0, nbin4 = 0, nbin5 = 0, npairs = 0, nopair = 0;
  // over the pairs below 0, which relatives do not reach, for the log: how
  // many there are, how far they spread, and how many reach -cutoff
  uint64 nmirror = 0, nneg = 0;
  double kneg2 = 0.0;
  std::vector<std::pair<uint32_t, uint32_t>> rel;  // the pairs >= kunrel
  const size_t relmax = size_t(1) << 25;          // 256 MB of them
  bool rel_overflow = false;
  const Eigen::Index CH = std::max<Eigen::Index>(8, 2 * omp_get_max_threads());
  std::vector<std::string> buf(CH);
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> relc(CH);
  Mat2D As;
  MatF MM;
  const double npairs_all = 0.5 * (double)N * (double)(N - 1);
  double pairs_done = 0.0, secs = 0.0;
  for (size_t st = 0; st < nstripes; ++st) {
    tick.clock();
    const Eigen::Index i0 = s0[st], nb = s0[st + 1] - i0, Np = N - i0;
    try {
      As.setZero(Np, nb);  // As(j - i0, i - i0) = S(j, i), j >= i0
      if (missing) MM.setZero(Np, nb);
    } catch (const std::bad_alloc&) {
      cao.error("evalAdmix: out of memory allocating a stripe of", elem * Np * nb / gib, "GB. lower -m");
    }
    // the stripe's own nb x nb block needs only its lower triangle: row panels
    // of 512 samples, each up to its diagonal, waste 256 / nb of it and, unlike
    // syrk_lower_add(), allocate nothing. Below it, one product.
    auto add_block = [&](const Mat2D& G, const MissingCalls& m) {
      for (Eigen::Index r0 = 0; r0 < nb; r0 += 512) {
        const Eigen::Index q = std::min<Eigen::Index>(512, nb - r0);
        As.block(r0, 0, q, r0 + q).noalias() += G.middleRows(i0 + r0, q) * G.middleRows(i0, r0 + q).transpose();
      }
      if (Np > nb)
        As.bottomRows(Np - nb).noalias() += G.middleRows(i0 + nb, Np - nb) * G.middleRows(i0, nb).transpose();
      if (missing) add_both_missing(MM, m, i0, nb);
    };
    if (!params.out_of_core) {
      add_block(data->G, ms);
    } else {
      double nf = 0.0;  // counted in pass 1 already
      uint64 nn = 0;
      data->check_file_offset_first_var();
      for (uint bi = 0; bi < data->nblocks; ++bi) {
        data->read_block_initial(data->start[bi], data->stop[bi], false);
        ms.clear();
        impute_missing(data->G, data->F.data() + data->start[bi], Q, nf, nn, record);
        add_block(data->G, ms);
      }
    }
    // A -> S -> diag(sb) S diag(sb) -> corres before the rescaling, as in
    // run_evaladmix() step 4. Only j > i is used
#pragma omp parallel for schedule(dynamic, 16)
    for (Eigen::Index c = 0; c < nb; ++c) {
      const Eigen::Index i = i0 + c;
      const double mbi = (double)M * b(i), sbi = F.sb(i);
      for (Eigen::Index rr = c + 1; rr < Np; ++rr)
        As(rr, c) = (As(rr, c) - mbi * b(i0 + rr)) * (F.sb(i0 + rr) * sbi);
    }
    As.noalias() -= F.L.middleRows(i0, Np) * F.R.middleRows(i0, nb).transpose();
    // rescale, clip and write, a chunk of the stripe's samples at a time
    for (Eigen::Index c0 = 0; c0 < nb; c0 += CH) {
      const Eigen::Index h = std::min(CH, nb - c0);
#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+ : npairs, nopair, nbin0, nbin1, nbin2, nbin3, nbin4, nbin5, nmirror, nneg, kneg2)
      for (Eigen::Index t = 0; t < h; ++t) {
        const Eigen::Index c = c0 + t, i = i0 + c;
        std::string& o = buf[t];
        o.clear();
        relc[t].clear();
        char num[48];
        for (Eigen::Index rr = c + 1; rr < Np; ++rr) {
          const Eigen::Index j = i0 + rr;
          double v = As(rr, c), nsnp = (double)M - miss[i] - miss[j];
          if (missing) {
            const double mm = MM(rr, c), nij = (double)ninf - miss[i] - miss[j] + mm;
            nsnp += mm;
            if (!(nij > 0)) {
              ++nopair;
              continue;
            }
            v *= std::sqrt(ni(j) * ni(i)) / nij;
          }
          v = std::min(1.0, std::max(-1.0, v));
          const double kin = v * 0.5;
          if (kin < 0) {
            ++nneg;
            kneg2 += kin * kin;
            nmirror += kin <= -kmin;
          }
          if (!(kin >= kmin)) continue;
          ++npairs;
          if (kin >= edges[0])
            ++nbin0;
          else if (kin >= edges[1])
            ++nbin1;
          else if (kin >= edges[2])
            ++nbin2;
          else if (kin >= edges[3])
            ++nbin3;
          else if (kin >= edges[4])
            ++nbin4;
          else
            ++nbin5;
          o += pre[i];
          o += pre[j];
          char* p = put_uint((uint64)(nsnp + 0.5), num);
          *p++ = '\t';
          p = put_fixed6(kin, p);
          *p++ = '\n';
          o.append(num, p - num);
          if (want_unrel && kin >= kunrel) relc[t].emplace_back((uint32_t)i, (uint32_t)j);
        }
      }
      for (Eigen::Index t = 0; t < h; ++t) {
        out.write(buf[t].data(), buf[t].size());
        if (rel_overflow) continue;
        rel.insert(rel.end(), relc[t].begin(), relc[t].end());
        if (rel.size() > relmax) {
          rel_overflow = true;
          std::vector<std::pair<uint32_t, uint32_t>>().swap(rel);
        }
      }
    }
    if (!out) cao.error("error writing " + fkin);
    const double t = tick.reltime();
    secs += t;
    pairs_done += (double)nb * (double)Np - 0.5 * (double)nb * (double)(nb + 1);
    cao.print(tick.date(), "evalAdmix: stripe", st + 1, "of", nstripes, "(samples", i0 + 1, "to", i0 + nb, ") in", t,
              "seconds;", npairs, "pair(s) written so far; about",
              (uint64)std::ceil(secs * (npairs_all - pairs_done) / std::max(1.0, pairs_done)), "seconds left");
  }
  As.resize(0, 0);
  MM.resize(0, 0);
  data->G.resize(0, 0);
  out.close();
  if (!out) cao.error("error writing " + fkin);
  if (nopair > 0)
    cao.warn(nopair, "pair(s) of samples are never genotyped at the same site; they have no estimate and are left out");
  cao.print(tick.date(), "evalAdmix:", npairs, "pair(s) with kinship >=", kmin, "saved to", fkin);
  // How far chance alone scatters the kinship of unrelated pairs, when the PCs
  // fit: corres is a correlation over the sites, each weighted by its genotype
  // variance, so its sd is about 1/sqrt(Meff), with Meff = (sum v)^2 / sum v^2
  // and v = f(1-f), scaled by the share of the sites a pair has in common.
  // Kinship halves it. The pairs scatter about 0 (the projection takes out the
  // mean), and those below 0 are free of relatives, so their root mean square
  // measures the scatter. Well above chance, it means the PCs leave structure,
  // whose residual correlation reads as relatedness in both tails; a cutoff
  // within a few sd of chance lets in unrelated pairs.
  {
    const double nall = npairs_all - (double)nopair, ksd = std::sqrt(kneg2 / std::max<double>(1, nneg));
    double s1 = 0.0, s2 = 0.0;
    for (uint si = 0; si < M; ++si) {
      const double v = data->F(si) * (1.0 - data->F(si));
      s1 += v, s2 += v * v;
    }
    const double meff = s2 > 0 ? s1 * s1 / s2 : 0.0, share = missing ? ni.mean() / (double)ninf : 1.0;
    const double sdchance = meff > 0 ? 0.5 / (share * std::sqrt(meff)) : 0.0;
    const double nchance = kmin > 0 && sdchance > 0 ? nall * 0.5 * std::erfc(kmin / sdchance / std::sqrt(2.0)) : 0.0;
    cao.print(tick.date(), "evalAdmix: unrelated pairs scatter with an sd of", ksd, "(from the", nneg,
              "pairs below 0); chance alone gives about", sdchance, "with", meff, "effective sites.", nmirror,
              "pair(s) are at or below", -kmin);
    if (sdchance > 0 && ksd > 1.3 * sdchance)
      cao.warn("evalAdmix: unrelated pairs scatter", ksd / sdchance, "times as much as chance alone: the PCs leave "
               "population structure (too few PCs, or fine-scale or founder groups), whose residual correlation "
               "reads as relatedness. see 'Getting the number of PCs wrong' in docs/evaladmix.md");
    else if (nchance > 0.1 * std::max(1.0, (double)npairs))
      cao.warn("evalAdmix: about", nchance, "unrelated pairs reach kinship", kmin, "by chance alone, against the",
               npairs, "written: the cutoff is within the noise of", M, "sites. raise it, or use more sites");
  }
  cao.print(tick.date(), "evalAdmix: of those, by KING degree: duplicate/MZ (>= 0.354):", nbin0,
            ", 1st (>= 0.177):", nbin1, ", 2nd (>= 0.0884):", nbin2, ", 3rd (>= 0.0442):", nbin3,
            ", 4th (>= 0.0221):", nbin4, ", lower:", nbin5);

  // ---- the unrelated set ----------------------------------------------------
  if (!want_unrel) {
    cao.print(tick.date(), "evalAdmix: no .unrelated set, which needs a positive kinship cutoff "
              "(--evaladmix-unrelated)");
    return;
  }
  if (rel_overflow) {
    cao.warn("evalAdmix: more than", relmax, "pairs have kinship >=", kunrel,
             "; no .unrelated set is written. raise --evaladmix-unrelated");
    return;
  }
  std::vector<char> keep = unrelated_set(N, rel, miss);
  std::vector<char> hasrel(N, 0);
  for (const auto& e : rel) hasrel[e.first] = hasrel[e.second] = 1;
  const Eigen::Index nrel = std::count(hasrel.begin(), hasrel.end(), 1);
  // a sample genotyped at none of the sites has no estimate, so no relatives
  Eigen::Index nempty = 0;
  if (missing)
    for (Eigen::Index i = 0; i < N; ++i)
      if (ni(i) == 0 && keep[i]) keep[i] = 0, ++nempty;
  if (nempty > 0)
    cao.warn(nempty, "sample(s) are genotyped at none of the sites and are left out of the .unrelated set");
  const Eigen::Index nkeep = std::count(keep.begin(), keep.end(), 1);
  const std::string fun = params.fileout + ".unrelated";
  std::ofstream ou(fun, std::ios::binary);
  if (!ou.is_open()) cao.error("can not open file for writing: " + fun);
  for (Eigen::Index i = 0; i < N; ++i)
    if (keep[i]) ou << (ids.has_fid ? ids.fid[i] + "\t" : "") << ids.iid[i] << "\n";
  ou.close();
  if (!ou) cao.error("error writing " + fun);
  cao.print(tick.date(), "evalAdmix:", nrel, "sample(s) are in the", rel.size(), "pair(s) with kinship >=", kunrel,
            "; dropping", N - nkeep - nempty, "of them leaves", nkeep, "unrelated sample(s), saved to", fun);
}
