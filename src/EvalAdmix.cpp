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
 * The PC scores are read with Utils::read_usv(). --evaladmix-k is a regression
 * test for its old row-major/column-major bug, which transposed any file with
 * more than one column.
 ******************************************************************************/
#include "EvalAdmix.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <new>
#include <sstream>
#include <string>
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
static void handle_missing(Mat2D& G, const double* f, const Mat2D& Q, Mat2D& Nobs, double& nfull, uint64& ninf) {
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
  if (Nobs.size() == 0) Nobs = Mat2D::Zero(N, N);
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
    syrk_lower_add(Nobs, O.leftCols(w));
    const Mat2D fit = Q * (Q.transpose() * X.leftCols(w));  // centred, since P1 = 1
#pragma omp parallel for schedule(static)
    for (Eigen::Index c = 0; c < w; ++c) {
      const Eigen::Index s = cols[c0 + c];
      for (Eigen::Index i = 0; i < N; ++i)
        if (O(i, c) == 0.0) G(i, s) = std::min(1.0 - f[s], std::max(-f[s], fit(i, c)));
    }
  }
}

void run_evaladmix(Data* data, const Param& params) {
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
             "matrix. this does not go down with --memory, which only bounds the genotype blocks. reduce the "
             "sample set if that is more than this machine has.");

  // ---- 1. principal component scores -------------------------------------
  // Always the .eigvecs this run just wrote, never params.fileU. The statistic
  // is the residual of THIS genotype matrix after projecting out THESE PCs, so
  // taking the scores from a -P/--USV prefix would silently pair one run's
  // residuals with another run's ancestry, and only a row-count mismatch would
  // ever catch it.
  const std::string fpcs = params.fileout + ".eigvecs";
  Mat2D U = read_usv(fpcs);  // N x kmax
  cao.print(tick.date(), "evalAdmix: read", U.rows(), "x", U.cols(), "PC scores from", fpcs);
  if (U.rows() != N)
    cao.error("evalAdmix:", fpcs, "has", U.rows(), "rows but the genotype file has", N, "samples");
  Eigen::Index k = params.evaladmix_k > 0 ? params.evaladmix_k : U.cols();
  if (k > U.cols())
    cao.error("--evaladmix-k is", k, "but only", U.cols(), "PCs were computed; raise -k or lower --evaladmix-k");
  cao.print(tick.date(), "evalAdmix: using", k, "PC(s) + intercept =", k + 1, "dimensions");

  // ---- 2. projection onto [PCs, intercept] -------------------------------
  // P = Q Q' with Q an orthonormal basis of span[PC_1..PC_k, 1], found by a
  // rank-revealing QR so that collinear columns are dropped as the
  // pseudo-inverse would drop them. P is never formed: with r = rank <= k+1,
  // every product with (I-P) below is a rank-r correction of an N x N matrix,
  // O(N^2 r) instead of the O(N^3) of multiplying by a dense N x N (I-P).
  Mat2D Q;
  {
    Mat2D V(N, k + 1);
    V.leftCols(k) = U.leftCols(k);
    V.col(k).setOnes();
    Eigen::ColPivHouseholderQR<Mat2D> qr(V);
    const Eigen::Index r = qr.rank();
    if (r < k + 1)
      cao.warn("evalAdmix: the", k, "PC(s) and the intercept span only", r,
               "dimensions; the collinear direction(s) are dropped");
    Q = qr.householderQ() * Mat2D::Identity(N, r);
  }

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
  // (-9), which poisons G*G' and makes d negative. See Main.cpp.
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
  // "zero" up to the rounding of (x - f) + f
  const Eigen::Index nohet = (d.array() <= 1e-9 * d.mean()).count();
  if (nohet > 0)
    cao.warn(nohet,
             "sample(s) have no heterozygous genotype. the evalAdmix null model takes each sample's variance from "
             "its heterozygosity, which assumes diploid genotypes; for these samples it is zero and their "
             "statistic is not meaningful");

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
  const Eigen::Index r = Q.cols();
  Mat2D T;
  T.noalias() = A * Q;
  const Mat2D QtT = Q.transpose() * T;
  Mat2D W = T;
  W.noalias() -= 0.5 * Q * QtT;
  T.resize(0, 0);
  const Mat2D dQ = d.asDiagonal() * Q;
  const Mat2D H = Q.transpose() * dQ;
  const Mat2D QH = Q * H;
  // cov2cor: 1/sqrt of each diagonal, floored as before
  Mat1D sb(N), sc(N);
  for (Eigen::Index i = 0; i < N; ++i) {
    const double bii = A(i, i) - 2.0 * Q.row(i).dot(W.row(i));
    const double cii = d(i) - 2.0 * dQ.row(i).dot(Q.row(i)) + Q.row(i).dot(QH.row(i));
    sb(i) = 1.0 / std::sqrt(std::max(bii, 1e-300));
    sc(i) = 1.0 / std::sqrt(std::max(cii, 1e-300));
  }
  Mat2D Lall(N, 4 * r), Rall(N, 4 * r);
  Lall << sb.asDiagonal() * Q, sb.asDiagonal() * W, sc.asDiagonal() * Q, sc.asDiagonal() * dQ;
  Rall << sb.asDiagonal() * W, sb.asDiagonal() * Q, sc.asDiagonal() * (QH - dQ), -(sc.asDiagonal() * Q);
#pragma omp parallel for schedule(static)
  for (Eigen::Index j = 0; j < N; ++j) A.col(j).array() *= sb.array() * sb(j);
  A.noalias() -= Lall * Rall.transpose();
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
