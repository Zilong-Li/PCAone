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
 * Memory is O(N^2), independent of the number of sites.
 *
 * Missing genotypes are imputed to the site mean, which is what keeps the
 * projection -- taken across individuals within a site -- well defined. The
 * identity above then holds exactly for the imputed matrix.
 *
 * Note PCAone codes genotypes on the 0..1 scale (BED2GENO = {1, NA, 0.5, 0}),
 * i.e. x = g/2. Both bhat and chat are correlation matrices, so the constant
 * factor between the x and g scales cancels and no rescaling is needed.
 *
 * The PC scores are read with Utils::read_usv(), which this branch also fixes:
 * it used to map a row-major buffer as column-major, transposing any file with
 * more than one column. --evaladmix-k is a direct regression test for that.
 ******************************************************************************/
#include "EvalAdmix.hpp"

#include <fstream>
#include <iomanip>
#include <new>
#include <sstream>
#include <string>
#include <vector>

#include "Utils.hpp"

// turn a symmetric PSD matrix into a correlation matrix, in place -- these are
// N x N, so returning a copy would cost another one of them
static void cov2cor(Mat2D& C) {
  const Eigen::Index n = C.rows();
  Mat1D inv(n);
  for (Eigen::Index i = 0; i < n; ++i) inv(i) = 1.0 / std::sqrt(std::max(C(i, i), 1e-300));
  for (Eigen::Index j = 0; j < n; ++j)  // column-major order
    for (Eigen::Index i = 0; i < n; ++i) C(i, j) *= inv(i) * inv(j);
}

static void write_matrix(const std::string& fn, const Mat2D& M, const std::vector<std::string>& ids,
                         double scale = 1.0) {
  std::ofstream ofs(fn);
  if (!ofs.is_open()) cao.error("can not open file for writing: " + fn);
  ofs << std::fixed << std::setprecision(6);
  if (!ids.empty()) {
    for (size_t i = 0; i < ids.size(); ++i) ofs << (i ? "\t" : "") << ids[i];
    ofs << "\n";
  }
  for (Eigen::Index i = 0; i < M.rows(); ++i) {
    for (Eigen::Index j = 0; j < M.cols(); ++j) ofs << (j ? "\t" : "") << M(i, j) * scale;
    ofs << "\n";
  }
}

// bytes of one N x N double matrix, in GiB
static double nn_gib(Eigen::Index N) { return (double)N * (double)N * 8.0 / 1073741824.0; }

void run_evaladmix(Data* data, const Param& params) {
  const Eigen::Index N = data->nsamples;
  const uint M = data->nsnps;

  // ---- 0. size -----------------------------------------------------------
  // Everything below is dense N x N: four such matrices at peak, and two dense
  // N x N text files on the way out. Both grow with the square of the sample
  // count and neither depends on the number of sites, so a cohort whose PCA
  // runs comfortably out-of-core can still be far out of reach here. Say so
  // before spending a pass over the genotypes rather than dying in the
  // allocator afterwards.
  const double ram = 4.0 * nn_gib(N);
  const double perfile = (double)N * (double)N * 9.0 / 1073741824.0;  // ~9 bytes per printed value
  cao.print(tick.date(), "evalAdmix:", N, "samples needs about", ram, "GB of RAM, and writes two files of about",
            perfile, "GB each");
  if (ram > 8.0)
    cao.warn("evalAdmix needs about ", ram, " GB of RAM for the ", N, " x ", N,
             " matrices. this does not go down with --memory, which only bounds the genotype blocks. reduce the "
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
    cao.error("evalAdmix: ", fpcs, " has ", U.rows(), " rows but the genotype file has ", N, " samples");
  Eigen::Index k = params.evaladmix_k > 0 ? params.evaladmix_k : U.cols();
  if (k > U.cols())
    cao.error("--evaladmix-k is ", k, " but only ", U.cols(), " PCs were computed; raise -k or lower --evaladmix-k");
  cao.print(tick.date(), "evalAdmix: using", k, "PC(s) + intercept =", k + 1, "dimensions");

  // ---- 2. one streaming pass: Gram matrix, mean genotype, heterozygosity --
  //
  // Both branches below see *centred* genotypes (x - f) with missing calls
  // imputed to the site mean (centred value 0). That is deliberate:
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
  // Mean imputation keeps the projection (I-P), which mixes individuals within
  // a site, well defined when a genotype is absent. It attenuates the residual
  // of a missing call towards zero, so heavily missing samples get a slightly
  // shrunk statistic; pairwise-complete counts would avoid that but cannot be
  // folded into the single O(N^2) pass.
  Mat2D A;
  try {
    A = Mat2D::Zero(N, N);  // G'G
  } catch (const std::bad_alloc&) {
    cao.error("evalAdmix: out of memory allocating the ", N, " x ", N, " Gram matrix (", nn_gib(N),
              " GB); the whole statistic needs about ", ram, " GB");
  }
  Mat1D b = Mat1D::Zero(N);  // sum_s g_s
  Mat1D d = Mat1D::Zero(N);  // sum_s g_s .* (1 - g_s)   (0..1 scale)
  tick.clock();
  if (!params.out_of_core) {
    // G is nsamples x nsnps, centred, not standardized.
    const Mat2D& G = data->G;
    A.noalias() = G * G.transpose();
    b = G.rowwise().sum();
    for (Eigen::Index i = 0; i < G.cols(); ++i) {
      const double f = data->F(i);
      d.array() += (G.col(i).array() + f) * (1.0 - (G.col(i).array() + f));
    }
  } else {
    // Out-of-core. read_block_initial() estimates F on the fly, so F and the
    // lookup table must be allocated here before the first block.
    data->F = Mat1D::Zero(data->nsnps);
    data->centered_geno_lookup = Arr2D::Zero(4, data->nsnps);
    data->check_file_offset_first_var();
    for (uint bi = 0; bi < data->nblocks; ++bi) {
      data->read_block_initial(data->start[bi], data->stop[bi], false);
      const Mat2D& G = data->G;
      A.noalias() += G * G.transpose();
      b += G.rowwise().sum();
      for (Eigen::Index i = 0; i < G.cols(); ++i) {
        const double f = data->F(data->start[bi] + i);
        d.array() += (G.col(i).array() + f) * (1.0 - (G.col(i).array() + f));
      }
    }
  }
  b /= (double)M;
  d /= (double)M;
  cao.print(tick.date(), "evalAdmix: accumulated summary statistics over", M, "sites in", tick.reltime(),
            "seconds");
  cao.print(tick.date(),
            "evalAdmix: missing genotypes, if any, were imputed to the site mean; high missingness shrinks the "
            "statistic towards zero");

  // ---- 3. projection onto [PCs, intercept] -------------------------------
  Mat2D V(N, k + 1);
  Mat2D IP;
  {
    Mat2D V(N, k + 1);
    V.leftCols(k) = U.leftCols(k);
    V.col(k).setOnes();
    const Mat2D VtV = V.transpose() * V;
    IP.noalias() = V * VtV.completeOrthogonalDecomposition().pseudoInverse() * V.transpose();
  }
  IP = -IP;                      // IP = I - P, without materialising I
  IP.diagonal().array() += 1.0;

  // ---- 4. bhat and chat --------------------------------------------------
  // Written to hold four N x N matrices at peak rather than eight: A is reused
  // for every intermediate and ends up holding corres, and chat is folded in as
  // soon as it exists. At N = 10,000 that is the difference between 3 GB and
  // 6 GB.
  try {
    Mat2D W(N, N), C(N, N);
    // chat = cov2cor( (I-P) Dhat (I-P) ), Dhat diagonal so the left product is
    // a column scaling rather than a matrix multiply
    W = IP * d.asDiagonal();
    C.noalias() = W * IP;
    cov2cor(C);
    // bhat = cov2cor( (I-P) [A - M gbar gbar'] (I-P) )
    A.noalias() -= (double)M * b * b.transpose();
    W.noalias() = IP * A;
    A.noalias() = W * IP;
    cov2cor(A);
    A -= C;  // corres = bhat - chat
  } catch (const std::bad_alloc&) {
    cao.error("evalAdmix: out of memory; the ", N, " x ", N, " matrices need about ", ram, " GB in total");
  }
  IP.resize(0, 0);
  // a correlation cannot leave [-1, 1]
  A = A.cwiseMax(-1.0).cwiseMin(1.0);
  A.diagonal().setZero();

  // ---- 5. output ---------------------------------------------------------
  std::vector<std::string> ids;
  if (params.file_t == FileType::PLINK) {
    std::ifstream ffam(params.filein + ".fam");
    std::string line, fid, iid;
    while (std::getline(ffam, line)) {
      std::istringstream iss(line);
      if (iss >> fid >> iid) ids.push_back(iid);
    }
  }
  if ((Eigen::Index)ids.size() != N) ids.clear();
  write_matrix(params.fileout + ".corres", A, ids);
  // the correlation of residuals estimates 2*phi, so halve it for kinship. done
  // at write time rather than in another N x N copy; the diagonal is already 0.
  write_matrix(params.fileout + ".kinship", A, ids, 0.5);
  cao.print(tick.date(), "evalAdmix: correlation of residuals saved to", params.fileout + ".corres");
  cao.print(tick.date(), "evalAdmix: kinship (corres/2) saved to", params.fileout + ".kinship");
}
