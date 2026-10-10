/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/InbredSamples.cpp
 * @author      Zilong Li
 * Copyright (C) 2025-2026. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#include "InbredSamples.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>

#include "Common.hpp"
#include "InbredSites.hpp"  // open_inbreed_target
#include "Utils.hpp"

// The model. For sample i at site j with individual allele frequency pi, an
// inbreeding coefficient F gives the genotype probabilities
//
//   P(0) = (1-pi)^2 + pi(1-pi) F,   P(1) = 2 pi(1-pi) (1-F),   P(2) = pi^2 + pi(1-pi) F.
//
// The estimator is PCAngsd's --inbreed-samples, the one --inbreed 1 uses per
// site with the roles of sites and samples swapped:
//
//   F_i = 1 - O_i / E_i,   E_i = sum_j 2 pi_ij (1 - pi_ij),
//
// O_i the heterozygotes observed and E_i those expected without inbreeding,
// both summed over the same sites. It is the moment estimator of plink --het,
// with each sample's own pi_ij in place of one frequency per site.
//
//  * called genotypes (-b, -p): O_i counts the heterozygous calls. That is
//    closed form, so a single pass over the data gives F; there is nothing to
//    iterate.
//  * genotype likelihoods (-G): O_i is the posterior expected number of
//    heterozygotes under the prior above, which depends on F_i, so
//    F = 1 - O(F) / E is solved by EM, accelerated by SQUAREM as in PCAngsd.
//    --maxiter and --tol-em apply.
//
// Two differences from PCAngsd, neither of which moves the fixed point:
//
//  * Missing data are left out of O and E. PCAngsd keeps a missing call, or a
//    flat GL, in both sums, where it adds its prior expectation, 2pi(1-pi)(1-F)
//    to O and 2pi(1-pi) to E. At the fixed point that cancels exactly:
//    F (E_obs + E_mis) = E_obs + E_mis - O_obs - (1-F) E_mis gives
//    F = 1 - O_obs / E_obs whatever the missing sites are. Keeping them only
//    slows the EM down, by the fraction missing, and for called genotypes turns
//    a closed form into an iteration. Leaving them out also makes N_SITES in
//    the output the number of sites that carry information.
//  * PCAngsd floors each prior at 1e-4 and renormalises. The renormalisation
//    cancels in the posterior; the floor inflated the homozygote prior of rare
//    variants (pi = 0.001: 1e-6 -> 1e-4) and so the heterozygote posterior of
//    every GL that leans that way. PROB_EPS only keeps a prior positive where
//    F < 0 would make it negative, as in InbredSites.cpp.
//
// F is bounded to [-1, 1] like PCAngsd and --inbreed 1, and not cut at 0: a
// negative F (excess heterozygosity: contamination, a sample mixture, structure
// the PCs miss) is a finding for quality control.

static constexpr double PROB_EPS = 1e-12;  // as in InbredSites.cpp
// a GL whose three values agree to this is taken as no data. ANGSD writes a site
// without reads as 0.333333 0.333333 0.333333, whose third is 0.333334 here
static constexpr double GL_FLAT = 1e-5;

namespace {

// sums over the informative sites of each sample
struct HetSums {
  Mat1D obs;                   // heterozygotes observed (called) or expected a posteriori (GL)
  Mat1D exp;                   // heterozygotes expected at F = 0: sum of 2 pi (1 - pi)
  std::vector<uint64> nsites;  // called genotypes, or GLs that are not flat

  void reset(Eigen::Index n) {
    obs = Mat1D::Zero(n);
    exp = Mat1D::Zero(n);
    nsites.assign(n, 0);
  }
};

using SampleName = std::array<std::string, 2>;  // FID, IID

}  // namespace

// The samples are split into tiles and each tile walks the sites in order, so
// a sample's sums are added up site by site in file order whatever the thread
// count, the tile size or -m: in-core and out-of-core runs, and runs with any
// -n, give the same bits. The samples of a tile are contiguous in each column
// of the column-major PI, G and L.
static Eigen::Index tile_width(Eigen::Index n) {
  const Eigen::Index t = std::max(1, omp_get_max_threads());
  Eigen::Index w = (n + 4 * t - 1) / (4 * t);  // about four tiles per thread
  w = (w + 7) / 8 * 8;                         // whole cache lines of doubles
  return std::min<Eigen::Index>(256, std::max<Eigen::Index>(8, w));
}

// called genotypes: G is N x m on the 0..1 scale of BED2GENO, where anything
// but 0, 0.5 and 1 is a missing call
static void add_calls(HetSums& s, const Mat2D& PI, const Mat2D& G, Eigen::Index m) {
  const Eigen::Index N = PI.rows(), w = tile_width(N), ntiles = (N + w - 1) / w;
#pragma omp parallel for schedule(static)
  for (Eigen::Index t = 0; t < ntiles; ++t) {
    const Eigen::Index i0 = t * w, n = std::min(w, N - i0);
    // this tile's running sums, held locally while it walks the sites
    std::vector<double> o(s.obs.data() + i0, s.obs.data() + i0 + n);
    std::vector<double> e(s.exp.data() + i0, s.exp.data() + i0 + n);
    std::vector<uint64> c(s.nsites.begin() + i0, s.nsites.begin() + i0 + n);
    for (Eigen::Index j = 0; j < m; ++j) {
      const double* pi = PI.col(j).data() + i0;
      const double* g = G.col(j).data() + i0;
      for (Eigen::Index i = 0; i < n; ++i) {
        const double x = g[i];
        if (x != BED2GENO[0] && x != BED2GENO[2] && x != BED2GENO[3]) continue;  // missing
        e[i] += 2.0 * pi[i] * (1.0 - pi[i]);
        if (x == BED2GENO[2]) o[i] += 1.0;
        c[i] += 1;
      }
    }
    std::copy(o.begin(), o.end(), s.obs.data() + i0);
    std::copy(e.begin(), e.end(), s.exp.data() + i0);
    std::copy(c.begin(), c.end(), s.nsites.begin() + i0);
  }
}

// genotype likelihoods: L is 2N x m, rows 2i and 2i+1 the first two GLs of
// sample i, the third being 1 minus both; F the current estimates
static void add_gls(HetSums& s, const Mat1D& F, const Mat2D& PI, const Mat2D& L, Eigen::Index m) {
  const Eigen::Index N = PI.rows(), w = tile_width(N), ntiles = (N + w - 1) / w;
#pragma omp parallel for schedule(static)
  for (Eigen::Index t = 0; t < ntiles; ++t) {
    const Eigen::Index i0 = t * w, n = std::min(w, N - i0);
    // this tile's running sums, held locally while it walks the sites
    std::vector<double> o(s.obs.data() + i0, s.obs.data() + i0 + n);
    std::vector<double> e(s.exp.data() + i0, s.exp.data() + i0 + n);
    std::vector<uint64> c(s.nsites.begin() + i0, s.nsites.begin() + i0 + n);
    const double* f = F.data() + i0;
    for (Eigen::Index j = 0; j < m; ++j) {
      const double* pi = PI.col(j).data() + i0;
      const double* l = L.col(j).data() + 2 * i0;
      for (Eigen::Index i = 0; i < n; ++i) {
        const double l0 = l[2 * i], l1 = l[2 * i + 1], l2 = std::max(0.0, 1.0 - l0 - l1);
        if (std::fabs(l0 - l1) <= GL_FLAT && std::fabs(l1 - l2) <= GL_FLAT) continue;  // no data
        const double p = pi[i], h = p * (1.0 - p);
        // the prior, up to a factor that cancels in the posterior
        const double q0 = l0 * std::max(PROB_EPS, (1.0 - p) * (1.0 - p) + h * f[i]);
        const double q1 = l1 * std::max(PROB_EPS, 2.0 * h * (1.0 - f[i]));
        const double q2 = l2 * std::max(PROB_EPS, p * p + h * f[i]);
        const double qs = q0 + q1 + q2;
        if (!(qs > 0.0) || !std::isfinite(qs)) continue;  // an unreadable GL
        o[i] += q1 / qs;
        e[i] += 2.0 * h;
        c[i] += 1;
      }
    }
    std::copy(o.begin(), o.end(), s.obs.data() + i0);
    std::copy(e.begin(), e.end(), s.exp.data() + i0);
    std::copy(c.begin(), c.end(), s.nsites.begin() + i0);
  }
}

// one pass over all sites; F is used only by the GL model
static void sum_sites(HetSums& s, const Mat1D& F, bool gl, Data* data, Data* Pi, const Param& params) {
  s.reset(Pi->nsamples);
  if (!params.out_of_core) {
    if (gl)
      add_gls(s, F, Pi->G, data->P, Pi->nsnps);
    else
      add_calls(s, Pi->G, data->G, Pi->nsnps);
    return;
  }
  data->check_file_offset_first_var();
  Pi->check_file_offset_first_var();
  for (uint b = 0; b < data->nblocks; b++) {
    data->read_block_initial(data->start[b], data->stop[b], false);
    Pi->read_block_initial(Pi->start[b], Pi->stop[b], false);
    const Eigen::Index m = Pi->stop[b] - Pi->start[b] + 1;
    if (gl)
      add_gls(s, F, Pi->G, data->P, m);
    else
      add_calls(s, Pi->G, data->G, m);
  }
}

// F = 1 - O/E within [-1, 1]. A sample without an informative site keeps 0
// here, so that it cannot upset the SQUAREM step, and is written as NA.
static void update_F(Mat1D& F, const HetSums& s) {
  for (Eigen::Index i = 0; i < F.size(); ++i)
    F(i) = (s.nsites[i] > 0 && s.exp(i) > 0.0) ? std::min(1.0, std::max(-1.0, 1.0 - s.obs(i) / s.exp(i))) : 0.0;
}

// FID and IID of each sample of the target, in its order: the .fam, the .psam
// (its #-header names the columns, and IID alone means FID 0 as in .eigvecs2;
// without a header it is laid out like a .fam), or the BEAGLE header
static std::vector<SampleName> read_sample_names(const Param& params) {
  std::vector<SampleName> ids;
  if (params.file_t == FileType::BEAGLE) {
    for (const auto& id : parse_beagle_samples(params.filein)) ids.push_back({id, id});
    return ids;
  }
  const std::string path = params.filein + (params.file_t == FileType::PGEN ? ".psam" : ".fam");
  std::ifstream ifs(path);
  if (!ifs.is_open()) cao.error("can not open " + path);
  int cfid = 0, ciid = 1;
  std::string line, tok;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    std::vector<std::string> t;
    std::istringstream iss(line[0] == '#' ? line.substr(1) : line);
    while (iss >> tok) t.push_back(tok);
    if (line[0] == '#') {
      if (line.size() > 1 && line[1] == '#') continue;  // ## lines of a .psam
      cfid = ciid = -1;
      for (int c = 0; c < (int)t.size(); ++c) {
        if (t[c] == "FID") cfid = c;
        if (t[c] == "IID") ciid = c;
      }
      if (ciid < 0) cao.error("the header line of " + path + " has no IID column");
      continue;
    }
    if ((int)t.size() <= std::max(cfid, ciid)) cao.error("too few columns in " + path + " at: " + line);
    ids.push_back({cfid >= 0 ? t[cfid] : std::string("0"), t[ciid]});
  }
  return ids;
}

static void write_inbred_samples(const std::string& fout,
                                 const std::vector<SampleName>& ids,
                                 const HetSums& s,
                                 const Mat1D& F) {
  std::ofstream ofs(fout);
  if (!ofs.is_open()) cao.error("can not open " + fout);
  ofs << "#FID\tIID\tN_SITES\tO_HET\tE_HET\tF\n" << std::setprecision(8);
  for (Eigen::Index i = 0; i < F.size(); ++i) {
    ofs << ids[i][0] << '\t' << ids[i][1] << '\t' << s.nsites[i] << '\t' << s.obs(i) << '\t' << s.exp(i) << '\t';
    if (s.nsites[i] > 0)
      ofs << F(i) << '\n';
    else
      ofs << "NA\n";
  }
  if (!ofs) cao.error("error writing " + fout);
}

// SQUAREM on the EM map T(F) = 1 - O(F)/E, as PCAngsd and --inbreed 1
static void solve_gl(Mat1D& F, HetSums& s, Data* data, Data* Pi, const Param& params) {
  const Eigen::Index N = F.size();
  Mat1D F0(N), F1(N), F2(N);
  const uint maxiter = std::max(1u, params.maxiter);  // --maxiter 0 would leave s empty
  for (uint it = 0; it < maxiter; it++) {
    F0 = F;
    sum_sites(s, F0, true, data, Pi, params);
    F1 = F0;
    update_F(F1, s);
    sum_sites(s, F1, true, data, Pi, params);
    F2 = F1;
    update_F(F2, s);
    const Mat1D D1 = F1 - F0, D2 = F2 - F1;
    const double sr2 = D1.squaredNorm(), sv2 = (D2 - D1).squaredNorm();
    double diff;
    if (sv2 == 0.0) {
      // T(F1) - F1 = F1 - F0: no curvature to extrapolate from, take the EM step
      F = F2;
      diff = rmse1d(F1, F2);
    } else {
      const double alpha = std::min(256.0, std::max(1.0, std::sqrt(sr2 / sv2)));
      if (params.verbose > 1) cao.print("alpha:", alpha, ", sr2:", sr2, ", sv2:", sv2);
      F = (F0 + 2.0 * alpha * D1 + alpha * alpha * (D2 - D1)).cwiseMax(-1.0).cwiseMin(1.0);
      // stabilization step; its sums are the ones written out, so the output
      // satisfies F = 1 - O_HET/E_HET exactly
      sum_sites(s, F, true, data, Pi, params);
      update_F(F, s);
      diff = rmse1d(F0, F);
    }
    cao.print(tick.date(), "Inbreeding coefficients estimated, iter =", it + 1, ", RMSE =", diff);
    if (diff < params.tolem) {
      cao.print(tick.date(), "EM inbreeding coefficient converged");
      return;
    }
  }
  cao.warn("EM inbreeding coefficient not converged in " + std::to_string(maxiter) + " iterations! raise --maxiter");
}

void run_inbred_samples(Data* Pi, const Param& params) {
  // before any work: a malformed sample file should not cost an EM run
  const std::vector<SampleName> ids = read_sample_names(params);
  std::unique_ptr<Data> data(open_inbreed_target(Pi, params));
  const Eigen::Index N = Pi->nsamples;
  if ((Eigen::Index)ids.size() != N)
    cao.error("found " + std::to_string(ids.size()) + " sample IDs for the " + std::to_string(N) +
              " samples of the input");

  const bool gl = params.file_t == FileType::BEAGLE;
  cao.print(tick.date(), "run inbreeding coefficient estimator per sample");
  HetSums s;
  Mat1D F = Mat1D::Zero(N);
  if (gl) {
    solve_gl(F, s, data.get(), Pi, params);
  } else {
    // called genotypes: closed form, one pass
    sum_sites(s, F, false, data.get(), Pi, params);
    update_F(F, s);
  }

  uint64 nempty = 0, nneg = 0;
  double fmin = std::numeric_limits<double>::infinity(), fmax = -fmin, fsum = 0.0;
  for (Eigen::Index i = 0; i < N; ++i) {
    if (s.nsites[i] == 0) {
      ++nempty;
      continue;
    }
    fsum += F(i);
    fmin = std::min(fmin, F(i));
    fmax = std::max(fmax, F(i));
    if (F(i) < 0.0) ++nneg;
  }
  if (nempty > 0)
    cao.warn(std::to_string(nempty) + " sample(s) have no " +
             (gl ? "informative genotype likelihood" : "called genotype") + " at any site; their F is NA");
  if (nempty < (uint64)N)
    cao.print(tick.date(), "per-sample F: mean =", fsum / (double)(N - nempty), ", min =", fmin, ", max =", fmax, ",",
              nneg, "sample(s) with F < 0");
  write_inbred_samples(params.fileout + ".inbred", ids, s, F);
  cao.print(tick.date(), "per-sample inbreeding coefficients saved to", params.fileout + ".inbred");
}
