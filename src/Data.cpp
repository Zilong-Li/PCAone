/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Data.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Data.hpp"

#include <omp.h>

#include "Cmd.hpp"
#include "LinAlg.hpp"
#include "TextIO.hpp"
#include "Utils.hpp"
#include "Variants.hpp"

using namespace std;

void Data::prepare() {
  if (nsamples > nsnps) nsamples_ge_nsnps = true;

  if (!params.dopca) {  // for projection, read F from reference set
    cao.print(tick.date(), "read allele frequency from .mbim file: " + params.filebim);
    F = read_frq(params.filebim);
    if (!keepRefSNPs.empty()) {
      Mat1D Fnew(keepRefSNPs.size());
      for (int i = 0; i < (int)keepRefSNPs.size(); ++i) {
        Fnew(i) = F(keepRefSNPs[i]);
      }
      F = Fnew;
      cao.print(tick.date(), "keep AF of", F.size(), " matched sites in the reference");
    }
    if (!flipSNPs.empty()) {
      for (int i : flipSNPs) F(i) = 1.0 - F(i);
    }
  }

  if (!params.out_of_core) {
    tick.clock();
    read_all();
    readtime += tick.reltime();
    cao.print(tick.date(), "done reading all data");
    return;
  }

  if (params.project > 0) {
    // out-of-core --project 1/2 reads only the sites matched to the reference,
    // in .bim/.pvar order (keepSNPs; all of them when identical), and sizes its
    // blocks itself from -m, K and N (run_projection_blockwise -> set_blocks)
    if (!keepSNPs.empty()) {
      nsnps_all = nsnps;
      nsnps = keepSNPs.size();
    }
    return;
  }

  // LD and evalAdmix keep dopca on to estimate F, but run no PCA.
  // Size their blocks like any other non-PCA run.
  const bool pca_blocks = params.dopca && !params.ld && !params.evaladmix;
  const bool exact = pca_blocks && params.svd_t == SvdType::FULL;
  if (exact) {
    // exact PCA (Exact.cpp): the N x N GRM and F, plus the loadings unless --no-loadings,
    // are held besides one block. Without -m, blocks of ~64 MB (and at least
    // 2048 sites) keep the GRM update a large matrix product.
    const double n = std::max(1u, nsamples);
    const double fixed = n * n + 6.0 * nsnps + (params.printv ? (double)nsnps * params.k : 0.0);
    double b = std::max(2048.0, 8388608.0 / n);
    if (params.memory > 0) {
      if (params.memory * 134217728 > 1.1 * fixed + 64 * n)
        b = (params.memory * 134217728 - fixed) / n;
      else
        cao.warn("the exact PCA needs at least", fixed / 134217728, " GB for the", nsamples, " x", nsamples,
                 " GRM. -m is too small, using blocks of", (uint64)std::ceil(b), " sites");
    }
    blocksize = (uint)std::min<double>(std::max(1.0, std::ceil(b)), std::max(1u, nsnps));
  } else if (pca_blocks) {
    if (params.svd_t == SvdType::IRAM) {
      // ram of arnoldi = n * b * 8 / 1024 kb
      blocksize = (uint)ceil((double)params.memory * 134217728 / nsamples);
    } else {
      // ram of halko = (3*n*l + 2*m*l + 5*m + n*b)*8/1024 Kb
      // in doubles: the unsigned products wrapped once nsnps * l reached 2^31
      // the RSVD of the EM iterations (--em-k) and of the final matrix (-k) run one at a time
      const double l = std::max((double)params.k + params.oversamples, (double)params.em_k + params.em_oversamples);
      const double fixed = 3.0 * nsamples * l + 2.0 * nsnps * l + 5.0 * nsnps;  // doubles held besides the block
      double m = fixed / 134217728;
      if (params.memory > 1.1 * m)
        m = 0;
      else
        cao.warn("minimum RAM required is ", m, " GB. trying to allocate more RAM.");
      blocksize = (unsigned int)ceil(((m + params.memory) * 134217728 - fixed) / nsamples);
    }
  } else if (params.evaladmix_pairs) {
    // --evaladmix-kin (EvalAdmix.cpp): the stripes of the N x N matrix take
    // most of -m. A block only has to be wide enough for an efficient product
    // with them: a fifth of -m, between 256 and 2048 sites. More would only take
    // room from the stripes, i.e. add passes over the genotypes.
    const double b = 0.2 * params.memory * 134217728 / std::max(1u, nsamples);
    blocksize = (uint)std::min<double>(std::max(1u, nsnps), std::max(256.0, std::min(2048.0, std::floor(b))));
  } else if (params.inbreed) {
    // --inbreed (InbredSites.cpp, InbredSamples.cpp) holds a block of the target
    // and the matching block of individual allele frequencies (N rows) at once.
    // The target is N rows of genotypes, or 2N rows of likelihoods for BEAGLE, so
    // -m is split over 2 or 3 such N-row blocks. Both Data objects size their
    // blocks here, alike.
    const double parts = params.file_t == FileType::BEAGLE ? 3.0 : 2.0;
    blocksize = (uint)ceil((double)params.memory / parts * 134217728 / nsamples);
  } else {
    // ram of non-pca run
    blocksize = (uint)ceil((double)params.memory * 134217728 / nsamples);
  }

  nblocks = (unsigned int)ceil((double)nsnps / blocksize);
  if (exact) {
    cao.print(tick.date(), "blocks for the exact PCA: blocksize =", blocksize, ", nblocks =", nblocks);
  } else {
    cao.print(tick.date(), "initial setting by -m/--memory: blocksize =", blocksize, ", nblocks =", nblocks,
              ", factor =", bandFactor);
    // the passes of --evaladmix-kin read the blocks again, however many there are
    if (nblocks == 1 && !params.evaladmix_pairs) cao.error("only one block exists. please remove -m option");
  }
  if (pca_blocks && params.svd_t == SvdType::PCAoneAlg2) {
    // decrease blocksize for the winSVD
    if (nblocks < params.bands) {
      blocksize = (unsigned int)ceil((double)nsnps / params.bands);
    } else {
      bandFactor = (unsigned int)ceil((double)nblocks / params.bands);
      blocksize = (unsigned int)ceil((double)nsnps / (params.bands * bandFactor));
    }
    nblocks = (unsigned int)ceil((double)nsnps / blocksize);
    cao.print(tick.date(), "after adjustment by PCAone: -w =", params.bands, ", blocksize =", blocksize,
              ", nblocks =", nblocks, ", factor =", bandFactor);
  }
  start.resize(nblocks);
  stop.resize(nblocks);
  for (uint i = 0; i < nblocks; i++) {
    start[i] = i * blocksize;
    stop[i] = start[i] + blocksize - 1;
    stop[i] = stop[i] >= nsnps ? nsnps - 1 : stop[i];
  }
}

void Data::set_blocks(uint b) {
  blocksize = std::max(1u, std::min(b, std::max(1u, nsnps)));
  nblocks = (nsnps + blocksize - 1) / blocksize;
  start.resize(nblocks);
  stop.resize(nblocks);
  for (uint i = 0; i < nblocks; i++) {
    start[i] = i * blocksize;
    stop[i] = std::min<uint64>((uint64)start[i] + blocksize, nsnps) - 1;
  }
}

// filter snps, update keepSNPs, reassign nsnps;
void Data::filter_snps_resize_F() {
  if (!params.filterSNP) return;
  if (!(params.maf > 0 && params.maf <= 0.5)) cao.error("--maf has to be between (0, 0.5)");

  nsnps_all = F.size();
  Mat1D Fnew(F.size());  // make a temp F
  int i, j;
  for (i = 0, j = 0; j < (int)F.size(); j++) {
    if (MAF(F(j)) > params.maf) {
      keepSNPs.push_back(j);  // keep track of index of element > maf
      Fnew(i++) = F(j);
    }
  }
  nsnps = keepSNPs.size();  // new number of SNPs
  cao.print(tick.date(), "number of SNPs after filtering by MAF >", params.maf, ":", nsnps);
  if (nsnps < 1) cao.error("no SNPs left after filtering!");
  // resize F
  F.noalias() = Fnew.head(nsnps);
}

// The variants of the .bim/.pvar, each with its allele frequency, in the order
// of the input; not written for CSV and BEAGLE, which have no .bim
void Data::save_snps_in_mbim() {
  const bool pvar = params.file_t == FileType::PGEN;
  const std::string bim_path = params.filein + (pvar ? ".pvar" : ".bim");
  std::ifstream ifs_bim(bim_path);
  if (!ifs_bim.is_open()) {
    // CSV has no variant metadata; now that .loadings are written by default,
    // warning on every CSV run would be noise
    if (params.file_t != FileType::CSV) cao.warn(params.filein + ".bim/.pvar not found; skipping mbim output");
    return;
  }
  const std::string fout = params.fileout + ".mbim";
  std::ofstream ofs_bim(fout, std::ios::binary);
  if (!ofs_bim.is_open()) cao.error("can not open file for writing: " + fout);
  // a line of the .mbim: the .bim line and the frequency
  auto put_line = [](TextBuf& buf, const std::string& meta, double f) {
    buf.append(meta);
    char* p = buf.reserve(TEXT_NUM_MAX + 2);
    *p++ = '\t';
    p = put_g6(f, p);
    *p++ = '\n';
    buf.commit(p);
  };
  std::string line;
  const bool metadata_is_permuted = metadata_permuted && perm.indices().size() == nsnps;
  const bool frequency_is_permuted = params.perm && params.out_of_core && perm.indices().size() == nsnps;

  // the logical site whose frequency goes with each line
  std::vector<Eigen::Index> original_to_logical(nsnps);
  for (Eigen::Index logical_idx = 0; logical_idx < nsnps; ++logical_idx) {
    Eigen::Index original_idx = frequency_is_permuted ? perm.indices()[logical_idx] : logical_idx;
    original_to_logical[original_idx] = logical_idx;
  }

  if (!pvar) {
    // the .bim at once, its lines found in the buffer, and the .mbim written
    // from it on every thread: reading it a line at a time into strings took
    // as long as an epoch of the winSVD
    ifs_bim.seekg(0, std::ios::end);
    std::string text((size_t)ifs_bim.tellg(), '\0');
    ifs_bim.seekg(0);
    ifs_bim.read(text.data(), (std::streamsize)text.size());
    if (!ifs_bim) cao.error("cannot read " + bim_path);
    std::vector<size_t> start{0};  // the lines as std::getline splits them
    for (const char *p = text.data(), *e = p + text.size(); (p = (const char*)std::memchr(p, '\n', e - p)); ++p)
      start.push_back((size_t)(p - text.data()) + 1);
    if (start.back() < text.size()) start.push_back(text.size() + 1);  // a last line without '\n'
    const Eigen::Index nlines = (Eigen::Index)start.size() - 1;
    // the source line of each site, by the order the metadata was read in
    std::vector<Eigen::Index> src(nsnps, -1);
    for (Eigen::Index j = 0; j < nsnps; ++j) {
      Eigen::Index s = j;
      if (metadata_is_permuted) {
        src[perm.indices()[j]] = j;
        continue;
      }
      if (params.filterSNP) s = keepSNPs[j];
      src[j] = s;
    }
    // with the streamed read, only as many lines as F has
    const Eigen::Index nout =
        (!params.filterSNP && !params.perm) ? std::min<Eigen::Index>(nsnps, nlines) : (Eigen::Index)nsnps;
    write_lines(ofs_bim, nout, 64, [&](Eigen::Index i, TextBuf& buf) {
      const Eigen::Index s = src[i];
      if (s >= 0 && s < nlines) buf.append(text.data() + start[s], start[s + 1] - 1 - start[s]);
      char* p = buf.reserve(TEXT_NUM_MAX + 2);
      *p++ = '\t';
      p = put_g6(F(original_to_logical[i]), p);
      *p++ = '\n';
      buf.commit(p);
    });
    ofs_bim.close();
    if (!ofs_bim) cao.error("error writing " + fout);
    if (params.filterSNP || params.perm)
      cao.print(tick.date(), "save matched sites in .mbim file and permutation mode is", params.perm);
    return;
  }

  if (!params.filterSNP && !params.perm) {  // the variants are in the order of F: stream them
    TextBuf buf;
    for (Eigen::Index j = 0; j < F.size() && read_bim_line(ifs_bim, pvar, bim_path, line); ++j) {
      put_line(buf, line, F(j));
      if (buf.size() >= (size_t(32) << 20)) buf.write_to(ofs_bim), buf.clear();
    }
    buf.write_to(ofs_bim);
    ofs_bim.close();
    if (!ofs_bim) cao.error("error writing " + fout);
    return;
  }

  std::vector<std::string> metadata(nsnps);
  if (metadata_is_permuted) {
    for (Eigen::Index permuted_idx = 0; permuted_idx < nsnps && read_bim_line(ifs_bim, pvar, bim_path, line);
         ++permuted_idx) {
      Eigen::Index original_idx = perm.indices()[permuted_idx];
      metadata[original_idx] = line;
    }
  } else if (params.filterSNP) {
    int kept = 0;
    for (int source_idx = 0; kept < (int)nsnps && read_bim_line(ifs_bim, pvar, bim_path, line); ++source_idx) {
      if (kept < (int)keepSNPs.size() && keepSNPs[kept] == source_idx) {
        metadata[kept++] = line;
      }
    }
  } else {
    for (Eigen::Index j = 0; j < nsnps && read_bim_line(ifs_bim, pvar, bim_path, line); ++j) {
      metadata[j] = line;
    }
  }

  write_lines(ofs_bim, nsnps, 64,
              [&](Eigen::Index i, TextBuf& buf) { put_line(buf, metadata[i], F(original_to_logical[i])); });
  ofs_bim.close();
  if (!ofs_bim) cao.error("error writing " + fout);
  cao.print(tick.date(), "save matched sites in .mbim file and permutation mode is", params.perm);
}

/**
  T = U'/s
  VT = (U'/s) * G = T * G
  V = G' * (U/s) // calculate V is not a good idea
 **/
void Data::calcu_vt_initial(const Mat2D& T, Mat2D& VT, bool standardize) {
  if (nblocks == 1) {
    cao.error("only one block exists. please use in-memory mode instead");
  }
  uint actual_block_size;
  check_file_offset_first_var();
  for (uint i = 0; i < nblocks; ++i) {
    actual_block_size = stop[i] - start[i] + 1;
    // G (nsamples, actual_block_size)
    read_block_initial(start[i], stop[i], standardize);
    VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
  }

  return;
}

void Data::calcu_vt_update(
    const Mat2D& T, const Mat2D& U, const Mat1D& svals, const Mat2D& VT, Mat2D& VTout, bool standardize) {
  if (nblocks == 1) {
    cao.error("only one block exists. please use in-memory mode instead");
  }
  uint actual_block_size;
  check_file_offset_first_var();
  for (uint i = 0; i < nblocks; ++i) {
    actual_block_size = stop[i] - start[i] + 1;
    // G (nsamples, actual_block_size)
    read_block_update(start[i], stop[i], U, svals, VT, standardize);
    VTout.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
  }

  return;
}

// S: signular values
// E: eigen values
void Data::set_svd_transform(bool standardized) {
  // standardize_E() is a no-op unless params.scale is the genetic default, so a
  // decomposition that "standardized" under any other --scale did not actually
  // get the sqrt(ploidy)/sd transform that -P knows how to undo.
  svd_scale = (standardized && params.scale == SCALE_STANDARDIZE_GENETIC) ? SCALE_STANDARDIZE_GENETIC : params.scale;
  if (!standardized && params.scale == SCALE_STANDARDIZE_GENETIC) svd_scale = 0;  // centred only
  // fit_with_pi() builds the pcangsd matrix as dosage - 2f, on the 0..2 scale,
  // rather than PCAone's usual 0..1 coding.
  svd_gscale = params.pcangsd ? 2 : 1;
}

void Data::write_eigs_files(const Mat1D& E, const Mat1D& S, const Mat2D& U, const Mat2D& V) {
  // key=value fields are appended after nsamples,nsnps; older PCAone parses
  // the first two with stoi and ignores the rest. See UsvTransform.
  const std::string sig = "#" + std::to_string(U.rows()) + "," + std::to_string(V.rows()) +
                          ",scale=" + std::to_string(svd_scale) + ",ploidy=" + std::to_string(params.ploidy) +
                          ",gscale=" + std::to_string(svd_gscale) + "\n";
  write_matrix(params.fileout + ".sigvals", S, {}, sig);
  write_matrix(params.fileout + ".eigvals", E);
  write_matrix(params.fileout + ".eigvecs", U);
  // the PCs again, with the FID and IID of each sample, as plink2's .eigenvec
  if (params.file_t == FileType::PLINK || params.file_t == FileType::PGEN)
    write_eigvecs2(params.fileout + ".eigvecs2", U, read_sample_ids(params));
  if (params.printv) {
    save_snps_in_mbim();
    // the loadings of a shuffled run are put back in the order of the input
    std::vector<Eigen::Index> order;
    if (params.perm && V.rows() == nsnps && perm.indices().size() == nsnps) {
      order.resize(nsnps);
      for (Eigen::Index logical_idx = 0; logical_idx < nsnps; ++logical_idx)
        order[perm.indices()[logical_idx]] = logical_idx;
    }
    write_matrix(params.fileout + ".loadings", V, {}, "", order.empty() ? nullptr : &order);
  }

  cao.print(tick.date(), "eigen vectors and values saved");
}

void Data::write_pcangsd_cov() const {
  Mat2D C = G * G.transpose();
  C.array() /= (double)nsnps;
  C.diagonal() = Dc.array() / (double)nsnps;
  TextFormat fmt;  // Eigen's default `os << C`, as it has always been written
  fmt.sep = ' ';
  fmt.align = true;
  write_matrix(params.fileout + ".cov", C, fmt);
  // C is symmetric: its eigenvectors are its singular vectors. JacobiSVD took
  // 21 s at N = 1000 and 215 s at N = 1500; this takes 0.7 s and 2.5 s.
  Eigen::SelfAdjointEigenSolver<Mat2D> es(C);
  if (es.info() != Eigen::Success) cao.error("failed eigendecomposition of the PCAngsd covariance matrix");
  const Eigen::Index k = std::min<Eigen::Index>(params.k, C.rows());
  Mat2D U = es.eigenvectors().rightCols(k).rowwise().reverse();  // largest first
  for (Eigen::Index j = 0; j < k; ++j) {                         // same sign rule as flip_UV
    Eigen::Index x;
    U.col(j).cwiseAbs().maxCoeff(&x);
    if (U(x, j) < 0) U.col(j) *= -1;
  }
  write_eigvecs2(params.fileout + ".eigvecs2", U, read_sample_ids(params));
}

void Data::fit_with_pi(const Mat2D& U, const Mat1D& svals, const Mat2D& VT) {
  if (params.verbose >= 3) cao.print(tick.date(), "call fit_with_pi");
  uint ks = svals.size();
  if (params.pcangsd) {  // for pcangsd with beagle input
    const Mat2D US = U * svals.asDiagonal();
    for_each_product_column(US, VT, 0, nsnps, [&](Eigen::Index jj, const auto& recon) {
      const uint j = (uint)jj;
      double p0, p1, p2;
      const uint original = unpermuted_snp_index(j);
      const double f = F(original);
      uint s = params.filterSNP ? keepSNPs[original] : original;
      for (uint i = 0; i < nsamples; ++i) {
        // Rescale individual allele frequencies
        double pt = recon(i);
        pt = (pt + 2.0 * f) / 2.0;
        pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
        // update E, which is G here
        p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        G(i, j) = (p1 + 2.0 * p2) / (p0 + p1 + p2) - 2.0 * f;
      }
    });
  }

  if (params.emu) {
#pragma omp parallel for
    for (uint i = 0; i < nsnps; ++i) {
      const uint original = unpermuted_snp_index(i);
      const double f = F(original);
      for (uint j = 0; j < nsamples; ++j) {
        if (C[(uint64)original * nsamples + j]) {  // sites need to be predicted
          G(j, i) = 0.0;
          for (uint k = 0; k < ks; ++k) {
            G(j, i) += U(j, k) * svals(k) * VT(k, i);
          }
          G(j, i) = fmin(fmax(G(j, i), -f), 1 - f);
        }
      }
    }
  }
}

bool Data::packed_pca() const {
  return params.file_t == FileType::PLINK && params.dopca && !params.nopacked &&
         (params.svd_t == SvdType::PCAoneAlg1 || params.svd_t == SvdType::PCAoneAlg2) && !params.missme &&
         !params.pcangsd && params.center && !params.ld && !params.evaladmix && params.project == 0 &&
         params.selection == 0 && params.inbreed == 0 &&
         (params.scale == SCALE_STANDARDIZE_GENETIC || params.scale == 0);
}

void Data::set_packed_sites(bool standardize, uint64 j0, uint64 n) {
  if (!packed) {  // out-of-core: the records come with each block
    const PackedCoding coding = packed_coding();
    packed = std::make_unique<PCAone::PackedGenotypes>(nullptr, nsamples, nsnps, coding.value, coding.missing);
  }
  if (n == (uint64)-1) n = nsnps - j0;
  Mat1D a(n), b(n);
#pragma omp parallel for
  for (uint64 i = 0; i < n; ++i) {
    const uint64 j = j0 + i;
    const double f = F(unpermuted_snp_index(j));
    double s = 1.0;
    if (standardize && params.scale == SCALE_STANDARDIZE_GENETIC) {
      const double sd = sqrt(f * (1.0 - f));
      if (sd > VAR_TOL) s = sqrt((double)params.ploidy) / sd;
    }
    a(i) = -s * f, b(i) = s;
  }
  packed->set_sites(j0, a, b, Mat1D::Zero(n));
}

void Data::standardize_E() {
  if (params.verbose >= 3) cao.print(tick.date(), "standardize the matrix");
  if (params.scale != -9) return;
#pragma omp parallel for
  for (uint i = 0; i < nsnps; ++i) {
    const double f = F(unpermuted_snp_index(i));
    const double sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) {  // in case denominator is too small.
      G.col(i) *= sqrt((double)params.ploidy) / sd;
    }
  }
}

// Decide how to scale a target genotype matrix so it matches the matrix a
// reference PCA decomposed.
//
// Before the transform was recorded in .sigvals, projection and selection
// scaled their G from this run's --scale, which silently mismatched whenever
// the reference had used anything else. That included the case where the user
// passed no --scale at all: a reference built with -D/--ld has standardisation
// turned off, and -D is a documented way to get the .mbim that -P needs.
bool Data::resolve_ref_scaling(const UsvTransform& t, const std::string& src, bool allow_dosage) const {
  if (!t.known) {
    cao.warn(src,
             " predates the recording of the PCA scaling, so the reference is assumed to have used this run's "
             "-C/--scale. that is wrong if the reference was built with -D/--ld, or with a non-default --scale; "
             "rerun it to record the transform.");
    return params.scale == SCALE_STANDARDIZE_GENETIC;
  }
  if (t.scale > 0)
    cao.error("the reference PCA used -C/--scale ", t.scale, ", which cannot be replayed on the target genotypes here");
  if (t.gscale != 1 && !(allow_dosage && t.gscale == 2))
    cao.error("the reference PCA decomposed dosages on the 0..2 scale (gscale=", t.gscale,
              "), which is not the 0..1 coding used here. this combination was silently mis-scaled before; "
              "rerun the reference on the same genotype coding");
  // mirrors FileUSV::check_transform(): standardize_E() and
  // pcangsd_standardize_E() disagree on what this matrix is, so it has no inverse
  if (t.gscale == 2 && t.scale == SCALE_STANDARDIZE_GENETIC)
    cao.error(
        "the reference PCA standardised a dosage-scale matrix (gscale=2, scale=-9), which cannot be "
        "replayed on the target here. rerun the reference PCA so it decomposes centred dosages");
  if (t.ploidy != params.ploidy)
    cao.error("the reference PCA used ploidy ", t.ploidy, " but this run uses ", params.ploidy,
              "; pass --haploid consistently");
  const bool standardize = (t.scale == SCALE_STANDARDIZE_GENETIC);
  cao.print(tick.date(), "scaling the target genotypes as the reference did:",
            standardize ? "standardized" : (t.gscale == 2 ? "centred dosages (0..2)" : "centred only"),
            "(scale =", t.scale, ", ploidy =", t.ploidy, ", gscale =", t.gscale, ")");
  if (params.scale != t.scale)
    cao.warn("-C/--scale ", params.scale, " is ignored here; the scaling is taken from the reference PCA");
  return standardize;
}

// standardize_E() but driven by the reference's transform rather than this
// run's --scale, and using the reference's ploidy.
void Data::standardize_E_ref(const UsvTransform& t) {
  const double rploidy = sqrt((double)t.ploidy);
#pragma omp parallel for
  for (uint i = 0; i < nsnps; ++i) {
    const double f = F(i);
    const double sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) G.col(i) *= rploidy / sd;
  }
}

// the same, for one out-of-core block already read unstandardized
void Data::standardize_block_ref(const UsvTransform& t, uint64 start_idx, uint block_cols) {
  const double rploidy = sqrt((double)t.ploidy);
#pragma omp parallel for
  for (uint i = 0; i < block_cols; ++i) {
    const double f = F(start_idx + i);
    const double sd = sqrt(f * (1.0 - f));
    if (sd > VAR_TOL) G.col(i) *= rploidy / sd;
  }
}

void Data::pcangsd_standardize_E(const Mat2D& U, const Mat1D& svals, const Mat2D& VT) {
  if (params.scale != -9) return;
  cao.print(tick.date(), "begin to standardize the matrix for pcangsd procedure");
  const Mat2D US = U * svals.asDiagonal();
  Mat2D diag_threads = Mat2D::Zero(nsamples, omp_get_max_threads());  // a column per thread
  for_each_product_column(US, VT, 0, nsnps, [&](Eigen::Index jj, const auto& recon) {
    auto diag_private = diag_threads.col(omp_get_thread_num());
    {
      const uint j = (uint)jj;
      double p0, p1, p2, pt, pSum, tmp;
      const uint original = unpermuted_snp_index(j);
      const double f = F(original);
      double norm = sqrt(2.0 * f * (1.0 - f));
      uint s = params.filterSNP ? keepSNPs[original] : original;
      for (uint i = 0; i < nsamples; i++) {
        // Rescale individual allele frequencies
        pt = recon(i);
        pt = (pt + 2.0 * f) / 2.0;
        pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
        // Update e
        p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
        p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
        pSum = p0 + p1 + p2;
        G(i, j) = (p1 + 2 * p2) / pSum - 2.0 * f;
        if (norm > VAR_TOL) G(i, j) /= norm;

        // Update diag
        tmp = (0.0 - 2.0 * f) * (0.0 - 2.0 * f) * (p0 / pSum);
        tmp += (1.0 - 2.0 * f) * (1.0 - 2.0 * f) * (p1 / pSum);
        tmp += (2.0 - 2.0 * f) * (2.0 - 2.0 * f) * (p2 / pSum);
        if (norm > VAR_TOL) diag_private(i) += tmp / (2.0 * f * (1.0 - f));  // f = 0 gave inf
      }
    }
  });
  Dc = diag_threads.rowwise().sum();  // in thread order, not in the order a critical section was entered
}

void Data::predict_missing_E(const Mat2D& U, uint64 start_idx, uint64 stop_idx) {
  // assume full size G is in memory
  // Eigen::ColPivHouseholderQR<Mat2D> qr(U);  // nsamples x nk
  const uint nk = U.cols();
// get beta by linear regresion on data without NAs
// Ux = g, x is the beta
#pragma omp parallel for
  for (uint j = start_idx; j <= stop_idx; j++) {
    // find non-missing samples
    Int1D idx, na_idx;
    for (uint i = 0; i < nsamples; i++) {
      if (!C((uint64)j * nsamples + i))
        idx.push_back(i);
      else
        na_idx.push_back(i);
    }
    if (idx.size() == 0) continue;
    Mat1D v = U(idx, Eigen::all).colPivHouseholderQr().solve(G(idx, j));
    // predict for NAs in E
    for (auto i : na_idx) {  // samples need to be predicted
      G(i, j) = 0.0;
      for (uint k = 0; k < nk; ++k) {
        G(i, j) += U(i, k) * v(k);
      }
      G(i, j) = fmin(fmax(G(i, j), -F(j)), 1 - F(j));
    }
  }
}

void Data::estep_reg(const Mat2D& U, Eigen::Ref<Mat2D> Z, uint64 j0) {
  const double t0 = omp_get_wtime();
  struct AddTime {
    double t0;
    double& total;
    ~AddTime() { total += omp_get_wtime() - t0; }
  } add_time{t0, esteptime};
  const Eigen::Index n = Z.cols(), N = Z.rows(), k = U.cols();
  // the missing samples of each column, as a CSR: in-core from C, built once in
  // G's column order; out-of-core from the NaN of this block
  std::vector<uint64> lptr;
  std::vector<uint> lidx;
  auto build = [&](std::vector<uint64>& ptr, std::vector<uint>& idx, Eigen::Index ncol, auto is_missing) {
    ptr.assign(ncol + 1, 0);
#pragma omp parallel for
    for (Eigen::Index j = 0; j < ncol; ++j) {
      uint64 c = 0;
      for (Eigen::Index i = 0; i < N; ++i) c += is_missing(i, j);
      ptr[j + 1] = c;
    }
    for (Eigen::Index j = 0; j < ncol; ++j) ptr[j + 1] += ptr[j];
    idx.resize(ptr[ncol]);
#pragma omp parallel for
    for (Eigen::Index j = 0; j < ncol; ++j) {
      uint64 t = ptr[j];
      for (Eigen::Index i = 0; i < N; ++i)
        if (is_missing(i, j)) idx[t++] = (uint)i;
    }
  };
  const uint64* ptr;
  const uint* idx;
  if (!params.out_of_core) {
    if (miss_ptr.empty())
      build(miss_ptr, miss_idx, nsnps,
            [&](Eigen::Index i, Eigen::Index j) { return C[(uint64)unpermuted_snp_index(j) * nsamples + i]; });
    ptr = miss_ptr.data() + j0;
    idx = miss_idx.data();
  } else {
    build(lptr, lidx, n, [&](Eigen::Index i, Eigen::Index j) { return std::isnan(Z(i, j)); });
    ptr = lptr.data();
    idx = lidx.data();
  }
  // per thread, buffers of the largest size, sliced per site: a heap
  // allocation per site in this loop cost more than the arithmetic
#pragma omp parallel
  {
    Mat2D Urows(N, k), A(k, k), Aclip(k, k);  // rows [0, nm): the missing samples
    Mat1D beta(k), c(k), y(k), zm(N), znew(N);
    Eigen::LLT<Mat2D> llt(k);
#pragma omp for schedule(dynamic, 256)
    for (Eigen::Index j = 0; j < n; ++j) {
      const uint64 nm = ptr[j + 1] - ptr[j];
      if (nm == 0) continue;
      const uint* m = idx + ptr[j];
      auto z = Z.col(j);
      for (uint64 t = 0; t < nm; ++t) z(m[t]) = 0.0;
      if (k == 0) continue;
      // the products of this loop are k wide: as dots and axpys of columns,
      // which vectorize. A BLAS call per site (GEMV), or one GEMM for the
      // block, cost more, and Eigen's lazy products did not vectorize
      for (Eigen::Index kk = 0; kk < k; ++kk) beta(kk) = U.col(kk).dot(z);  // U'x, the observed calls only
      const double lo = -F(unpermuted_snp_index(j0 + j)), hi = 1.0 + lo;
      auto Um = Urows.topRows(nm);
      for (Eigen::Index kk = 0; kk < k; ++kk)
        for (uint64 t = 0; t < nm; ++t) Um(t, kk) = U(m[t], kk);
      // A = U_o'U_o = I - U_m'U_m. (From U_o when most calls are missing cost
      // more: a branch per sample. The cancellation loses ~1e-16 / rho of the
      // observed fraction rho, harmless for any site with calls.)
      A.setIdentity();  // the lower triangle is what LLT reads
      for (Eigen::Index a = 0; a < k; ++a)
        for (Eigen::Index b = 0; b <= a; ++b) A(a, b) -= Um.col(a).dot(Um.col(b));
      // fewer calls than PCs leave A singular: the tiny ridge picks the
      // solution with no component the observed calls cannot see
      A.diagonal().array() += 1e-10;
      // EMU's fixed point z_m = clip(U_m U'z) by an active set: the calls at a
      // bound count as observed there, the others are fitted by least squares
      // on all the rest, and a call moves to or off its bound until none does.
      // (EMU's own step contracts by the missing fraction per step.)
      auto v = zm.head(nm);
      auto state = znew.head(nm);  // 0 free, -1 at lo, 1 at hi
      state.setZero();
      Mat2D& Ac = Aclip;
      for (int it = 0; it < 10; ++it) {
        Ac = A;
        y = beta;
        for (uint64 t = 0; t < nm; ++t) {
          if (state(t) == 0) continue;
          for (Eigen::Index a = 0; a < k; ++a)
            for (Eigen::Index b = 0; b <= a; ++b) Ac(a, b) += Um(t, a) * Um(t, b);
          y.noalias() += (state(t) < 0 ? lo : hi) * Um.row(t).transpose();
        }
        c = llt.compute(Ac).solve(y);
        v = c(0) * Um.col(0);
        for (Eigen::Index kk = 1; kk < k; ++kk) v += c(kk) * Um.col(kk);
        bool changed = false;
        for (uint64 t = 0; t < nm; ++t) {
          const double st = state(t);
          if (st == 0 && v(t) < lo) state(t) = -1, changed = true;
          else if (st == 0 && v(t) > hi) state(t) = 1, changed = true;
          else if ((st < 0 && v(t) > lo) || (st > 0 && v(t) < hi)) state(t) = 0, changed = true;
        }
        if (!changed) break;
      }
      for (uint64 t = 0; t < nm; ++t) z(m[t]) = std::fmin(std::fmax(v(t), lo), hi);
    }
  }
}

void Data::estep_pcangsd(const Mat2D& U, Eigen::Ref<Mat2D> Z, uint64 j0) {
  const double t0 = omp_get_wtime();
  const Eigen::Index n = Z.cols(), N = Z.rows(), k = U.cols();
  if (k == 0) return;
  // inner steps per site and pass: on 1,000 simulated samples at low depth,
  // 1 step took 18 passes to converge, 3 steps 12, and 10 steps 12 too, at
  // the same EM time as 1 and about 1.7 times that of 3
  const int inner = 3;
#pragma omp parallel
  {
    Mat1D c(k), r(N);
#pragma omp for schedule(dynamic, 64)
    for (Eigen::Index j = 0; j < n; ++j) {
      const uint original = unpermuted_snp_index(j0 + j);
      const double f = F(original);
      const uint s = params.filterSNP ? keepSNPs[original] : original;
      auto e = Z.col(j);
      for (int it = 0; it < inner; ++it) {
        // r = U U'e as dots and axpys of the k columns, as in estep_reg()
        for (Eigen::Index kk = 0; kk < k; ++kk) c(kk) = U.col(kk).dot(e);
        r = c(0) * U.col(0);
        for (Eigen::Index kk = 1; kk < k; ++kk) r += c(kk) * U.col(kk);
        double delta = 0;
        for (Eigen::Index i = 0; i < N; ++i) {
          const double pt = fmin(fmax((r(i) + 2.0 * f) / 2.0, 1e-4), 1.0 - 1e-4);
          const double p0 = P(2 * i + 0, s) * (1.0 - pt) * (1.0 - pt);
          const double p1 = P(2 * i + 1, s) * 2 * pt * (1.0 - pt);
          const double p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * pt * pt;
          const double ei = (p1 + 2.0 * p2) / (p0 + p1 + p2) - 2.0 * f;
          delta = fmax(delta, fabs(ei - e(i)));
          e(i) = ei;
        }
        if (delta < 1e-6) break;
      }
    }
  }
  esteptime += omp_get_wtime() - t0;
}
