/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FileBgen.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "FileBgen.hpp"

#include <thread>

using namespace std;

void FileBgen::read_all() {
  uint i, j, k, gc;
  double gs, af;
  if (!params.pcangsd) {
    F = Mat1D::Zero(nsnps);
    G = Mat2D::Zero(nsamples, nsnps);
    if (params.impute) C = ArrBool::Zero(nsnps * nsamples);
    for (j = 0, k = 0; j < nsnps; j++) {
      try {
        auto var = bg->next_var();
        dosages = var.minor_allele_dosage();
        gc = 0;
        gs = 0.0;
// calculate allele frequency
#pragma omp parallel for reduction(+ : gc) reduction(+ : gs)
        for (i = 0; i < nsamples; i++) {
          if (!std::isnan(dosages[i])) {
            gs += dosages[i] / 2.0;  // map to [0, 1];
            gc += 1;
          }
        }
        if (gc == 0)
          af = 0.0;
        else
          af = (double)gs / gc;
        if (af > params.maf)
          F(k) = af;
        else
          continue;
// do centering and initialing
#pragma omp parallel for
        for (i = 0; i < nsamples; i++) {
          if (std::isnan(dosages[i])) {
            if (params.impute) C[k * nsamples + i] = 1;
            G(i, k) = 0;
          } else {
            if (params.impute) C[k * nsamples + i] = 0;
            G(i, j) = dosages[i] / 2.0 - F(k);  // map to [0, 1];
          }
        }
        k++;
      } catch (const std::out_of_range &e) {
        throw e.what();
      }
    }
    if (k == 0)
      cao.error("the number of SNPs after filtering should be 0!");
    else
      cao.print(tick.date(), "number of SNPs after filtering by MAF >", params.maf, ":", k);
    // resize G, F, C;
    nsnps = k;  // resize nsnps;
    G.conservativeResize(Eigen::NoChange, nsnps);
    F.conservativeResize(nsnps);
    C.conservativeResize(nsnps * nsamples);
  } else {
    // read all GP data into P;
    for (j = 0; j < nsnps; j++) {
      try {
        auto var = bg->next_var();
        probs1d = var.probs_1d();
#pragma omp parallel for
        for (i = 0; i < nsamples; i++) {
          P(i * 2 + 0, j) = probs1d[i * 3 + 0];
          P(i * 2 + 1, j) = probs1d[i * 3 + 1];
          // no need to parse probs1d[i * 3 + 2]
        }
      } catch (const std::out_of_range &e) {
        throw e.what();
      }
    }
    assert(j == nsnps);
    cao.print(tick.date(), "begin to estimate allele frequencies using GP");
    Mat1D Ft(nsnps);
    F = Mat1D::Constant(nsnps, 0.25);
    // run EM to estimate allele frequencies
    double diff;
    for (uint it = 0; it < params.maxiter; it++) {
#pragma omp parallel for
      for (uint j = 0; j < nsnps; j++) {
        Ft(j) = F(j);
        double p0, p1, p2, pt = 0.0;
        for (uint i = 0; i < nsamples; i++) {
          p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
          p1 = P(2 * i + 1, j) * 2 * F(j) * (1.0 - F(j));
          p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
          pt += (p1 + 2 * p2) / (2 * (p0 + p1 + p2));
        }
        F(j) = pt / (double)nsamples;
      }
      // calculate differences between iterations
      diff = sqrt((F - Ft).array().square().sum() / nsnps);
      // Check for convergence
      if (diff < params.tolmaf) {
        cao.print(tick.date(), "EM (MAF) converged at iteration:", it + 1);
        break;
      } else if (it == (params.maxiter - 1)) {
        cao.print(tick.date(), "EM (MAF) not converged");
      }
    }
    filter_snps_resize_F();
    // initial E which is G
    G = Mat2D::Zero(nsamples, nsnps);
#pragma omp parallel for
    for (j = 0; j < nsnps; j++) {
      double p0, p1, p2;
      uint s = params.keepsnp ? keepSNPs[j] : j;
      for (i = 0; i < nsamples; i++) {
        p0 = P(2 * i + 0, s) * (1.0 - F(j)) * (1.0 - F(j));
        p1 = P(2 * i + 1, s) * 2 * F(j) * (1.0 - F(j));
        p2 = (1 - P(2 * i + 0, s) - P(2 * i + 1, s)) * F(j) * F(j);
        G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
      }
    }
  }
}

void FileBgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) {
  read_bgen_block(G, F, bg, dosages, probs1d, frequency_was_estimated, nsamples, nsnps, blocksize, start_idx,
                  stop_idx, standardize);
}

void read_bgen_block(Mat2D &G, Mat1D &F, bgen::CppBgenReader *bg, float *dosages, float *probs1d,
                     bool &frequency_was_estimated, uint64 nsamples, uint64 nsnps, uint blocksize,
                     uint64 start_idx, uint64 stop_idx, bool standardize) {
  uint actual_block_size = stop_idx - start_idx + 1;
  uint i, j, snp_idx;
  if (G.cols() < blocksize || (actual_block_size < blocksize)) {
    G = Mat2D::Zero(nsamples, actual_block_size);
  }
  if (frequency_was_estimated) {
    for (i = 0; i < actual_block_size; ++i) {
      snp_idx = start_idx + i;
      auto var = bg->next_var();
      dosages = var.minor_allele_dosage();
#pragma omp parallel for
      for (j = 0; j < nsamples; j++) {
        if (std::isnan(dosages[j])) {
          G(j, i) = 0;
        } else {
          G(j, i) = dosages[j] / 2.0 - F(snp_idx);
        }
        if (standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
          G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
      }
    }
  } else {
    uint gc;
    double gs;
    for (i = 0; i < actual_block_size; ++i) {
      snp_idx = start_idx + i;
      auto var = bg->next_var();
      dosages = var.minor_allele_dosage();
      gc = 0;
      gs = 0.0;
#pragma omp parallel for reduction(+ : gc) reduction(+ : gs)
      for (j = 0; j < nsamples; j++) {
        if (std::isnan(dosages[j])) {
          G(j, i) = 0;
        } else {
          G(j, i) = dosages[j] / 2.0;
          gs += G(j, i);
          gc += 1;
        }
      }
      if (gc == 0) cao.error("the allele frequency should not be 0. do filtering first");
      F(snp_idx) = (double)gs / gc;
// do centering
#pragma omp parallel for
      for (j = 0; j < nsamples; j++) {
        if (std::isnan(dosages[j])) {
          G(j, i) = 0;
        } else {
          G(j, i) -= F(snp_idx);
        }
        if (standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
          G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
      }
    }
    if (stop_idx + 1 == nsnps) frequency_was_estimated = true;
  }
}

// this would be fast
int shuffle_bgen_to_bin(std::string &fin, std::string fout, uint gb, bool standardize) {
  cao.print(tick.date(), "begin to permute BGEN into BINARY file");
  bgen::CppBgenReader *bg = new bgen::CppBgenReader(fin, "", true);
  uint nsamples = bg->header.nsamples;
  uint nsnps = bg->header.nvariants;
  const uint ibyte = 4;
  uint64 bytes_per_snp = nsamples * ibyte;
  uint blocksize = 1073741824 * gb / bytes_per_snp;
  uint nblocks = (nsnps + blocksize - 1) / blocksize;
  std::ofstream ofs(fout + ".perm.bin", std::ios::binary);
  std::ofstream ofs2(fout + ".perm.txt");
  ofs.write((char *)&nsnps, ibyte);
  ofs.write((char *)&nsamples, ibyte);
  uint magic = ibyte * 2;
  float *dosages = nullptr;
  float *probs1d = nullptr;
  bool frequency_was_estimated = false;
  std::vector<uint> perm(nsnps);
  std::iota(perm.begin(), perm.end(), 0);
  auto rng = std::default_random_engine{};
  std::shuffle(perm.begin(), perm.end(), rng);
  Mat2D G;
  Mat1D F(nsnps);
  uint64 idx, cur = 0;
  for (uint i = 0; i < nblocks; i++) {
    auto start_idx = i * blocksize;
    auto stop_idx = start_idx + blocksize - 1;
    stop_idx = stop_idx >= nsnps ? nsnps - 1 : stop_idx;
    read_bgen_block(G, F, bg, dosages, probs1d, frequency_was_estimated, nsamples, nsnps, blocksize,
                    start_idx, stop_idx, standardize);
    for (Eigen::Index p = 0; p < G.cols(); p++, cur++) {
      ofs2 << perm[cur] << "\n";
      idx = magic + perm[cur] * bytes_per_snp;
      ofs.seekp(idx, std::ios_base::beg);
      ofs.write((char *)G.col(p).data(), bytes_per_snp);
    }
  }
  delete bg;
  fin = fout + ".perm.bin";
  return (nsnps == cur);
}

void permute_bgen_thread(uint nsamples, std::vector<int> idx, std::string fin, std::string fout,
                         int ithread) {
  fout = fout + ".perm." + to_string(ithread) + ".bgen";
  uint geno_len, compress_flag = 2, layout = 2, ploidy_n = 2;  // 1:zlib, 2:zstd
  float *probs = nullptr;
  bool phased = false;
  uint8_t bit_depth = 8;
  string metadata;
  vector<string> sampleids;
  bgen::CppBgenWriter bw(fout, nsamples, metadata, compress_flag, layout, sampleids);
  bgen::CppBgenReader br(fin, "", true);
  br.parse_all_variants();
  for (auto i : idx) {
    auto var = br.variants[i];
    probs = var.probs_1d();
    geno_len = nsamples * var.probs_per_sample();
    bw.write_variant_header(var.varid, var.rsid, var.chrom, var.pos, var.alleles, var.n_samples);
    bw.add_genotype_data(var.alleles.size(), probs, geno_len, ploidy_n, phased, bit_depth);
  }
}

PermMat permute_bgen(std::string &fin, std::string fout, int nthreads) {
  cao.print(tick.date(), "begin to permute BGEN file");
  uint nsamples, nsnps;
  {
    bgen::CppBgenReader br(fin, "", true);
    nsamples = br.header.nsamples;
    nsnps = br.header.nvariants;
  }
  vector<int> perm(nsnps);
  std::iota(perm.begin(), perm.end(), 0);
  auto rng = std::default_random_engine{};
  std::shuffle(perm.begin(), perm.end(), rng);
  vector<std::thread> threads;
  uint tn = (nsnps + nthreads - 1) / nthreads;  // evenly spread index
  for (int i = 0; i < nthreads; i++) {
    vector<int> idx(perm.begin() + tn * i, i == nthreads - 1 ? perm.end() : perm.begin() + tn * (i + 1));
    threads.emplace_back(permute_bgen_thread, nsamples, idx, fin, fout, i);
  }
  // Wait for all threads to finish execution
  for (auto &t : threads) t.join();
  // now cat all bgen files into big one
  uint compress_flag = 2, layout = 2;  // 1:zlib, 2:zstd
  string metadata;
  vector<string> sampleids;
  string out = fout + ".perm.bgen";
  bgen::CppBgenWriter bw(out, nsamples, metadata, compress_flag, layout, sampleids);
  std::ostreambuf_iterator<char> outIt(bw.handle);
  for (int i = 0; i < nthreads; i++) {
    fin = fout + ".perm." + to_string(i) + ".bgen";
    bgen::CppBgenReader br(fin, "", true);
    br.handle.seekg(br.offset);
    std::istreambuf_iterator<char> inIt(br.handle);
    std::istreambuf_iterator<char> endIt;
    std::copy(inIt, endIt, outIt);  // copy everything
    bw.n_variants += br.header.nvariants;
    std::remove(fin.c_str());  // now delete the temp file
  }
  fin = out;  // point to the new file
  PermMat P;
  P.indices() = Eigen::Map<Eigen::VectorXi>(perm.data(), perm.size());
  return P;
}
