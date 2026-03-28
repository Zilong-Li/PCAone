/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/FilePgen.hpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef PCAONE_FILEPGEN_
#define PCAONE_FILEPGEN_

#include "Data.hpp"
#include "Utils.hpp"
#include "pgenlib/pgenlibr.h"

class FilePgen : public Data {
 public:
  FilePgen(const Param& params_) : Data(params_) {
    cao.print(tick.date(), "start parsing PLINK2 PGEN format");
    std::string fpsam = params.filein + ".psam";
    std::string fpgen = params.filein + ".pgen";
    // Count samples from .psam (data lines do not start with '#')
    {
      std::ifstream fin(fpsam);
      if (!fin.is_open()) cao.error("Cannot open psam file: " + fpsam);
      std::string line;
      while (std::getline(fin, line))
        if (!line.empty() && line[0] != '#') ++nsamples;
    }
    reader_threads = std::max(1u, params.threads);
    reader.Load(fpgen, nsamples, {}, reader_threads);
    nsnps = reader.GetVariantCt();
    dosage_mode = (!params.hardcall) && reader.DosagePresent();
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps, ". dosage_mode:", dosage_mode);
    
    snpmajor = true;
    thread_bufs.resize(reader_threads, std::vector<double>(nsamples));
    if (params.center) centered_geno_lookup = Arr2D::Zero(4, nsnps);
    if (params.dopca) F = Mat1D::Zero(nsnps);  // initial F

  }

  ~FilePgen() override = default;

  void read_all() final;
  // PgenReader supports random access; no file-offset state to reset.
  void check_file_offset_first_var() final {}
  void read_block_initial(uint64, uint64, bool) final;
  void read_block_update(uint64, uint64, const Mat2D&, const Mat1D&, const Mat2D&, bool) final;

 private:
  PgenReader reader;
  uint reader_threads = 1;
  std::vector<std::vector<double>> thread_bufs;
  bool frequency_was_estimated = false;
  bool dosage_mode = false;
  // ReadHardcalls returns 0.0/1.0/2.0/-3.0; map to lookup index 0/1/2/3
  static int pgen_code(double v) { return (v == -3.0) ? 3 : static_cast<int>(v); }
};

/// Compute a band-interleaving permutation for PGEN logical (no-I/O) permutation.
/// Equivalent to the index math in permute_plink but without any file writes.
PermMat compute_pgen_perm(uint nsnps, uint nbands);

#endif  // PCAONE_FILEPGEN_
