#ifndef PCAONE_FILEBGEN_
#define PCAONE_FILEBGEN_

#include "bgen/reader.h"
#include "BgenBlock.hpp"
#include "Data.hpp"

#include <omp.h>

#include <chrono>
#include <memory>
#include "Utils.hpp"

// const double GENOTYPE_THRESHOLD = 0.9;
// const double BGEN_MISSING_VALUE = -9;
// const double BGEN2GENO[4] = {0, 0.5, 1, BGEN_MISSING_VALUE};

class FileBgen : public Data {
 public:
  // using Data::Data;
  FileBgen(const Param& params_)
      : Data(params_) {
    cao.print(tick.date(), "start parsing BGEN format");
    bg = new bgen::CppBgenReader(params.filein, "", true);
    nsamples = bg->header.nsamples;
    nsnps = bg->header.nvariants;
    if (params.dopca) F = Mat1D::Zero(nsnps);  // initial F
    cao.print(tick.date(), "N(#samples) =", nsamples, ", M(#SNPs) =", nsnps);
    cao.print(tick.date(), "the layout is", bg->header.layout, ", compressed by",
              bg->header.compression == 2 ? "zstd" : (bg->header.compression == 1 ? "zlib" : "none"));
    if (!params.pcangsd) {
      reader_threads = std::max(1, omp_get_max_threads());
      const auto t0 = std::chrono::steady_clock::now();
      try {
        reader = std::make_unique<PCAone::BgenBlockReader>(params.filein, bg->header.layout, bg->header.compression,
                                                           nsamples, bg->header.offset + 4, nsnps, reader_threads);
      } catch (const std::exception& e) {
        cao.error(e.what());
      }
      thread_dosages.resize(reader_threads, std::vector<float>(nsamples));
      cao.print(tick.date(), "indexed", reader->bytes(0, nsnps), "bytes of BGEN variants in",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), "seconds");
    }
  }

  ~FileBgen() override {
    if (reader) reader->cancel();
    delete bg;
  }

  void read_all() final;

  // for blockwise
  void check_file_offset_first_var() final { bg->offset = bg->header.offset + 4; }

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D&, const Mat1D&, const Mat2D&, bool) final {}

 private:
  bgen::CppBgenReader* bg;
  std::vector<float> probs1d;
  bool frequency_was_estimated = false;
  // the dosages of all variants but --pcangsd's probabilities: each thread
  // reads and decodes variants of its own, and the out-of-core blocks request
  // the next block in the background (--no-prefetch: off)
  std::unique_ptr<PCAone::BgenBlockReader> reader;
  int reader_threads = 1;
  std::vector<std::vector<float>> thread_dosages;
  // the variants of a block in file order, and the column of each
  struct ReadRequest {
    uint32_t variant, column;
  };
  std::vector<ReadRequest> requests;
  std::vector<uint32_t> next_variants;
  uint64 last_end = 0, first_count = 0;
  void begin_block(uint64 start_idx, uint64 stop_idx);
  // the file indices of the logical SNPs [first, first + count), sorted
  void variants_of(uint64 first, uint64 count, std::vector<uint32_t>& out) const;
};

// The out-of-core winSVD permutation of the BGEN variants: a shuffle by --seed.
// The variants are read from the input in this order; nothing is rewritten.
PermMat compute_bgen_perm(uint nsnps, int seed);

#endif  // PCAONE_FILEBGEN_
