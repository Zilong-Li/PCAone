#ifndef PCAONE_FILEPLINK_
#define PCAONE_FILEPLINK_

#include <memory>

#include "Data.hpp"
#include "Prefetch.hpp"
#include "Utils.hpp"

class FileBed : public Data {
 public:
  //
  FileBed(const Param& params_)
      : Data(params_) {
    cao.print(tick.date(), "start parsing PLINK format");
    std::string fbim = params.filein + ".bim";
    std::string ffam = params.filein + ".fam";
    nsamples = count_lines(ffam);
    nsnps = count_lines(fbim);
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps);
    if (!nsamples || !nsnps) cao.error("BED input must contain at least one sample and one SNP.");
    snpmajor = true;
    bed_bytes_per_snp = (nsamples + 3) >> 2;
    std::string fbed = params.filein + ".bed";
    bed_path = fbed;
    bed_ifstream.open(fbed, std::ios::in | std::ios::binary);
    if (!bed_ifstream.is_open()) cao.error("Cannot open bed file.");
    // check magic number of bed file
    uchar header[3];
    bed_ifstream.read(reinterpret_cast<char*>(&header[0]), 3);
    if (!bed_ifstream || (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
      cao.error("Incorrect magic number in plink bed file.");
    if (params.center) centered_geno_lookup = Arr2D::Zero(4, nsnps);
    if (params.dopca) F = Mat1D::Zero(nsnps);  // initial F
  }

  ~FileBed() override;

  void read_all() final;
  // Called after prepare(), before any genotype reads or frequency estimates.
  // Shuffles the SNPs for out-of-core winSVD. By default the blocks read the
  // shuffled SNPs from the input; with --bed-copy they are written to
  // <out>.perm.*, which replaces the input. Returns true if it wrote the copy.
  bool apply_permutation(Param& config);
  // for blockwise
  void check_file_offset_first_var() final;

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D&, const Mat1D&, const Mat2D&, bool) final;

 private:
  std::ifstream bed_ifstream;
  uint64 bed_bytes_per_snp;
  bool frequency_was_estimated = false;
  uint64 nmono_seen = 0;  // sites with MAF=0 met while estimating F block by block
  std::vector<uchar> inbed;
  std::string bed_path;  // the BED read by the blocks: the input, or its permuted copy
  // the source SNP of each logical SNP, when the blocks read them from the input
  std::vector<uint32_t> bed_order;
  // out-of-core: reads the next block while the current one is used (--no-prefetch: off)
  std::unique_ptr<PCAone::RecordPrefetcher> prefetcher;
  // the packed records of SNPs [start_idx, start_idx + count), valid until the next call
  const uchar* read_records(uint64 start_idx, uint64 count);
};

// Shuffles the SNPs into random buckets of `bucket` SNPs, source order kept
// within each, and writes <fout>.perm.{bed,bim,fam}. fin becomes <fout>.perm.
PermMat permute_plink(std::string& fin, const std::string& fout, uint gb, uint64 bucket, int seed);

// The order permute_plink() writes, without writing it: the SNPs shuffled at
// random (--seed) into buckets of `bucket` SNPs, each in source order.
std::vector<uint32_t> bed_bucket_order(uint64 nsnps, uint64 bucket, int seed);

// 1 if the file is on a rotating disk, 0 if not, -1 if unknown (not Linux, or
// a network or virtual file system)
int on_rotating_disk(const std::string& path);

#endif  // PCAONE_FILEPLINK_
