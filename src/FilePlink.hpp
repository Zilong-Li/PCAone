#ifndef PCAONE_FILEPLINK_
#define PCAONE_FILEPLINK_

#include "Data.hpp"
#include "Utils.hpp"

class FileBed : public Data {
 public:
  //
  FileBed(const Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing PLINK format");
    std::string fbim = params.filein + ".bim";
    std::string ffam = params.filein + ".fam";
    nsamples = count_lines(ffam);
    nsnps = count_lines(fbim);
    cao.print(tick.date(), "N (# samples):", nsamples, ", M (# SNPs):", nsnps);
    snpmajor = true;
    bed_bytes_per_snp = (nsamples + 3) >> 2;
    std::string fbed = params.filein + ".bed";
    bed_ifstream.open(fbed, std::ios::in | std::ios::binary);
    if (!bed_ifstream.is_open()) cao.error("Cannot open bed file.");
    // check magic number of bed file
    uchar header[3];
    bed_ifstream.read(reinterpret_cast<char *>(&header[0]), 3);
    if ((header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
      cao.error("Incorrect magic number in plink bed file.");
    if (params.center) {
      centered_geno_lookup = Arr2D::Zero(4, nsnps);
      F = Mat1D::Zero(nsnps);
    }
    if (params.project > 0) {
      // read frq from original set
      cao.print(tick.date(), "read frequency of SNPs from the extended bim (.mbim)");
      F = read_frq(params.filebim);
      if (F.size() != nsnps) cao.error("the number of sites doesn't match each other");
    }
  }

  ~FileBed() override = default;

  void read_all() final;
  // for blockwise
  void check_file_offset_first_var() final;

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final;

 private:
  std::ifstream bed_ifstream;
  uint64 bed_bytes_per_snp;
  bool frequency_was_estimated = false;
  std::vector<uchar> inbed;
};

PermMat permute_plink(std::string &fin, const std::string &fout, uint gb, uint nbands);

#endif  // PCAONE_FILEPLINK_
