#ifndef FILEBINARY_H_
#define FILEBINARY_H_

#include "Data.hpp"

class FileBin : public Data {
 public:
  FileBin(Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing binary format");
    ifs_bin.open(params.filein, std::ios::in | std::ios::binary);
    is_zstd = isZstdCompressed(params.filein.c_str());
    if (is_zstd) {
      cao.error("does not support parsing zstd compressed binary fille yet");
    } else {
      ifs_bin.read((char *)&nsnps, ibyte);
      ifs_bin.read((char *)&nsamples, ibyte);
      cao.print(tick.date(), "shape of input matrix (features x samples) is", nsnps, "x", nsamples);
      bytes_per_snp = nsamples * ibyte;
    }
  }

  ~FileBin() {}

  virtual void read_all();
  // for blockwise
  virtual void check_file_offset_first_var();

  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize);

  virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D &U, const Mat1D &svals,
                                 const Mat2D &VT, bool standardize) {}

 private:
  std::ifstream ifs_bin;
  const uint ibyte = 4;
  uint64 bytes_per_snp;
  uint64 magic = ibyte * 2;
  bool is_zstd = false;
};

#endif  // FILEBINARY_H_
