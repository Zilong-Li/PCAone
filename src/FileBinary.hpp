#ifndef FILEBINARY_H_
#define FILEBINARY_H_

#include "Data.hpp"

class FileBin : public Data {
 public:
  FileBin(Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing binary format");
    is_zstd = isZstdCompressed(params.filein.c_str());
    ifs_bin.open(params.filein, std::ios::in | std::ios::binary);
    if (is_zstd) {
      // cao.error("does not support parsing zstd compressed binary fille yet");
      char inBuf[16 * 1024] = {0};  // should be enough for zstd magic + nsamples + nsnps
      ifs_bin.read(inBuf, 16 * 1024);
      size_t bytesRead = ifs_bin.gcount();
      if (bytesRead == 0) cao.error("file is empty!");
      ZSTD_DCtx *dctx = ZSTD_createDCtx();
      ZSTD_inBuffer input = {inBuf, bytesRead, 0};
      // Buffer for decompressed output - only need 8 bytes for two uint32_t
      rawBuf.resize(8);
      ZSTD_outBuffer output = {rawBuf.data(), rawBuf.size(), 0};
      size_t ret = ZSTD_decompressStream(dctx, &output, &input);
      if (ZSTD_isError(ret)) cao.error("Error: ZSTD decompression failed");
      //  Check if we got enough decompressed data
      if (output.pos >= 8) {
        std::memcpy(&nsnps, rawBuf.data(), sizeof(uint32_t));
        std::memcpy(&nsamples, rawBuf.data() + sizeof(uint32_t), sizeof(uint32_t));
      } else {
        ZSTD_freeDCtx(dctx);
        cao.error("// FIXME: !");
      }
      ZSTD_freeDCtx(dctx);
      ifs_bin.close();
      zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
    } else {
      ifs_bin.read((char *)&nsnps, ibyte);
      ifs_bin.read((char *)&nsamples, ibyte);
    }

    cao.print(tick.date(), "shape of input matrix (features x samples) is", nsnps, "x", nsamples);
    bytes_per_snp = nsamples * ibyte;
  }

  ~FileBin() override = default;

  void read_all() final;

  // for blockwise
  void check_file_offset_first_var() final;

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final {}

 private:
  ZstdDS zbuf;
  std::vector<char> rawBuf;      // decompressed bytes, will keep appending chars

  std::ifstream ifs_bin;
  const uint ibyte = 4;
  uint64 bytes_per_snp;
  bool is_zstd = false;
};

#endif  // FILEBINARY_H_
