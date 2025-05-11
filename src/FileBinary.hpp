#ifndef FILEBINARY_H_
#define FILEBINARY_H_

#include "Data.hpp"
#include "zstdpp/zstdpp.hpp"

using Resources = zstdpp::stream::Resources;
using Context = zstdpp::stream::Context;

class FileBin : public Data {
 public:
  FileBin(Param &params_) : Data(params_) {
    is_zstd = isZstdCompressed(params.filein.c_str());
    cao.print(tick.date(), "start parsing binary format with zstd =", is_zstd);
    ifs.open(params.filein, std::ios::binary);
    if (is_zstd) {
      char inBuf[32 * 1024] = {0};  // should be enough for zstd magic + nsamples + nsnps
      ifs.read(inBuf, 32 * 1024);
      size_t bytesRead = ifs.gcount();
      if (bytesRead == 0) cao.error("file is empty!");
      ZSTD_DCtx *dctx = ZSTD_createDCtx();
      ZSTD_inBuffer input = {inBuf, bytesRead, 0};
      std::vector<char> rawBuf(8);
      ZSTD_outBuffer output = {rawBuf.data(), rawBuf.size(), 0};
      size_t ret = ZSTD_decompressStream(dctx, &output, &input);
      if (ZSTD_isError(ret)) cao.error("Error: ZSTD decompression failed");
      //  check if we got enough decompressed data
      if (output.pos >= 8) {
        std::memcpy(&nsnps, rawBuf.data(), sizeof(uint32_t));
        std::memcpy(&nsamples, rawBuf.data() + sizeof(uint32_t), sizeof(uint32_t));
      } else {
        ZSTD_freeDCtx(dctx);
        cao.error("// FIXME: edge cases for zstd !");
      }
      ZSTD_freeDCtx(dctx);

    } else {
      ifs.read((char *)&nsnps, ibyte);
      ifs.read((char *)&nsamples, ibyte);
    }

    if (params.nsnps > 0 && params.nsamples > 0) {
      cao.warn(tick.date(), "specify the numer of samples and features in your command is not needed!");
      nsamples = params.nsamples;
      nsnps = params.nsnps;
    }

    cao.print(tick.date(), "shape of input matrix (features x samples) is", nsnps, " x", nsamples);
    bytes_per_snp = nsamples * ibyte;
  }

  ~FileBin() override = default;

  void read_all() final;

  // for blockwise
  void check_file_offset_first_var() final;

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final {}

 private:
  Resources res{};
  Context ctx{};
  std::ifstream ifs;
  int isEmpty{0};                       // check if ifstream return nothing, which means EOF
  size_t lastRet{0};                    // check the return value of decompress function
  std::vector<char> decompressedBytes;  // store decompressed bytes

  const uint ibyte = 4;
  uint64 bytes_per_snp;
  bool is_zstd = false;
};

#endif  // FILEBINARY_H_
