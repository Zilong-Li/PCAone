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
      // size_t read;
      // read = res.readFrom(ifs);
      // isEmpty = read == 0;
      // ZSTD_inBuffer input = {res.getRawInData(), read, 0};
      // if (input.pos < input.size) { // decompress once instead of loop
      //   ZSTD_outBuffer output = {res.getRawOutData(), res.getToWrite(), 0};
      //   lastRet = ctx(input, output);  // perform decompression
      //   if (ZSTD_isError(lastRet)) cao.error("Error: ZSTD decompression failed. lastRet:", lastRet);
      //   std::memcpy(&nsnps, (char *)res.getRawOutData(), sizeof(uint32_t));
      //   std::memcpy(&nsamples, (char*)res.getRawOutData() + sizeof(uint32_t), sizeof(uint32_t));
      // } else {
      //   cao.error("Error: can not allocate enough buffer!");
      // }
      // not eof, but end of processing, so reset
      lastRet = 0;
      isEmpty = 0;
    } else {
      ifs.read((char *)&nsnps, ibyte);
      ifs.read((char *)&nsamples, ibyte);
    }

    if (params.nsnps > 0 && params.nsamples > 0) {
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
  zstdpp::stream::Resources res{};
  zstdpp::stream::Context ctx{};
  std::ifstream ifs;
  int isEmpty{0};     // check if ifstream return nothing, which means EOF
  size_t lastRet{0};  // check the return value of decompress function
  
  const uint ibyte = 4;
  uint64 bytes_per_snp;
  bool is_zstd = false;
};

#endif  // FILEBINARY_H_
