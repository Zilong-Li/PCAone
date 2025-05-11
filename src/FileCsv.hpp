#ifndef PCAONE_FILECSV_
#define PCAONE_FILECSV_

#include "Data.hpp"
#include "zstdpp/zstdpp.hpp"

using Resources = zstdpp::stream::Resources;
using Context = zstdpp::stream::Context;

void parse_csvzstd(const std::string &fin, uint32_t &nsamples, uint32_t &nsnps, uint scale,
                   std::vector<double> &libsize, std::vector<size_t> &tidx, double &median_libsize);

void read_csvzstd_block(std::ifstream &ifs, Resources &res, Context &ctx, std::string &buffCur, int &isEmpty,
                        size_t &lastRet, uint blocksize, uint64_t start_idx, uint64_t stop_idx, Mat2D &G,
                        uint32_t nsamples, std::vector<double> &libsize, std::vector<size_t> &tidx,
                        double median_libsize, uint scale);

PermMat normCSV2BIN(std::string &, std::string, uint, uint, bool);

// for other types, assume data is already noralized only do centering
class FileCsv : public Data {
 public:
  FileCsv(const Param &params_) : Data(params_) {
    cao.print(tick.date(), "start parsing CSV format compressed by ZSTD");
    ifs.open(params.filein, std::ios::in | std::ios::binary);

    if (params.nsnps > 0 && params.nsamples > 0 && !params.cpmed) {
      if (params.verbose > 1) cao.print(tick.date(), "use numer of samples and features given by user.");
      nsamples = params.nsamples;
      nsnps = params.nsnps;
    } else {
      parse_csvzstd(params.filein, nsamples, nsnps, params.scale, libsize, tidx, median_libsize);
    }
    cao.print(tick.date(), "shape of input matrix (features x samples) is", nsnps, " x", nsamples);
  }

  ~FileCsv() override = default;

  void read_all() final;

  // for blockwise
  void check_file_offset_first_var() final;

  void read_block_initial(uint64, uint64, bool) final;

  void read_block_update(uint64, uint64, const Mat2D &, const Mat1D &, const Mat2D &, bool) final {}

 private:
  Resources res{};
  Context ctx{};
  std::ifstream ifs;
  std::string buffCur{""};
  int isEmpty{0};     // check if ifstream return nothing, which means EOF
  size_t lastRet{0};  // check the return value of decompress function
  std::vector<size_t> tidx;
  std::vector<double> libsize;
  double median_libsize;
};

#endif  // PCAONE_FILECSV_
