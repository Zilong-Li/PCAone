#ifndef PCAONE_FILECSV_
#define PCAONE_FILECSV_

#include "Data.hpp"

void parse_csvzstd(ZstdDS& zbuf, uint& nsamples, uint& nsnps, uint scale, std::vector<double>& libsize,
                   std::vector<size_t>& tidx, double& median_libsize);

void read_csvzstd_block(ZstdDS& zbuf, int blocksize, uint64 start_idx, uint64 stop_idx, Mat2D& G,
                        uint nsamples, std::vector<double>& libsize, std::vector<size_t>& tidx,
                        double median_libsize, uint scale);

PermMat shuffle_csvzstd_to_bin(std::string& fin, std::string fout, uint gb, uint scale);

// for other types, assume data is already noralized only do centering
class FileCsv : public Data {
 public:
  FileCsv(const Param& params_) : Data(params_) {
    cao.print(tick.date(), "start parsing CSV format compressed by ZSTD");

    if (params.nsnps > 0 && params.nsamples > 0 && !params.cpmed) {
      cao.print(tick.date(), "use nsamples and nsnps given by user.");
      nsamples = params.nsamples;
      nsnps = params.nsnps;
    } else {
      zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
      parse_csvzstd(zbuf, nsamples, nsnps, params.scale, libsize, tidx, median_libsize);
    }
    cao.print(tick.date(), "shape of input matrix (features x samples) is", nsnps, " x", nsamples);
  }

  ~FileCsv() {}

  virtual void read_all();
  // for blockwise
  virtual void check_file_offset_first_var();

  virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);

  virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const Mat2D& U, const Mat1D& svals,
                                 const Mat2D& VT, bool standardize = false) {}

 private:
  ZstdDS zbuf;
  std::vector<size_t> tidx;
  std::vector<double> libsize;
  double median_libsize;
};

#endif  // PCAONE_FILECSV_
