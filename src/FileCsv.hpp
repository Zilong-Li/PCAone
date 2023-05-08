#ifndef PCAONE_FILECSV_
#define PCAONE_FILECSV_

#include "Data.hpp"
#include "zstd.h"

struct ZstdBuffer
{
    ZstdBuffer()
    {
        buffInTmp.reserve(buffInSize);
        buffOutTmp.reserve(buffOutSize);
    }
    ~ZstdBuffer()
    {
        ZSTD_freeDCtx(dctx);
        fcloseOrDie(fin);
    }
    FILE * fin = nullptr;
    size_t const buffInSize = ZSTD_DStreamInSize();
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    ZSTD_DCtx * const dctx = ZSTD_createDCtx();
    size_t lastRet = 1;
    std::string buffCur = "";
    std::string buffLine, buffInTmp, buffOutTmp;
};

void parse_csvzstd(ZstdBuffer & zbuf,
                   uint & nsamples,
                   uint & nsnps,
                   uint scale,
                   std::vector<double> & libsize,
                   std::vector<size_t> & tidx,
                   double & median_libsize);

void read_csvzstd_block(ZstdBuffer & zbuf,
                        int blocksize,
                        uint64 start_idx,
                        uint64 stop_idx,
                        MyMatrix & G,
                        uint nsamples,
                        std::vector<double> & libsize,
                        std::vector<size_t> & tidx,
                        double median_libsize,
                        uint scale);

int shuffle_csvzstd_to_bin(std::string & fin, std::string fout, uint gb, uint scale);

// assume data is already noralized
// only do centering
class FileCsv : public Data
{
  public:
    FileCsv(Param & params_) : Data(params_)
    {
        llog << timestamp() << "start parsing CSV format compressed by ZSTD" << std::endl;
        if(params.nsnps > 0 && params.nsamples > 0 && !params.cpmed)
        {
            llog << timestamp() << "use nsamples and nsnps given by user." << std::endl;
            nsamples = params.nsamples;
            nsnps = params.nsnps;
        }
        else
        {
            zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
            parse_csvzstd(zbuf, nsamples, nsnps, params.scale, libsize, tidx, median_libsize);
        }
        llog << timestamp() << "shape of input matrix(samples x features) is (" << nsamples << ", " << nsnps
             << ")" << std::endl;
    }

    ~FileCsv() {}

    virtual void read_all();
    // for blockwise
    virtual void check_file_offset_first_var();

    virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);

    virtual void read_block_update(uint64 start_idx,
                                   uint64 stop_idx,
                                   const MyMatrix & U,
                                   const MyVector & svals,
                                   const MyMatrix & VT,
                                   bool standardize = false)
    {
    }

  private:
    ZstdBuffer zbuf;
    std::vector<size_t> tidx;
    std::vector<double> libsize;
    double median_libsize;
};

#endif // PCAONE_FILECSV_
