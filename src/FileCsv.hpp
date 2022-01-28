#ifndef _FILECSV_H_
#define _FILECSV_H_

#include "Data.hpp"
#include "zstd.h"

// assume data is already noralized
class FileCsv : public Data
{
public:
    FileCsv(const Param &params_) : Data(params_)
    {
        size_t const buffInSize = ZSTD_DStreamInSize();
        buffInTmp.reserve(buffInSize);
        auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));

        auto buffOutSize = ZSTD_DStreamOutSize();
        buffOutTmp.reserve(buffOutSize);
        auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));

        ZSTD_DCtx *const dctx = ZSTD_createDCtx();

        size_t const toRead = buffInSize;
        size_t lastRet = 0;
        size_t read, p, ncol, lastCol;
        int isEmpty = 1;
        fin = fopenOrDie(params.csvfile.c_str(), "rb");
        nsamples = 0;
        buffCur = "";
        while ((read = freadOrDie(buffIn, toRead, fin))) {
            isEmpty = 0;
            ZSTD_inBuffer input = {buffIn, read, 0};
            while (input.pos < input.size) {
                ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
                size_t const ret = ZSTD_decompressStream(dctx, &output, &input);
                lastRet = ret;
                buffCur += std::string((char *)buffOut, output.pos);
                while (( p = buffCur.find("\n") ) != std::string::npos) {
                    nsamples++;
                    buffLine = buffCur.substr(0, p);
                    buffCur.erase(0, p + 1);
                    lastCol = ncol;
                    ncol = 1;
                    while (( p = buffLine.find(",") ) != std::string::npos) {
                        ncol++;
                        buffLine.erase(0, p + 1);
                    }
                    if (nsamples > 2 && (lastCol != ncol)) {
                        throw std::invalid_argument("the csv file has unaligned columns\n");
                    }
                }
            }
        }

        if (isEmpty) {
            throw std::invalid_argument("input file is empty.\n");
        }

        if (lastRet != 0) {
            throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");
        }

        nsnps = ncol;
        ZSTD_freeDCtx(dctx);
        fcloseOrDie(fin);

        std::cout << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
    }

    ~FileCsv() {}

    virtual void read_all_and_centering();
    // for blockwise
    virtual void read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {}
    virtual void read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize = false) {}
    virtual void check_file_offset_first_var() {}

private:
    FILE *fin;
    std::string buffCur, buffLine, buffInTmp, buffOutTmp;
};

#endif // FILECSV_H_
