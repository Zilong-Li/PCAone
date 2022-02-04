#ifndef _FILECSV_H_
#define _FILECSV_H_

#include "Data.hpp"
#include "zstd.h"

// assume data is already noralized
// only do centering
class FileCsv : public Data
{
public:
    FileCsv(const Param &params_) : Data(params_)
        {
            std::cout << timestamp() << "start parsing CSV format" << std::endl;
            buffInTmp.reserve(buffInSize);
            buffOutTmp.reserve(buffOutSize);
            buffCur = "";

            if (params.nsnps > 0 && params.nsamples > 0) {
                std::cout << timestamp() << "use nsamples and nsnps given by user." << std::endl;
                nsamples = params.nsamples;
                nsnps = params.nsnps;

            } else {
                auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));
                auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));

                size_t read, i, p, ncol = 0, lastCol = 0;
                int isEmpty = 1;
                nsnps = 0;
                fin = fopenOrDie(params.csvfile.c_str(), "rb");
                while ((read = freadOrDie(buffIn, buffInSize, fin))) {
                    isEmpty = 0;
                    ZSTD_inBuffer input = {buffIn, read, 0};
                    while (input.pos < input.size) {
                        ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
                        lastRet = ZSTD_decompressStream(dctx, &output, &input);
                        buffCur += std::string((char *)buffOut, output.pos);
                        while ( (p = buffCur.find("\n")) != std::string::npos ) {
                            nsnps++;
                            buffLine = buffCur.substr(0, p);
                            buffCur.erase(0, p+1);
                            lastCol = ncol;
                            ncol = 1;
                            for (i = 0; i < buffLine.size(); i++) {
                                if (buffLine[i] == ',') ncol++;
                            }
                            if (nsnps > 2 && (lastCol != ncol)) {
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

                lastRet = 1;
                nsamples = ncol;
                tidx.resize(nsamples + 1);  // tidx[0] = 0;
            }

            std::cout << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
        }


    ~FileCsv() {
        ZSTD_freeDCtx(dctx);
        fcloseOrDie(fin);
    }

    virtual void read_all_and_centering();
    // for blockwise
    virtual void check_file_offset_first_var();

    virtual void read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);

    virtual void read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize = false) {}

private:
    FILE *fin = nullptr;
    size_t const buffInSize = ZSTD_DStreamInSize();
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    ZSTD_DCtx *const dctx = ZSTD_createDCtx();
    std::string buffCur, buffLine, buffInTmp, buffOutTmp;
    size_t lastRet = 1;
    std::vector<size_t> tidx;
};

#endif // FILECSV_H_
