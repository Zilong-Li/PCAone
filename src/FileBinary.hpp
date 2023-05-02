#ifndef FILEBINARY_H_
#define FILEBINARY_H_

#include "Data.hpp"
#define ZSTD_STATIC_LINKING_ONLY
#include "zstd.h"

bool isZstdCompressed(const char * filename);

class FileBin : public Data
{
  public:
    FileBin(const Param & params_) : Data(params_)
    {
        llog << timestamp() << "start parsing binary format" << std::endl;
        ifs_bin.open(params.filein, std::ios::in | std::ios::binary);
        is_zstd = isZstdCompressed(params.filein.c_str());
        if(is_zstd)
        {
            inbuf.resize(buffInSize);
            outbuf.resize(buffOutSize);
            ifs_bin.read(inbuf.data(), buffInSize);
            bytesRead = ifs_bin.gcount();
            ZSTD_inBuffer input = {inbuf.data(), bytesRead, 0};
            while(input.pos < input.size)
            {

                ZSTD_outBuffer output = {outbuf.data(), buffOutSize, 0};
                const size_t ret = ZSTD_decompressStream(dctx, &output, &input);
                if(ZSTD_isError(ret))
                {
                    std::cerr << "Decompression error: " << ZSTD_getErrorName(ret) << std::endl;
                    exit(1);
                }
                if(!is_magic_read)
                {
                    uint * magic_ = reinterpret_cast<uint *>(outbuf.data());
                    nsamples = magic_[0];
                    nsnps = magic_[1];
                    is_magic_read = true;
                    dat = reinterpret_cast<float *>(outbuf.data() + sizeof(uint) * 2);
                }
                else
                {
                    dat = reinterpret_cast<float *>(outbuf.data());
                }
            }
        }
        else
        {
            ifs_bin.read((char *)&nsamples, ibyte);
            ifs_bin.read((char *)&nsnps, ibyte);
            llog << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
            bytes_per_snp = nsamples * ibyte;
        }
    }

    ~FileBin()
    {
        ZSTD_freeDCtx(dctx);
    }

    virtual void read_all();
    // for blockwise
    virtual void check_file_offset_first_var();

    virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize);

    virtual void read_block_update(uint64 start_idx,
                                   uint64 stop_idx,
                                   const MyMatrix & U,
                                   const MyVector & svals,
                                   const MyMatrix & VT,
                                   bool standardize)
    {
    }

  private:
    ZSTD_DCtx * const dctx = ZSTD_createDCtx();
    size_t const buffInSize = ZSTD_DStreamInSize();
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    size_t bytesRead;
    std::ifstream ifs_bin;
    std::vector<char> inbuf;
    std::vector<char> outbuf;
    const uint ibyte = 4;
    uint64 bytes_per_snp;
    uint64 magic = ibyte * 2;
    // uint in, im;
    bool is_zstd = false, is_magic_read = false;
    float * dat = nullptr;
};

#endif // FILEBINARY_H_
