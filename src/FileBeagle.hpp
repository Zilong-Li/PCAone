#ifndef PCAONE_FILEBEAGLE_
#define PCAONE_FILEBEAGLE_

#include "Data.hpp"
#include <zlib.h>

int tgets(gzFile gz, char** buf, uint64* l);

class FileBeagle : public Data
{
public:
    FileBeagle(const Param& params_) : Data(params_)
    {
        llog << timestamp() << "start parsing BEAGLE format" << std::endl;
        original = buffer = (char*)calloc(bufsize, sizeof(char));

        fp = gzopen(params.beagle.c_str(), "r");
        tgets(fp, &buffer, &bufsize);
        int nCol = 1;
        if (buffer != original)
            original = buffer;
        strtok_r(buffer, delims, &buffer);
        while (strtok_r(NULL, delims, &buffer))
            nCol++;
        if (nCol % 3)
        {
            throw std::runtime_error("number of columns should be a multiple of 3.\n");
        }
        nsamples = nCol / 3 - 1;
        // continue getting the number of sites
        // assume the number of columns of each line is the same. should check it first.
        buffer = original;
        nsnps = 0;
        while (tgets(fp, &buffer, &bufsize))
        {
            if (buffer != original)
                original = buffer;
            bgls.emplace_back((std::string)buffer);
            buffer = original;
            nsnps++;
        }
        gzclose(fp);

        if (params.nsnps > 0 && params.nsamples > 0 && (params.nsamples != nsamples || params.nsnps != nsnps))
        {
            llog << timestamp() << "nsamples and nsnps given by user does not match the file!" << std::endl;
            throw std::runtime_error("something wrong!\n");
        }
        P = MyMatrix::Zero(nsamples * 2, nsnps); // MyMatrix is column major
        llog << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
    }

    ~FileBeagle()
    {
    }

    virtual void read_all();
    // below are for blockwise, remain for future.
    virtual void check_file_offset_first_var()
    {
    }

    virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false)
    {
    }

    virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize)
    {
    }


private:
    gzFile fp = nullptr;
    char *original, *buffer;
    uint64 bufsize = (uint64)128 * 1024 * 1024;
    const char* delims = "\t \n";
    std::vector<std::string> bgls;
};



#endif // PCAONE_FILEBEAGLE_
