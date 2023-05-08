#ifndef PCAONE_FILEBGEN_
#define PCAONE_FILEBGEN_

#include "Data.hpp"
#include "bgen/reader.h"
#include "bgen/writer.h"

// const double GENOTYPE_THRESHOLD = 0.9;
// const double BGEN_MISSING_VALUE = -9;
// const double BGEN2GENO[4] = {0, 0.5, 1, BGEN_MISSING_VALUE};

void read_bgen_block(MyMatrix & G,
                     MyVector & F,
                     bgen::CppBgenReader * bg,
                     float * dosages,
                     float * probs1d,
                     bool & frequency_was_estimated,
                     uint64 nsamples,
                     uint64 nsnps,
                     uint blocksize,
                     uint64 start_idx,
                     uint64 stop_idx,
                     bool standardize);

int shuffle_bgen_to_bin(std::string bgenfile, std::string binfile, uint gb, bool standardize);

void permute_bgen(std::string & fin, std::string fout);

class FileBgen : public Data
{
  public:
    // using Data::Data;
    FileBgen(Param & params_) : Data(params_)
    {
        llog << timestamp() << "start parsing BGEN format" << std::endl;
        bg = new bgen::CppBgenReader(params.filein, "", true);
        nsamples = bg->header.nsamples;
        nsnps = bg->header.nvariants;
        llog << timestamp() << "the layout of bgen file is " << bg->header.layout << ". N samples is "
             << nsamples << ". M snps is " << nsnps << std::endl;
    }

    ~FileBgen()
    {
        delete bg;
    }

    virtual void read_all();
    // for blockwise
    virtual void check_file_offset_first_var()
    {
        bg->offset = bg->header.offset + 4;
    }

    virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);

    virtual void read_block_update(uint64 start_idx,
                                   uint64 stop_idx,
                                   const MyMatrix & U,
                                   const MyVector & svals,
                                   const MyMatrix & VT,
                                   bool standardize)
    {
    }

  private:
    bgen::CppBgenReader * bg;
    float * dosages = nullptr;
    float * probs1d = nullptr;
    bool frequency_was_estimated = false;
};

#endif // PCAONE_FILEBGEN_
