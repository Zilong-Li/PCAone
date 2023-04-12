#ifndef FILEBINARY_H_
#define FILEBINARY_H_

#include "Data.hpp"

class FileBin : public Data
{
  public:
    FileBin(const Param & params_) : Data(params_)
    {
        llog << timestamp() << "start parsing binary format" << std::endl;
        ifs_bin.open(params.binfile, std::ios::in | std::ios::binary);
        ifs_bin.read((char *)&nsamples, sizeof(uint64));
        ifs_bin.read((char *)&nsnps, sizeof(uint64));
        llog << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
        // ifs_bin.read((char*)&nimpute, sizeof(uint64));
        bytes_per_snp = nsamples * sizeof(double);
        // for (uint64 i = 0; i < nimpute; i++)
        //     ifs_bin.read((char*)&missing_points[i], sizeof(uint64));
    }

    ~FileBin() {}

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
    std::ifstream ifs_bin;
    uint64 bytes_per_snp;
    std::vector<double> inbed;
    // std::vector<uint64> missing_points;
    uint64 magic = sizeof(uint64) * 2;
};

#endif // FILEBINARY_H_
