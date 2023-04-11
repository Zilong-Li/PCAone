#ifndef PCAONE_FILEPLINK_
#define PCAONE_FILEPLINK_

#include "Data.hpp"

/**
 * Recode genotype codes to allelic dosages of first allele in .bim file (A1),
 * similarly to .raw files generated with '--recode A' in PLINK. A coding for
 * the missing value needs to be provided in 'na_value'.
 * 00 ->  2 (copies of A1)
 * 10 ->  1 (copy of A1)
 * 11 ->  0 (copy of A1)
 * 01 ->  3 (missing)
 */
const double BED_MISSING_VALUE = -9;
const double BED2GENO[4] = {1, BED_MISSING_VALUE, 0.5, 0};

class FileBed : public Data
{
public:
    //
    FileBed(const Param& params_) : Data(params_)
    {
        llog << timestamp() << "start parsing PLINK format" << std::endl;
        std::string fbim = params.bed_prefix + ".bim";
        std::string ffam = params.bed_prefix + ".fam";
        nsamples = count_lines(ffam);
        nsnps = count_lines(fbim);
        snpmajor = true;
        llog << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
        bed_bytes_per_snp = (nsamples + 3) >> 2; // get ceiling(nsamples/4)
        std::string fbed = params.bed_prefix + ".bed";
        bed_ifstream.open(fbed, std::ios::in | std::ios::binary);
        if (!bed_ifstream.is_open())
        {
            throw std::invalid_argument(colerror + "Cannot open bed file.\n");
        }
        // check magic number of bed file
        uchar header[3];
        bed_ifstream.read(reinterpret_cast<char*>(&header[0]), 3);
        if ((header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
        {
            throw std::invalid_argument(colerror + "Incorrect magic number in plink bed file.\n");
        }
    }

    ~FileBed()
    {
    }

    virtual void read_all();
    // for blockwise
    virtual void check_file_offset_first_var();

    virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize);

    virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize);

private:
    std::ifstream bed_ifstream;
    uint64 bed_bytes_per_snp;
    bool frequency_was_estimated = false;
    std::vector<uchar> inbed;
};

void permute_plink(std::string& fin, const std::string& fout, uint gb, uint nbands);

#endif // PCAONE_FILEPLINK_
