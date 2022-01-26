#ifndef __FilePlink__
#define __FilePlink__

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
            std::string fbim = params.bed_prefix + ".bim";
            std::string ffam = params.bed_prefix + ".fam";
            nsamples = count_lines(ffam);
            nsnps = count_lines(fbim);
            snpmajor = true;
            std::cout << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
            bed_bytes_per_snp = (nsamples+3)>>2; // get ceiling(nsamples/4)
            std::string fbed = params.bed_prefix + ".bed";
            bed_ifstream.open(fbed, std::ios::in | std::ios::binary);
            if (!bed_ifstream.is_open()) {
                throw std::invalid_argument("ERROR: Cannot open bed file.\n");
            }
            // check magic number of bed file
            uchar header[3];
            bed_ifstream.read( reinterpret_cast<char *> (&header[0]), 3);
            if ( (header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01) ) {
                throw std::invalid_argument("ERROR: Incorrect magic number in plink bed file.\n");
            }
            
        }

    ~FileBed() {}

    virtual void read_all_and_centering();
    // for blockwise
    virtual void read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false);
    virtual void read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize = false);


    virtual void check_file_offset_first_var();

private:
    std::ifstream bed_ifstream;
    uint64 bed_bytes_per_snp;
    bool frequency_was_estimated = false;
    std::vector<uchar> inbed;
};

#endif
