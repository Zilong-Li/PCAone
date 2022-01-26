#ifndef __FileBeagle__
#define __FileBeagle__

#include "Data.hpp"
#include <zlib.h>

int tgets(gzFile gz,char**buf,int *l);

class FileBeagle : public Data
{
public:
    FileBeagle(const Param& params_) : Data(params_)
        {
            fp = gzopen(params.beagle.c_str(), "r");
            original = buffer =(char *) calloc(bufsize, sizeof(char));
            tgets(fp, &buffer, &bufsize);
            int nCol = 1;
            if(buffer!=original) original=buffer;
            strtok_r(buffer,delims,&buffer);
            while(strtok_r(NULL,delims,&buffer))
                nCol++;
            if(nCol % 3){
                throw std::runtime_error("Number of columns should be a multiple of 3.\n");
            }
            nsamples = nCol/3-1;
            // continue getting the number of sites
            // assume the number of columns of each line is the same. should check it first.
            buffer = original;
            int nSites = 0;
            while(tgets(fp, &buffer, &bufsize)) {
                nSites++;
            }
            nsnps = nSites;
            gzclose(fp);

            std::cout << timestamp() << "N samples is " << nsamples << ". M snps is " << nsnps << std::endl;
            P = MatrixXd(nsnps, nsamples * 3);
        }

    ~FileBeagle() {}

    virtual void read_all_and_centering();
    // below are for blockwise, remain for future.
    virtual void read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize = false) {}
    virtual void read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MatrixXd& U, const VectorXd& svals, const MatrixXd& VT, bool standardize = false) {}


    virtual void check_file_offset_first_var() {}

private:
    gzFile fp;
    char *original, *buffer;
    int bufsize = 128000;
    const char* delims = "\t \n";

};



#endif
