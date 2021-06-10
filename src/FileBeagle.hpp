#ifndef __FileBeagle__
#define __FileBeagle__

#include "Data.hpp"
#include <zlib.h>

// inspired by Angsd
int tgets(gzFile gz,char**buf,int *l);

class FileBeagle : public Data
{
public:
    //
    using Data::Data;

    ~FileBeagle() {}

    virtual void get_matrix_dimensions();
    virtual void read_all_and_centering();
    // below are for blockwise, remain for future.
    // void estimate_F() {}
    virtual void read_snp_block_initial(uint start_idx, uint stop_idx, bool standardize = false) {}
    virtual void read_snp_block_update(uint start_idx, uint stop_idx, const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT, bool standardize = false) {}


    virtual void open_check_file() {}
    virtual void close_check_file() {}

private:
    gzFile fp;
    char *original, *buffer;
    int bufsize = 128000;
    const char* delims = "\t \n";

};



#endif