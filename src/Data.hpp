#ifndef PCAONE_DATA_
#define PCAONE_DATA_

#include "Utils.hpp"


class Logger
{
public:
    std::ofstream clog;
    bool isstdout = true;
    template <class S>
    Logger& operator<<(const S& val)
    {
        clog << val;
        if (isstdout) std::cout << val;
        return *this;
    }

    Logger& operator<<(std::ostream& (*pfun)(std::ostream&))
    {
        pfun(clog);
        if (isstdout) pfun(std::cout);
        return *this;
    };

    Logger(void);
    ~Logger(void);
};

class Data
{
public:
    Data(const Param& params_) : params(params_)
    {
        // start logging
        llog.clog.open(std::string(params.outfile + ".log").c_str(), std::ios::out | std::ios::trunc);
        if (params.printu) llog.isstdout = false;
    }

    virtual ~Data()
    {
    }

    virtual void read_all() = 0;
    // for blockwise
    virtual void check_file_offset_first_var() = 0;
    virtual void read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize) = 0;
    virtual void read_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize) = 0;

    void prepare(uint& blocksize);
    void standardize_E();
    void filterSNPs_resizeF();
    void pcangsd_standardize_E(const MyMatrix& U, const MyVector& svals, const MyMatrix& VT);
    void update_batch_E(const MyMatrix& U, const MyVector& svals, const MyMatrix& VT);
    void write_eigs_files(const MyVector& S, const MyMatrix& U, const MyMatrix& V);
    // for blockwise
    void calcu_vt_initial(const MyMatrix& T, MyMatrix& VT, bool standardize);
    void calcu_vt_update(const MyMatrix& T, const MyMatrix& U, const MyVector& svals, MyMatrix& VT, bool standardize);

    const Param& params;
    Logger llog;
    bool snpmajor = true;
    bool nsamples_ge_nsnps = false; // if nsamples greater than or equal to nsnps
    bool initialFonly = false;
    uint64 nsamples, nsnps; // prevent from potential integer overflow
    double readtime = 0;
    uint nblocks = 1;
    uint bandFactor = 1;
    uint nops = 0;
    std::vector<uint> start, stop;
    MyMatrix G;  // genotype matrix, can be initial E or centered E, which is nsamples x nsnps;
    MyMatrix P;  // genotype probability, nsamples x 3 x nsnps.
    MyVector F;  // observed or estimated population allele frequency
    MyVector Dc; // diagnal vector of covariance matrix
    ArrayXb C;   // boolean array indicates if a ind's snp is missing and need to be predicted.
    MyArrayX centered_geno_lookup;
    std::vector<Eigen::Index> keepSNPs; // store index of SNPs to keep
};


#endif // PCAONE_DATA_
