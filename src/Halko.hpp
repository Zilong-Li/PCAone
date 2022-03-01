#ifndef __HALKO__
#define __HALKO__

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData
{
public:
    using Index = Eigen::Index;
    bool update = false, standardize = false, verbose = false;
    MyMatrix U, V;
    MyVector S;

public:
    virtual ~RsvdOpData() {}

    virtual Index rows() const = 0;
    virtual Index cols() const = 0;
    virtual Index ranks() const = 0;
    virtual Index oversamples() const = 0;

    virtual void computeGandH(MyMatrix& G, MyMatrix& H, int pi) = 0;

    inline void setFlags(bool is_update, bool is_standardize, bool is_verbose)
        {update = is_update; standardize = is_standardize; verbose = is_verbose;}

    void computeUSV(int p, double tol);

};


class NormalRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    uint64 actual_block_size, start_idx, stop_idx;
    MyMatrix Omg;
public:

    NormalRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {}

    ~NormalRsvdOpData() {}

    Index rows() const { return data->nsnps; } // for snpmajor input
    Index cols() const { return data->nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MyMatrix& G, MyMatrix& H, int pi=0);

};

class FancyRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    uint64 band, blocksize, actual_block_size, start_idx, stop_idx;
    MyMatrix Omg, Omg2;

public:

    FancyRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {}

    ~FancyRsvdOpData() {}

    Index rows() const { return data->nsnps; }
    Index cols() const { return data->nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MyMatrix& G, MyMatrix& H, int pi = 0);
};

// void print_summary_table(const MyMatrix& Upre, const MyMatrix& Ucur);
void run_pca_with_halko(Data* data, const Param& params);

#endif
