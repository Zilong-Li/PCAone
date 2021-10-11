#ifndef __EMU_HALKO__
#define __EMU_HALKO__

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData
{
public:
    using Index = Eigen::Index;
    bool update = false, standardize = false, verbose = false;
    MatrixXd U, V;
    VectorXd S;

public:
    virtual ~RsvdOpData() {}

    virtual Index rows() const = 0;
    virtual Index cols() const = 0;
    virtual Index ranks() const = 0;
    virtual Index oversamples() const = 0;

    virtual void computeGandH(MatrixXd& G, MatrixXd& H, int p) = 0;
    inline void setFlags(bool is_update, bool is_standardize, bool is_verbose)
        {update = is_update; standardize = is_standardize; verbose = is_verbose;}

};

class NormalRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    bool stop = false;
    uint64 actual_block_size, start_idx, stop_idx;

public:

    NormalRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {}

    ~NormalRsvdOpData() {}

    Index rows() const { return data->nsnps; } // for snpmajor input
    Index cols() const { return data->nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MatrixXd& G, MatrixXd& H, int p);

};

class FancyRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    bool stop = false;
    uint64 actual_block_size, start_idx, stop_idx;

public:

    FancyRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {}

    ~FancyRsvdOpData() {}

    Index rows() const { return data->nsnps; }
    Index cols() const { return data->nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MatrixXd& G, MatrixXd& H, int p);
};

bool check_if_halko_converge(int pi, double tol, MatrixXd& Upre, MatrixXd& Ucur, MatrixXd& G, MatrixXd& H, int k, int nrow, int ncol, int size, bool verbose);
void print_summary_table(const MatrixXd& Upre, const MatrixXd& Ucur);
void run_pca_with_halko(Data* data, const Param& params);

#endif
