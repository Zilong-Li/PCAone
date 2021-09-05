#ifndef __EMU_HALKO__
#define __EMU_HALKO__

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData
{
public:
    using Index = Eigen::Index;
    bool update = false, standardize = false, verbose = false;
    MatrixXf U, V;
    VectorXf S;

public:
    virtual ~RsvdOpData() {}

    virtual Index rows() const = 0;
    virtual Index cols() const = 0;
    virtual Index ranks() const = 0;
    virtual Index oversamples() const = 0;

    virtual void computeGandH(MatrixXf& G, MatrixXf& H, int p) = 0;
    inline void setFlags(bool is_update, bool is_standardize, bool is_verbose)
        {update = is_update; standardize = is_standardize; verbose = is_verbose;}

};

class NormalRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    MatrixXf Omg, Upre, Ucur;
    bool stop = false;
    uint64 actual_block_size, start_idx, stop_idx;

public:

    NormalRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {
            auto rng = std::default_random_engine {};
            Omg = StandardNormalRandom<MatrixXf, std::default_random_engine>(data->nsamples, size, rng);
        }

    ~NormalRsvdOpData() {}

    Index rows() const { return data->nsnps; } // for snpmajor input
    Index cols() const { return data->nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MatrixXf& G, MatrixXf& H, int p);

};

class FancyRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    MatrixXf Omg, Upre, Ucur;
    bool stop = false;
    uint64 actual_block_size, start_idx, stop_idx;

public:

    FancyRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {
            auto rng = std::default_random_engine {};
            Omg = UniformRandom<MatrixXf, std::default_random_engine>(data->nsamples, size, rng);
        }

    ~FancyRsvdOpData() {}

    Index rows() const { return data->nsnps; }
    Index cols() const { return data->nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    void computeGandH(MatrixXf& G, MatrixXf& H, int p);
};

bool check_if_halko_converge(int pi, double tol, MatrixXf& Upre, MatrixXf& Ucur, const MatrixXf& G, const MatrixXf& H, int k, int nrow, int ncol, int size, bool verbose);
void print_summary_table(const MatrixXf& Upre, const MatrixXf& Ucur, const string& outfile);
void run_pca_with_halko(Data* data, const Param& params);

#endif
