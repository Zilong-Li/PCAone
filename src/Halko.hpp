#ifndef __EMU_HALKO__
#define __EMU_HALKO__

#include "Data.hpp"
#include "RSVD.hpp"

class RsvdOpData
{
public:
    using Index = Eigen::Index;
    bool update = false, standardize = false;
    MatrixXf U, V;
    VectorXf S;

public:
    virtual ~RsvdOpData() {}

    virtual Index rows() const = 0;
    virtual Index cols() const = 0;
    virtual Index ranks() const = 0;
    virtual Index oversamples() const = 0;

    virtual void computeGandH(MatrixXf& G, MatrixXf& H, int p) = 0;
    inline void setFlags(bool is_update, bool is_standardize)
        {update = is_update; standardize = is_standardize;}

};

// halko should work on tall matrix otherwise it may lose accuracy.
// for normal halko we always make the matrix tall, ie. nrow > ncol
class NormalRsvdOpData : public RsvdOpData
{
private:
    Data* data;
    const Index nk, os, size;
    MatrixXf Omg, Upre, Ucur;
    bool stop = false;

public:

    NormalRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {
            std::mt19937_64 randomEngine{};
            randomEngine.seed(111);
            Omg = StandardNormalRandom<MatrixXf, std::mt19937_64>(data->nsamples, size, randomEngine);
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

public:

    FancyRsvdOpData(Data* data_, int k_, int os_ = 10) :
        data(data_), nk(k_), os(os_), size(k_ + os_)
        {
            std::mt19937_64 randomEngine{};
            randomEngine.seed(111);
            Omg = StandardNormalRandom<MatrixXf, std::mt19937_64>(data->nsamples, size, randomEngine);
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
