#ifndef __EMU_HALKO__
#define __EMU_HALKO__

#include "Data.hpp"
#include "RSVD.hpp"

// note:for blockwise algorithm, nsnps must be larger than nsamples, ie. matrix should be tall or it may lose some accuracy
class NormalRsvdOpData : public RsvdOpOnePass<MatrixXf>
{
private:
    using Index = Eigen::Index;

    Data* data;
    bool batch = false, update = false, standardize = false, pcangsd = false;
    const Index nsnps, nsamples, nk, os, size;
    MatrixXf Omg;

public:

    NormalRsvdOpData(Data* data_, bool batch_, int k_, int os_ = 10) :
        data(data_), batch(batch_), nsnps(data_->nsnps), nsamples(data_->nsamples), nk(k_), os(os_), size(k_ + os_)
        {
            std::mt19937_64 randomEngine{};
            randomEngine.seed(111);
            Omg = StandardNormalRandom<MatrixXf, std::mt19937_64>(nsamples, size, randomEngine);
        }

    ~NormalRsvdOpData() {}

    Index rows() const { return nsnps; }
    Index cols() const { return nsamples; }
    Index ranks() const { return nk; }
    Index oversamples() const { return os; }

    inline void setFlags(bool is_update, bool is_standardize, bool is_pcangsd) {update = is_update; standardize = is_standardize; pcangsd = is_pcangsd;}

    void computeGandH(MatrixXf& G, MatrixXf& H, int p);

    MatrixXf U, V;
    VectorXf S;
};


void run_pca_with_halko(Data* data, const Param& params);

#endif
