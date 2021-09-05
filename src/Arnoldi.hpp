#ifndef __EMU_Arnoldi__
#define __EMU_Arnoldi__

#include "Data.hpp"

class ArnoldiOpData
{
public:
    ArnoldiOpData(Data* data_) : data(data_), n(data_->nsamples)
        { nops = 1; }

    ~ArnoldiOpData() {}

    inline uint64 rows() const { return n; }
    inline uint64 cols() const { return n; }
    // y = G * G' * x ; data.G is n x m;
    void perform_op(const float *x_in, float* y_out);
    inline void setFlags(bool is_update, bool is_standardize, bool is_pcangsd) {update = is_update; standardize = is_standardize; pcangsd = is_pcangsd;}

    MatrixXf U, VT;
    VectorXf S;
    uint nops;

private:
    Data* data;
    const uint64 n;
    bool update = false, standardize = false, pcangsd = false;
};


void run_pca_with_arnoldi(Data* data, const Param& params);

#endif
