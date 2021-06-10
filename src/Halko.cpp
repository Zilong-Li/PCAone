#include "Halko.hpp"

void NormalRsvdOpData::computeGandH(MatrixXf& G, MatrixXf& H, int p)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
        cerr << "Error: the size of G or H doesn't match.\n";
        exit(EXIT_FAILURE);
    }
    if (batch)
    {
        if (update)
        {
            data->update_batch_E(U, S, V.transpose());
        }
        if (standardize)
        {
            if (pcangsd) {
                data->pcangsd_standardize_E(U, S, V.transpose());
            } else {
                data->standardize_E();
            }
        }
        G.noalias() = data->G.transpose() * Omg;
        H.noalias() = data->G * G;
        if (p > 0) {
            for (int i=0; i<p; ++i) {
                // Eigen::ColPivHouseholderQR<Eigen::Ref<MatrixXf>> qr(H);
                Eigen::HouseholderQR<Eigen::Ref<MatrixXf>> qr(H);
                H.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                G.noalias() = data->G.transpose() * H;
                H.noalias() = data->G * G;
            }
        }
    } else {
        // the bigger size of block is, the faster program is
        uint actual_block_size, start_idx, stop_idx;
        H = MatrixXf::Zero(cols(), size);
        data->open_check_file();
        for (uint i = 0 ; i < data->nblocks ; ++i)
        {
            start_idx = data->start[i];
            stop_idx = data->stop[i];
            actual_block_size = stop_idx - start_idx + 1;
            if (update)
            {
                data->read_snp_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
            } else {
                data->read_snp_block_initial(start_idx, stop_idx, standardize);
            }
            G.block(start_idx, 0, actual_block_size, size) = data->G.transpose() * Omg;
            H.noalias() = H +  data->G * G.block(start_idx, 0, actual_block_size, size);
        }
        data->close_check_file();
        if (p > 0)
        {
            MatrixXf Hpi;
            for (int pi=0; pi < p; ++pi)
            {
                Hpi = H;
                // Eigen::ColPivHouseholderQR<Eigen::Ref<MatrixXf>> qr(Hpi);
                Eigen::HouseholderQR<Eigen::Ref<MatrixXf>> qr(Hpi);
                Hpi.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                H = MatrixXf::Zero(cols(), size);
                data->open_check_file();
                for (uint i = 0 ; i < data->nblocks ; ++i)
                {
                    start_idx = data->start[i];
                    stop_idx = data->stop[i];
                    actual_block_size = stop_idx - start_idx + 1;
                    if (update)
                    {
                        data->read_snp_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
                    } else {
                        data->read_snp_block_initial(start_idx, stop_idx, standardize);
                    }
                    G.block(start_idx, 0, actual_block_size, size) = data->G.transpose() * Hpi;
                    H.noalias() = H +  data->G * G.block(start_idx, 0, actual_block_size, size);
                }
                data->close_check_file();
            }
        }
    }
}

void run_pca_with_halko(Data* data, const Param& params)
{
    cerr << timestamp() << "begin to run_pca_with_halko\n";
    MatrixXf Vpre;
    VectorXf S;
    NormalRsvdOpData op(data, params.batch, params.k);
    RsvdOnePass< MatrixXf, NormalRsvdOpData > rsvd(op);
    // if (!params.truefile.empty()) {
    //     TM = load_csv<MatrixXf>(params.truefile);
    // }
    if (params.maxiter == 0)
    {
        cerr << timestamp() << "begin to do non-em PCA.\n";
        op.setFlags(false, true, false);
        rsvd.computeUSV(params.p);
        op.V = rsvd.matrixU();
        op.U = rsvd.matrixV();
        // flip_UV(U, V, false);
        op.S = rsvd.singularValues().array().square() / data->nsnps;
        cerr << timestamp() << "begin to save eigenvecs and eigenvals.\n";
        data->write_eigs_files(op.S, op.U);
    } else {
        rsvd.computeUSV(params.p);
        op.V = rsvd.matrixU();
        op.U = rsvd.matrixV();
        op.S = rsvd.singularValues();
        // flip_UV(Upre, V, false);
        cerr << timestamp() << "begin to do EM\n";
        double diff;
        // std::ofstream out_vecs(params.outfile + ".eigvecs.alliters");
        op.setFlags(true, false, params.pcangsd);
        for (uint i = 0; i < params.maxiter; ++i)
        {
            Vpre = op.V;
            rsvd.computeUSV(params.p);
            op.V = rsvd.matrixU();
            op.U = rsvd.matrixV();
            op.S = rsvd.singularValues();
            // flip_UV(U, V, false);
            diff = rmse(op.V, Vpre);
            cerr << timestamp() << "Individual allele frequencies estimated (iter=" << i+1 << "), RMSE=" << diff <<".\n";
            if (diff < params.tol)
            {
                cerr << timestamp() << "Come to convergence!\n";
                break;
            }
        }
        cerr << timestamp() << "Begin to standardize the matrix\n";
        op.setFlags(true, true, params.pcangsd);
        rsvd.computeUSV(params.p);
        op.U = rsvd.matrixV();
        op.S = rsvd.singularValues().array().square() / data->nsnps;
        // flip_UV(U, V, false);
        cerr << timestamp() << "begin to save eigenvecs and eigenvals.\n";
        data->write_eigs_files(op.S, op.U);

        // if pcangsd, estimate GRM.
        if (params.pcangsd) {
            MatrixXf C = data->G * data->G.transpose();
            C.array() /= (double) data->nsnps;
            C.diagonal() = data->Dc.array() / (double) data->nsnps;
            std::ofstream out_cov(params.outfile + ".cov");
            if (out_cov.is_open()) {
                out_cov << C << "\n";
            }
        }
    }

}
