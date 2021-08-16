#include "Halko.hpp"

void NormalRsvdOpData::computeGandH(MatrixXf& G, MatrixXf& H, int p)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
        cerr << "Error: the size of G or H doesn't match.\n";
        exit(EXIT_FAILURE);
    }
    if (data->params.batch)
    {
        cout << timestamp() << "running in batch mode with one-pass halko.\n";
        if (update)
        {
            data->update_batch_E(U, S, V.transpose());
        }
        if (standardize)
        {
            if (data->params.pcangsd) {
                data->pcangsd_standardize_E(U, S, V.transpose());
            } else {
                data->standardize_E();
            }
        }
        if (data->snpmajor || true) {
            for (int pi = 0; pi <= p; ++pi) {
                if (pi > 0) {
                    Eigen::HouseholderQR<Eigen::Ref<MatrixXf>> qr(H);
                    Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                }
                G.noalias() = data->G.transpose() * Omg;
                H.noalias() = data->G * G;
                stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size);
                if (stop || pi == p) {
                    cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                    print_summary_table(Upre, Ucur, data->params.outfile);
                    break;
                }
                Upre = Ucur;
            }
        }
    } else {
        // for block version
        cout << timestamp() << "running in block mode with one-pass halko.\n";
        uint actual_block_size, start_idx, stop_idx;
        // data->G is always nsamples x nsnps;
        if (data->snpmajor || true) {
            // for nsnps > nsamples
            for (int pi = 0; pi <= p ; ++pi) {
                if (pi >= 1) {
                    Eigen::HouseholderQR<Eigen::Ref<MatrixXf>> qr(H);
                    Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                }
                H = MatrixXf::Zero(cols(), size);
                data->check_file_offset_first_var();
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
                stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size);
                if (stop || pi == p) {
                    cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                    print_summary_table(Upre, Ucur, data->params.outfile);
                    break;
                }
                Upre = Ucur;
            }
        }
    }
}

void FancyRsvdOpData::computeGandH(MatrixXf& G, MatrixXf& H, int p)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
        cerr << "Error: the size of G or H doesn't match.\n";
        exit(EXIT_FAILURE);
    }
    if (data->params.batch)
    {
        cout << timestamp() << "running in batch mode with fancy halko.\n";
        if (update)
        {
            data->update_batch_E(U, S, V.transpose());
        }
        if (standardize)
        {
            if (data->params.pcangsd) {
                data->pcangsd_standardize_E(U, S, V.transpose());
            } else {
                data->standardize_E();
            }
        }
        // permute snps of G, see https://stackoverflow.com/questions/15858569/randomly-permute-rows-columns-of-a-matrix-with-eigen
        PermutationMatrix<Dynamic,Dynamic> perm(data->G.cols());
        perm.setIdentity();
        // std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
        auto rng = std::default_random_engine {};
        std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
        data->G = data->G * perm; // permute columns in-place

        uint B=data->params.bands;
        uint band=4;
        uint blocksize = ceil(data->nsnps / B);
        uint actual_block_size, start_idx, stop_idx;
        MatrixXf H1, H2, H3, H4;
        for (int pi=0; pi <= p; ++pi)
        {
            band = fmin(band * pow(2, pi), B);
            H1 = MatrixXf::Zero(cols(), size);
            H2 = MatrixXf::Zero(cols(), size);
            H3 = MatrixXf::Zero(cols(), size);
            H4 = MatrixXf::Zero(cols(), size);
            for (uint b = 0, i = 1; b < B; ++b, ++i) {
                start_idx = b * blocksize;
                stop_idx = (b + 1) * blocksize >= data->nsnps ? data->nsnps - 1 : (b + 1) * blocksize - 1 ;
                actual_block_size = stop_idx - start_idx + 1;
                G.block(start_idx, 0, actual_block_size, size) = data->G.block(start_idx, 0, data->G.rows(), actual_block_size).transpose() * Omg;
                if (i <= band / 4) {
                    H1 = H1 +  data->G.block(start_idx, 0, data->G.rows(), actual_block_size) * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > band / 4 && i <= band / 2) {
                    H2 = H2 +  data->G.block(start_idx, 0, data->G.rows(), actual_block_size) * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > band / 2 && i <= 3 * band / 4) {
                    H3 = H3 +  data->G.block(start_idx, 0, data->G.rows(), actual_block_size) * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > 3 * band / 4 && i <= band) {
                    H4 = H4 +  data->G.block(start_idx, 0, data->G.rows(), actual_block_size) * G.block(start_idx, 0, actual_block_size, size);
                }
                if( (b+1) >= band ) {
                    if (i == band) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H1 = MatrixXf::Zero(cols(), size);
                        i = 0;
                    }else if (i == band / 4) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H2 = MatrixXf::Zero(cols(), size);
                    }else if (i == band * 2 / 4) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H3 = MatrixXf::Zero(cols(), size);
                    }else if (i == band * 3 / 4) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H4 = MatrixXf::Zero(cols(), size);
                    }else if( (b+1) == data->nblocks) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                    }
                }
            }
            stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size);
            if (stop || pi == p) {
                cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                print_summary_table(Upre, Ucur, data->params.outfile);
                break;
            }
            Upre = Ucur;
        }
    } else {
        uint actual_block_size, start_idx, stop_idx;
        uint band = 4 * data->bandFactor;
        MatrixXf H1, H2, H3, H4;
        cout << timestamp() << "running in block mode with fancy halko.\n";
        for (int pi=0; pi <= p; ++pi)
        {
            band = fmin(band * pow(2, pi), data->nblocks);
            H1 = MatrixXf::Zero(cols(), size);
            H2 = MatrixXf::Zero(cols(), size);
            H3 = MatrixXf::Zero(cols(), size);
            H4 = MatrixXf::Zero(cols(), size);
            data->check_file_offset_first_var();
            for (uint b = 0, i = 1 ; b < data->nblocks ; ++b, ++i)
            {
                start_idx = data->start[b];
                stop_idx = data->stop[b];
                actual_block_size = stop_idx - start_idx + 1;
                if (update)
                {
                    data->read_snp_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
                } else {
                    data->read_snp_block_initial(start_idx, stop_idx, standardize);
                }
                G.block(start_idx, 0, actual_block_size, size) = data->G.transpose() * Omg;
                if (i <= band / 4) {
                    H1 = H1 +  data->G * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > band / 4 && i <= band / 2) {
                    H2 = H2 +  data->G * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > band / 2 && i <= 3 * band / 4) {
                    H3 = H3 +  data->G * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > 3 * band / 4 && i <= band) {
                    H4 = H4 +  data->G * G.block(start_idx, 0, actual_block_size, size);
                }
                if( (b+1) >= band ) {
                    if (i == band) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H1 = MatrixXf::Zero(cols(), size);
                        i = 0;
                    }else if (i == band / 4) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H2 = MatrixXf::Zero(cols(), size);
                    }else if (i == band * 2 / 4) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H3 = MatrixXf::Zero(cols(), size);
                    }else if (i == band * 3 / 4) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                        H4 = MatrixXf::Zero(cols(), size);
                    }else if( (b+1) == data->nblocks) {
                        H.noalias() = H1 + H2 + H3 + H4;
                        Eigen::HouseholderQR<MatrixXf> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXf::Identity(cols(), size);
                    }
                }
            }
            stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size);
            if (stop || pi == p) {
                cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                print_summary_table(Upre, Ucur, data->params.outfile);
                break;
            }
            Upre = Ucur;
        }
    }
}

bool check_if_halko_converge(int pi, double tol, MatrixXf& Upre, MatrixXf& Ucur, const MatrixXf& G, const MatrixXf& H, int k, int nrow, int ncol, int size){
    MatrixXf Q(nrow, size), B(size, ncol), R(size, size), Rt(size, size);
    Eigen::HouseholderQR<MatrixXf> qr(G);
    Q.noalias() = qr.householderQ() * MatrixXf::Identity(nrow, size);
    R.noalias() = MatrixXf::Identity(size, nrow) * qr.matrixQR().triangularView<Eigen::Upper>();
    Eigen::HouseholderQR<MatrixXf> qr2(Q);
    Q.noalias() = qr2.householderQ() * MatrixXf::Identity(nrow, size);
    Rt.noalias() = MatrixXf::Identity(size, nrow) * qr2.matrixQR().triangularView<Eigen::Upper>();
    R = Rt * R;
    // R.T * B = H.T
    B.noalias() = R.transpose().householderQr().solve(H.transpose());
    Eigen::JacobiSVD<MatrixXf> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Ucur = svd.matrixV().leftCols(k);
    if (pi == 0) {
        return false;
    } else {
        double diff_rmse, diff_mev;
        diff_mev = 1 - mev(Upre, Ucur);
        diff_rmse = rmse(Upre, Ucur);
        cout << timestamp() << "running of epoch=" << pi << ", RMSE=" << diff_rmse << ", 1-MEV=" << diff_mev << ".\n";
        if (diff_rmse <= tol || diff_mev <= tol) {
            return true;
        } else {
            return false;
        }
    }
}

void print_summary_table(const MatrixXf& Upre, const MatrixXf& Ucur, const string& outfile)
{
    VectorXd rmse = rmse_byk(Upre, Ucur);
    VectorXd mev  = mev_byk(Upre, Ucur);
    std::ofstream outlog(outfile + ".log");
    string out = "summary:";
    for (int i=0; i < Ucur.cols(); i++) {
        out += " PC1-" + std::to_string(i+1);
    }
    outlog << out << "\n"
           << "RMSE:" << rmse.transpose() << ".\n"
           << "1-MEV:" << mev.transpose() << ".\n";
}

void run_pca_with_halko(Data* data, const Param& params)
{
    cout << timestamp() << "begin to run_pca_with_halko\n";
    MatrixXf Vpre;
    VectorXf S;
    RsvdOpData* op;
    if (params.fast) {
        op = new FancyRsvdOpData(data, params.k);
    } else {
        op = new NormalRsvdOpData(data, params.k);
    }
    RsvdOnePass< MatrixXf, RsvdOpData >* rsvd = new RsvdOnePass< MatrixXf, RsvdOpData >(*op);
    if (params.maxiter == 0)
    {
        cout << timestamp() << "begin to do non-EM PCA.\n";
        op->setFlags(false, true);
        rsvd->computeUSV(params.p);
        op->U = rsvd->matrixU(data->snpmajor);
        op->V = rsvd->matrixV(data->snpmajor);
        // flip_UV(U, V, false);
        op->S = rsvd->singularValues().array().square() / data->nsnps;
        cout << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";
        data->write_eigs_files(op->S, op->U);
    } else {
        // for EM iteration
        op->setFlags(false, false);
        rsvd->computeUSV(params.p);
        op->U = rsvd->matrixU(data->snpmajor);
        op->V = rsvd->matrixV(data->snpmajor);
        op->S = rsvd->singularValues();
        // flip_UV(Upre, V, false);
        cout << timestamp() << "begin to do EM PCA.\n";
        double diff;
        op->setFlags(true, false);
        for (uint i = 0; i < params.maxiter; ++i)
        {
            Vpre = op->V;
            rsvd->computeUSV(params.p);
            op->U = rsvd->matrixU(data->snpmajor);
            op->V = rsvd->matrixV(data->snpmajor);
            op->S = rsvd->singularValues();
            // flip_UV(U, V, false);
            diff = rmse(op->V, Vpre);
            cout << timestamp() << "Individual allele frequencies estimated (iter=" << i+1 << "), RMSE=" << diff <<".\n";
            if (diff < params.tol)
            {
                cout << timestamp() << "Come to convergence!\n";
                break;
            }
        }
        cout << timestamp() << "Begin to standardize the matrix\n";
        op->setFlags(true, true);
        rsvd->computeUSV(params.p);
        op->U = rsvd->matrixU(data->snpmajor);
        op->S = rsvd->singularValues().array().square() / data->nsnps;
        // flip_UV(U, V, false);
        cout << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";
        data->write_eigs_files(op->S, op->U);

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

    delete op;
    delete rsvd;

}
