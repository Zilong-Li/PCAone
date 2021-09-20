#include "Halko.hpp"

void NormalRsvdOpData::computeGandH(MatrixXd& G, MatrixXd& H, int p)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
        throw std::runtime_error("Error: the size of G or H doesn't match.\n");
    }
    if (data->params.batch)
    {
        verbose && cout << timestamp() << "running in batch mode with one-pass halko.\n";
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
                    Eigen::HouseholderQR<Eigen::Ref<MatrixXd>> qr(H);
                    Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                }
                G.noalias() = data->G.transpose() * Omg;
                H.noalias() = data->G * G;
                stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size, verbose);
                if (stop || pi == p) {
                    if (verbose) {
                        cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                        print_summary_table(Upre, Ucur);
                    }
                    break;
                }
                Upre = Ucur;
            }
        }
    } else {
        // for block version
        verbose && cout << timestamp() << "running in blockwise mode with one-pass halko.\n";
        // data->G is always nsamples x nsnps;
        if (data->snpmajor || true) {
            // for nsnps > nsamples
            for (int pi = 0; pi <= p ; ++pi) {
                if (pi >= 1) {
                    Eigen::HouseholderQR<Eigen::Ref<MatrixXd>> qr(H);
                    Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                }
                H = MatrixXd::Zero(cols(), size);
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
                    G.block(start_idx, 0, actual_block_size, size).noalias() = data->G.transpose() * Omg;
                    H.noalias() += data->G * G.block(start_idx, 0, actual_block_size, size);
                }
                stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size, verbose);
                if (stop || pi == p) {
                    if (verbose) {
                        cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                        print_summary_table(Upre, Ucur);
                    }
                    break;
                }
                Upre = Ucur;
            }
        }
    }
}

void FancyRsvdOpData::computeGandH(MatrixXd& G, MatrixXd& H, int p)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
        throw std::runtime_error("Error: the size of G or H doesn't match.\n");
    }
    if (data->params.batch)
    {
        verbose && cout << timestamp() << "running in batch mode with fancy halko.\n";
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
        auto rng = std::default_random_engine {};
        std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
        data->G = data->G * perm; // permute columns in-place

        uint band=4;  // maybe 2
        uint blocksize = (unsigned int)ceil((double)data->nsnps / data->params.bands);
        for (int pi=0; pi <= p; ++pi)
        {
            // 4, 8, 32, 256
            band = fmin(band * pow(2, pi), data->params.bands);
            H1 = MatrixXd::Zero(cols(), size);
            H2 = MatrixXd::Zero(cols(), size);
            for (uint b = 0, i = 1; b < data->params.bands; ++b, ++i) {
                start_idx = b * blocksize;
                stop_idx = (b + 1) * blocksize >= data->nsnps ? data->nsnps - 1 : (b + 1) * blocksize - 1 ;
                actual_block_size = stop_idx - start_idx + 1;
                G.block(start_idx, 0, actual_block_size, size).noalias() = data->G.block(0, start_idx, data->G.rows(), actual_block_size).transpose() * Omg;
                if (i <= band / 2) {
                    H1.noalias() += data->G.block(0, start_idx, data->G.rows(), actual_block_size) * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > band / 2 && i <= band) {
                    H2.noalias() += data->G.block(0, start_idx, data->G.rows(), actual_block_size) * G.block(start_idx, 0, actual_block_size, size);
                }
                if( (b+1) >= band ) {
                    if (i == band) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Y(Omg2, Omg);
                        Omg2 = Omg;
                        H1 = MatrixXd::Zero(cols(), size);
                        i = 0;
                    }else if (i == band / 2) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Y(Omg2, Omg);
                        Omg2 = Omg;
                        H2 = MatrixXd::Zero(cols(), size);
                    }else if( (b+1) == data->nblocks) {
                        // shouldn't go here if the bands is proper.
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Y(Omg2, Omg);
                        Omg2 = Omg;
                    }
                }
            }
            // band = fmin(band * 2, data->params.bands);
            stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size, verbose);
            if (stop || pi == p) {
                if (verbose) {
                    cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                    print_summary_table(Upre, Ucur);
                }
                break;
            }
            Upre = Ucur;
        }
    } else {
        uint band = 4 * data->bandFactor;
        verbose && cout << timestamp() << "running in blockwise mode with fancy halko.\n";
        for (int pi=0; pi <= p; ++pi)
        {
            // band : 4, 8, 32, 256
            band = fmin(band * pow(2, pi), data->nblocks);
            H1 = MatrixXd::Zero(cols(), size);
            H2 = MatrixXd::Zero(cols(), size);
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
                G.block(start_idx, 0, actual_block_size, size).noalias() = data->G.transpose() * Omg;
                if (i <= band / 2) {
                    H1.noalias() += data->G * G.block(start_idx, 0, actual_block_size, size);
                }else if (i > band / 2 && i <= band) {
                    H2.noalias() += data->G * G.block(start_idx, 0, actual_block_size, size);
                }
                if( (b+1) >= band ) {
                    if (i == band) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Y(Omg2, Omg);
                        Omg2 = Omg;
                        H1 = MatrixXd::Zero(cols(), size);
                        i = 0;
                    }else if (i == band / 2) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Y(Omg2, Omg);
                        Omg2 = Omg;
                        H2 = MatrixXd::Zero(cols(), size);
                    }else if( (b+1) == data->nblocks) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Y(Omg2, Omg);
                        Omg2 = Omg;
                    }
                }
            }
            stop = check_if_halko_converge(pi, data->params.tol_halko, Upre, Ucur, G, H, nk, rows(), cols(), size, verbose);
            if (stop || pi == p) {
                if (verbose) {
                    cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                    print_summary_table(Upre, Ucur);
                }
                break;
            }
            Upre = Ucur;
        }
    }
}

bool check_if_halko_converge(int pi, double tol, MatrixXd& Upre, MatrixXd& Ucur, const MatrixXd& G, const MatrixXd& H, int k, int nrow, int ncol, int size, bool verbose = false)
{
    MatrixXd Q(nrow, size), B(size, ncol), R(size, size), Rt(size, size);
    Eigen::HouseholderQR<MatrixXd> qr(G);
    Q.noalias() = qr.householderQ() * MatrixXd::Identity(nrow, size);
    R.noalias() = MatrixXd::Identity(size, nrow) * qr.matrixQR().triangularView<Eigen::Upper>();
    Eigen::HouseholderQR<MatrixXd> qr2(Q);
    Q.noalias() = qr2.householderQ() * MatrixXd::Identity(nrow, size);
    Rt.noalias() = MatrixXd::Identity(size, nrow) * qr2.matrixQR().triangularView<Eigen::Upper>();
    R = Rt * R;
    // R.T * B = H.T
    B.noalias() = R.transpose().householderQr().solve(H.transpose());
    Eigen::JacobiSVD<MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Ucur = svd.matrixV().leftCols(k);
    if (pi == 0) {
        return false;
    } else {
        double diff_rmse, diff_mev;
        diff_mev = 1 - mev(Upre, Ucur);
        diff_rmse = rmse(Upre, Ucur);
        verbose && cout << timestamp() << "running of epoch=" << pi << ", RMSE=" << diff_rmse << ", 1-MEV=" << diff_mev << ".\n";
        if (diff_rmse <= tol || diff_mev <= tol) {
            return true;
        } else {
            return false;
        }
    }
}

void print_summary_table(const MatrixXd& Upre, const MatrixXd& Ucur)
{
    string out = "summary:";
    for (int i=0; i < Ucur.cols(); i++) {
        out += " PC1-" + std::to_string(i+1);
    }
    VectorXd Vrmse = VectorXd::Zero(Ucur.cols());
    VectorXd Vmev  = VectorXd::Zero(Ucur.cols());
    mev_rmse_byk(Upre, Ucur, Vmev, Vrmse);
    cout   << out << "\n"
           << "RMSE:" << Vrmse.transpose() << ".\n"
           << "1-MEV:" << Vmev.transpose() << ".\n";
}

void run_pca_with_halko(Data* data, const Param& params)
{
    if (params.batch) {
        cout << timestamp() << "begin to run_pca_with_halko batch mode\n";
    } else {
        cout << timestamp() << "begin to run_pca_with_halko blockwise mode\n";
    }
    MatrixXd Vpre;
    VectorXd S;
    RsvdOpData* op;
    if (!params.runem && params.fast) {
        op = new FancyRsvdOpData(data, params.k, fmax(params.k, params.oversamples));
    } else {
        op = new NormalRsvdOpData(data, params.k, fmax(params.k, params.oversamples));
    }
    RsvdOnePass< MatrixXd, RsvdOpData >* rsvd = new RsvdOnePass< MatrixXd, RsvdOpData >(*op);
    if (!params.runem)
    {
        cout << timestamp() << "begin to do non-EM PCA.\n";
        op->setFlags(false, true, params.verbose);
        rsvd->computeUSV(params.p);
        op->U = rsvd->matrixU(data->snpmajor);
        op->V = rsvd->matrixV(data->snpmajor);
        // flip_UV(U, V, false);
        op->S = rsvd->singularValues().array().square() / data->nsnps;
        data->write_eigs_files(op->S, op->U);
    } else {
        // for EM iteration
        op->setFlags(false, false, false);
        rsvd->computeUSV(params.p);
        op->U = rsvd->matrixU(data->snpmajor);
        op->V = rsvd->matrixV(data->snpmajor);
        op->S = rsvd->singularValues();
        // flip_UV(Upre, V, false);
        double diff;
        op->setFlags(true, false, false);
        cout << timestamp() << "begin to do EM PCA.\n";
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
        cout << timestamp() << "Begin to standardize the matrix.\n";

        // if pcangsd, estimate GRM.
        if (params.pcangsd) {
            data->pcangsd_standardize_E(op->U, op->S, op->V.transpose());
            MatrixXd C = data->G * data->G.transpose();
            C.array() /= (double) data->nsnps;
            C.diagonal() = data->Dc.array() / (double) data->nsnps;
            std::ofstream out_cov(params.outfile + ".cov");
            if (out_cov.is_open()) {
                out_cov << C << "\n";
            }
            // use Eigen::JacobiSVD to get eigenvecs
            Eigen::JacobiSVD< MatrixXd > svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
            data->write_eigs_files(svd.singularValues().head(params.k), svd.matrixU().leftCols(params.k));
        } else {
            op->setFlags(true, true, false);
            rsvd->computeUSV(params.p);
            op->U = rsvd->matrixU(data->snpmajor);
            op->S = rsvd->singularValues().array().square() / data->nsnps;
            // flip_UV(U, V, false);
            data->write_eigs_files(op->S, op->U);
        }
    }

    delete op;
    delete rsvd;

}
