#include "Halko.hpp"

void RsvdOpData::computeUSV(int p, double tol)
{
    const Index size{ranks() + oversamples()};
    const Index k{ranks()};
    const Index nrow{rows()};
    const Index ncol{cols()};
    MatrixXd Upre, H(ncol, size), G(nrow, size), R(size, size), Rt(size, size), B(size, ncol);
    double diff_mev;
    for (int pi=0; pi <= p; ++pi) {
        computeGandH(G, H, pi);
        // check if converged
        {
            Eigen::HouseholderQR<Ref<MatrixXd>> qr(G);
            R.noalias() = MatrixXd::Identity(size, nrow) * qr.matrixQR().triangularView<Eigen::Upper>();
            G.noalias() = qr.householderQ() * MatrixXd::Identity(nrow, size);
        }
        {
            Eigen::HouseholderQR<Ref<MatrixXd>> qr(G);
            Rt.noalias() = MatrixXd::Identity(size, nrow) * qr.matrixQR().triangularView<Eigen::Upper>();
            G.noalias() = qr.householderQ() * MatrixXd::Identity(nrow, size);
        }
        R = Rt * R;
        // B is size x ncol
        // R.T * B = H.T
        B.noalias() = R.transpose().colPivHouseholderQr().solve(H.transpose());
        Eigen::JacobiSVD<MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
        U = svd.matrixV().leftCols(k);
        if (pi > 0) {
            diff_mev = 1 - mev(Upre, U);
            verbose && cout << timestamp() << "running of epoch=" << pi << ", 1-MEV=" << diff_mev << endl;
            if (diff_mev < tol || pi == p) {
                V.noalias() = G * svd.matrixU().leftCols(k);
                S = svd.singularValues().head(k);
                verbose && cout << timestamp() << "stops at epoch=" << pi + 1 << ".\n";
                break;
            } else {
                Upre = U;
            }
        } else {
            Upre = U;
        }
    }
}

void NormalRsvdOpData::computeGandH(MatrixXd& G, MatrixXd& H, int pi)
{
    // check size of G and H first;
    if (pi == 0) {
        auto rng = std::default_random_engine {};
        Omg = StandardNormalRandom<MatrixXd, std::default_random_engine>(data->nsamples, size, rng);
    }
    if (data->params.batch)
    {
        if (pi == 0) {
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
        }
        if (data->snpmajor || true) {
            if (pi > 0) {
                Eigen::HouseholderQR<Eigen::Ref<MatrixXd>> qr(H);
                Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
            }
            G.noalias() = data->G.transpose() * Omg;
            H.noalias() = data->G * G;
        }
    } else {
        // for block version
        if (pi==0) verbose && cout << timestamp() << "running in blockwise mode with one-pass halko.\n";
        // data->G is always nsamples x nsnps;
        if (data->snpmajor || true) {
            // for nsnps > nsamples
            if (pi > 0) {
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
                G.middleRows(start_idx, actual_block_size).noalias() = data->G.transpose() * Omg;
                H.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
            }
        }
    }
}

void FancyRsvdOpData::computeGandH(MatrixXd& G, MatrixXd& H, int pi)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows()) {
        throw std::runtime_error("Error: the size of G or H doesn't match.\n");
    }
    MatrixXd H1, H2;
    if (pi == 0) {
        auto rng = std::default_random_engine {};
        Omg = UniformRandom<MatrixXd, std::default_random_engine>(data->nsamples, size, rng);
        Omg2 = Omg;
    }
    if (data->params.batch)
    {
        if (pi == 0) {
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
            band = 2;
            blocksize = (unsigned int)ceil((double)data->nsnps / data->params.bands);

        }
        {
            // band : 4, 8, 16, 32, 64, 128
            band = fmin(band * 2, data->params.bands);
            H1 = MatrixXd::Zero(cols(), size);
            H2 = MatrixXd::Zero(cols(), size);
            for (uint b = 0, i = 1; b < data->params.bands; ++b, ++i) {
                start_idx = b * blocksize;
                stop_idx = (b + 1) * blocksize >= data->nsnps ? data->nsnps - 1 : (b + 1) * blocksize - 1 ;
                actual_block_size = stop_idx - start_idx + 1;
                G.middleRows(start_idx, actual_block_size).noalias() = data->G.middleCols(start_idx, actual_block_size).transpose() * Omg;
                if (i <= band / 2) {
                    H1.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                }else if (i > band / 2 && i <= band) {
                    H2.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                }
                if( (b+1) >= band ) {
                    if (i == band) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H1 = MatrixXd::Zero(cols(), size);
                        i = 0;
                    }else if (i == band / 2) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H2 = MatrixXd::Zero(cols(), size);
                    }else if( (b+1) == data->nblocks) {
                        // shouldn't go here if the bands is proper.
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                    }
                }
            }
        }
    } else {
        if (pi == 0) {
            verbose && cout << timestamp() << "running in blockwise mode with fancy halko.\n";
            band = 2 * data->bandFactor;
        }
        {
            // band : 4, 8, 16, 32, 64, 128
            band = fmin(band * 2, data->nblocks);
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
                G.middleRows(start_idx, actual_block_size).noalias() = data->G.transpose() * Omg;
                if (i <= band / 2) {
                    H1.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
                }else if (i > band / 2 && i <= band) {
                    H2.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
                }
                if( (b+1) >= band ) {
                    if (i == band) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H1 = MatrixXd::Zero(cols(), size);
                        i = 0;
                    }else if (i == band / 2) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H2 = MatrixXd::Zero(cols(), size);
                    }else if( (b+1) == data->nblocks) {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MatrixXd> qr(H);
                        Omg.noalias() = qr.householderQ() * MatrixXd::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                    }
                }
            }
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
    RsvdOpData* rsvd;
    if (!params.runem && params.fast) {
        rsvd = new FancyRsvdOpData(data, params.k, params.oversamples);
    } else {
        rsvd = new NormalRsvdOpData(data, params.k, params.oversamples);
    }
    if (!params.runem)
    {
        cout << timestamp() << "begin to do non-EM PCA.\n";
        rsvd->setFlags(false, true, params.verbose);
        rsvd->computeUSV(params.maxp, params.tol_halko);
        data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->U, rsvd->V);
    } else {
        // for EM iteration
        rsvd->setFlags(false, false, false);
        rsvd->computeUSV(params.maxp, params.tol_halko);
        // flip_UV(Upre, V, false);
        double diff;
        rsvd->setFlags(true, false, false);
        cout << timestamp() << "begin to do EM PCA.\n";
        for (uint i = 0; i < params.maxiter; ++i)
        {
            Vpre = rsvd->V;
            rsvd->computeUSV(params.maxp, params.tol_halko);
            // flip_UV(U, V, false);
            diff = rmse(rsvd->V, Vpre);
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
            data->pcangsd_standardize_E(rsvd->U, rsvd->S, rsvd->V.transpose());
            MatrixXd C = data->G * data->G.transpose();
            C.array() /= (double) data->nsnps;
            C.diagonal() = data->Dc.array() / (double) data->nsnps;
            std::ofstream out_cov(params.outfile + ".cov");
            if (out_cov.is_open()) {
                out_cov << C << "\n";
            }
            // use Eigen::JacobiSVD to get eigenvecs
            Eigen::JacobiSVD< MatrixXd > svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
            data->write_eigs_files(svd.singularValues().head(params.k), svd.matrixU().leftCols(params.k), svd.matrixU().leftCols(params.k));
        } else {
            rsvd->setFlags(true, true, false);
            rsvd->computeUSV(params.maxp, params.tol_halko);
            data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->U, rsvd->V);
        }
    }

    delete rsvd;

}
