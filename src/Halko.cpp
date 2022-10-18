#include "Halko.hpp"

using namespace std;

void RsvdOpData::computeUSV(int p, double tol)
{
    const Index size{ranks() + oversamples()};
    const Index k{ranks()};
    const Index nrow{rows()};
    const Index ncol{cols()};
    MyMatrix Upre, H(ncol, size), G(nrow, size), R(size, size), Rt(size, size), B(size, ncol);
    double diff;
    for (int pi = 0; pi <= p; ++pi)
    {
        computeGandH(G, H, pi);
        // check if converged
        {
            Eigen::HouseholderQR<Eigen::Ref<MyMatrix>> qr(G);
            R.noalias() = MyMatrix::Identity(size, nrow) * qr.matrixQR().triangularView<Eigen::Upper>(); // get R1
            G.noalias() = qr.householderQ() * MyMatrix::Identity(nrow, size);                            // hold Q1 in G
        }
        {
            Eigen::HouseholderQR<Eigen::Ref<MyMatrix>> qr(G);
            Rt.noalias() = MyMatrix::Identity(size, nrow) * qr.matrixQR().triangularView<Eigen::Upper>(); // get R2
            G.noalias() = qr.householderQ() * MyMatrix::Identity(nrow, size);                             // hold Q2 in G
        }
        R = Rt * R; // get R = R1R2;
        // B is size x ncol
        // R.T * B = H.T
        B.noalias() = R.transpose().fullPivHouseholderQr().solve(H.transpose());
        Eigen::JacobiSVD<MyMatrix> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
        U = svd.matrixV().leftCols(k);
        if (data->params.printu)
        {
            std::ofstream ulog(std::string(data->params.outfile + ".epoch." + std::to_string(pi) + ".eigvecs").c_str());
            ulog << U;
        }
        if (pi > 0)
        {
            if (data->params.mev)
                diff = 1 - mev(U, Upre);
            else
                diff = minSSE(U, Upre).sum() / Upre.cols();
            if (verbose)
                data->llog << timestamp() << "running of epoch=" << pi << ", diff=" << diff << endl;
            if (diff < tol || pi == p)
            {
                V.noalias() = G * svd.matrixU().leftCols(k);
                S = svd.singularValues().head(k);
                if (verbose)
                    data->llog << timestamp() << "stops at epoch=" << pi + 1 << endl;
                if (data->params.fast && std::pow(2, pi + 1) < data->params.bands)
                    data->llog
                        << timestamp() << colwarn
                        << "the algorithm2 converged but the window size is not reaching the whole data. the eigen values doesn't mean anything but the PCs still reflect the structure." +
                               colend
                        << endl;
                break;
            }
            else
            {
                Upre = U;
            }
        }
        else
        {
            Upre = U;
        }
    }
}

void NormalRsvdOpData::computeGandH(MyMatrix& G, MyMatrix& H, int pi)
{
    // check size of G and H first;
    if (pi == 0)
    {
        auto rng = std::default_random_engine{};
        if (data->params.rand)
            Omg = StandardNormalRandom<MyMatrix, std::default_random_engine>(data->nsamples, size, rng);
        else
            Omg = UniformRandom<MyMatrix, std::default_random_engine>(data->nsamples, size, rng);
    }
    if (data->params.batch)
    {
        if (pi == 0)
        {
            if (verbose)
                data->llog << timestamp() << "running in batch mode with one-pass rsvd." << endl;
            if (update)
            {
                data->update_batch_E(U, S, V.transpose());
            }
            if (standardize)
            {
                if (data->params.pcangsd)
                {
                    data->pcangsd_standardize_E(U, S, V.transpose());
                }
                else
                {
                    data->standardize_E();
                }
            }
        }
        if (data->snpmajor || true)
        { // only work with snpmajor input data now.
            if (pi > 0)
            {
                Eigen::HouseholderQR<Eigen::Ref<MyMatrix>> qr(H);
                Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size); // hold H in Omega
            }
            G.noalias() = data->G.transpose() * Omg;
            H.noalias() = data->G * G;
        }
    }
    else
    {
        // for block version
        if (pi == 0 && verbose)
            data->llog << timestamp() << "running in blockwise mode with one-pass rsvd." << endl;
        // data->G is always nsamples x nsnps;
        if (data->snpmajor || true)
        {
            // for nsnps > nsamples
            if (pi > 0)
            {
                Eigen::HouseholderQR<Eigen::Ref<MyMatrix>> qr(H);
                Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
            }
            H = MyMatrix::Zero(cols(), size);
            data->check_file_offset_first_var();
            for (int i = 0; i < data->nblocks; ++i)
            {
                start_idx = data->start[i];
                stop_idx = data->stop[i];
                actual_block_size = stop_idx - start_idx + 1;
                auto t1 = std::chrono::steady_clock::now();
                if (update)
                {
                    data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
                }
                else
                {
                    data->read_block_initial(start_idx, stop_idx, standardize);
                }
                auto t2 = std::chrono::steady_clock::now();
                data->readtime +=
                    std::chrono::duration<double>(t2 - t1).count() * std::chrono::duration<double>::period::num / std::chrono::duration<double>::period::den;
                G.middleRows(start_idx, actual_block_size).noalias() = data->G.transpose() * Omg;
                H.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
            }
        }
    }
}

void FancyRsvdOpData::computeGandH(MyMatrix& G, MyMatrix& H, int pi)
{
    // check size of G and H first;
    if (H.cols() != size || H.rows() != cols() || G.cols() != size || G.rows() != rows())
    {
        throw std::runtime_error("Error: the size of G or H doesn't match.\n");
    }
    MyMatrix H1, H2;
    if (pi == 0)
    {
        auto rng = std::default_random_engine{};
        if (data->params.rand)
            Omg = StandardNormalRandom<MyMatrix, std::default_random_engine>(data->nsamples, size, rng);
        else
            Omg = UniformRandom<MyMatrix, std::default_random_engine>(data->nsamples, size, rng);
        Omg2 = Omg;
    }
    if (data->params.batch)
    {
        if (pi == 0)
        {
            data->llog << timestamp() << "running in memory mode with PCAone (algorithm2)." << endl;
            if (update)
            {
                data->update_batch_E(U, S, V.transpose());
            }
            if (standardize)
            {
                if (data->params.pcangsd)
                {
                    data->pcangsd_standardize_E(U, S, V.transpose());
                }
                else
                {
                    data->standardize_E();
                }
            }
            band = 1;
            blocksize = ceil((double)data->nsnps / data->params.bands);
            if (blocksize < data->params.bands)
                throw std::runtime_error("data and blocksize is too small. please consider IRAM method by using -a option");
            // permute snps of G, see https://stackoverflow.com/questions/15858569/randomly-permute-rows-columns-of-a-matrix-with-eigen
            if (!data->params.noshuffle)
            {
                perm = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(data->G.cols());
                perm.setIdentity();
                auto rng = std::default_random_engine{};
                std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(), rng);
                data->G = data->G * perm; // permute columns in-place
            }
        }
        {
            H1 = MyMatrix::Zero(cols(), size);
            H2 = MyMatrix::Zero(cols(), size);
            // band : 2, 8, 16, 32, 64, 128
            band = fmin(band * 2, data->params.bands);
            for (int b = 0, i = 1, j = 0; b < data->params.bands; ++b, ++i, ++j)
            {
                start_idx = b * blocksize;
                stop_idx = (b + 1) * blocksize >= data->nsnps ? data->nsnps - 1 : (b + 1) * blocksize - 1;
                actual_block_size = stop_idx - start_idx + 1;
                G.middleRows(start_idx, actual_block_size).noalias() = data->G.middleCols(start_idx, actual_block_size).transpose() * Omg;
                if (pi > 0 && j <= std::pow(2, pi - 1) && std::pow(2, pi) < data->params.bands)
                {
                    H1.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                    // additional complementary power iteration for last read
                    if (j == std::pow(2, pi - 1))
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        // start over reading the data for next power iteration using new Omg
                        b = -1;
                        i = 0;
                    }
                }
                else if (i <= band / 2)
                {
                    H1.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                }
                else if (i > band / 2 && i <= band)
                {
                    H2.noalias() += data->G.middleCols(start_idx, actual_block_size) * G.middleRows(start_idx, actual_block_size);
                }
                if ((b + 1) >= band)
                {
                    if (i == band)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H1.setZero();
                        i = 0;
                    }
                    else if (i == band / 2)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H2.setZero();
                    }
                    else if ((b + 1) == data->nblocks)
                    {
                        // shouldn't go here if the bands is proper, ie. 2^{x}
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                    }
                }
            }
        }
    }
    else
    {
        if (pi == 0)
        {
            if (verbose)
                data->llog << timestamp() << "running in out-of-core mode with PCAone (algorithm 2)." << endl;
            band = data->bandFactor;
        }
        {
            H1 = MyMatrix::Zero(cols(), size);
            H2 = MyMatrix::Zero(cols(), size);
            data->check_file_offset_first_var();
            // band : 2, 4, 8, 16, 32, 64
            band = fmin(band * 2, data->nblocks);
            for (int b = 0, i = 1, j = 0; b < data->nblocks; ++b, ++i, ++j)
            {
                start_idx = data->start[b];
                stop_idx = data->stop[b];
                actual_block_size = stop_idx - start_idx + 1;
                auto t1 = std::chrono::steady_clock::now();
                if (update)
                {
                    data->read_block_update(start_idx, stop_idx, U, S, V.transpose(), standardize);
                }
                else
                {
                    data->read_block_initial(start_idx, stop_idx, standardize);
                }
                auto t2 = std::chrono::steady_clock::now();
                data->readtime +=
                    std::chrono::duration<double>(t2 - t1).count() * std::chrono::duration<double>::period::num / std::chrono::duration<double>::period::den;
                G.middleRows(start_idx, actual_block_size).noalias() = data->G.transpose() * Omg;
                if (pi > 0 && j <= std::pow(2, pi - 1) * data->bandFactor && std::pow(2, pi) < data->params.bands)
                {
                    H1.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
                    // additional complementary power iteration for last read
                    if (j == std::pow(2, pi - 1))
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        // start over reading the data for next power iteration using new Omg
                        data->check_file_offset_first_var();
                        b = -1;
                        i = 0;
                    }
                }
                else if (i <= band / 2)
                {
                    H1.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
                }
                else if (i > band / 2 && i <= band)
                {
                    H2.noalias() += data->G * G.middleRows(start_idx, actual_block_size);
                }
                if ((b + 1) >= band)
                {
                    if (i == band)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H1 = MyMatrix::Zero(cols(), size);
                        i = 0;
                    }
                    else if (i == band / 2)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                        H2 = MyMatrix::Zero(cols(), size);
                    }
                    else if ((b + 1) == data->nblocks)
                    {
                        H = H1 + H2;
                        Eigen::HouseholderQR<MyMatrix> qr(H);
                        Omg.noalias() = qr.householderQ() * MyMatrix::Identity(cols(), size);
                        flip_Omg(Omg2, Omg);
                    }
                }
            }
        }
    }
}


void run_pca_with_halko(Data* data, const Param& params)
{
    if (params.batch)
    {
        data->llog << timestamp() << "begin to run PCAone RSVD in memory mode" << endl;
    }
    else
    {
        data->llog << timestamp() << "begin to PCAone RSVD in out-of-core mode" << endl;
    }
    MyMatrix Vpre;
    MyVector S;
    RsvdOpData* rsvd;
    if (!params.runem && params.fast)
    {
        rsvd = new FancyRsvdOpData(data, params.k, params.oversamples);
    }
    else
    {
        rsvd = new NormalRsvdOpData(data, params.k, params.oversamples);
    }
    if (!params.runem)
    {
        data->llog << timestamp() << "begin to do non-EM PCA." << endl;
        if (params.intype == FileType::CSV)
        {
            rsvd->setFlags(false, false, params.verbose);
        }
        else
        {
            rsvd->setFlags(false, true, params.verbose);
        }
        rsvd->computeUSV(params.maxp, params.tol);
        if (params.fast && params.batch)
        {
            // recover original order for V
            data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->U, rsvd->perm * rsvd->V);
        }
        else
        {
            data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->U, rsvd->V);
        }
    }
    else
    {
        // for EM iteration
        rsvd->setFlags(false, false, false);
        rsvd->computeUSV(params.maxp, params.tol);
        // flip_UV(rsvd->U, rsvd->V, false);
        double diff;
        rsvd->setFlags(true, false, false);
        data->llog << timestamp() << "begin to do EM PCA." << endl;
        for (uint i = 0; i < params.maxiter; ++i)
        {
            Vpre = rsvd->V;
            rsvd->computeUSV(params.maxp, params.tol);
            // flip_UV(rsvd->U, rsvd->V, false);
            if (params.mev)
                diff = 1.0 - mev(rsvd->V, Vpre);
            else
                diff = minSSE(rsvd->V, Vpre).sum() / Vpre.cols();
            data->llog << timestamp() << "individual allele frequencies estimated (iter=" << i + 1 << "), diff=" << diff << endl;
            if (diff < params.tolem)
            {
                data->llog << timestamp() << "come to convergence!" << endl;
                break;
            }
        }
        data->llog << timestamp() << "begin to standardize the matrix." << endl;

        // if pcangsd, estimate GRM.
        if (params.pcangsd)
        {
            data->pcangsd_standardize_E(rsvd->U, rsvd->S, rsvd->V.transpose());
            MyMatrix C = data->G * data->G.transpose();
            C.array() /= (double)data->nsnps;
            C.diagonal() = data->Dc.array() / (double)data->nsnps;
            std::ofstream out_cov(params.outfile + ".cov");
            if (out_cov.is_open())
            {
                out_cov << C << "\n";
            }
            // use Eigen::JacobiSVD to get eigenvecs
            Eigen::JacobiSVD<MyMatrix> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
            data->write_eigs_files(svd.singularValues().head(params.k), svd.matrixU().leftCols(params.k), svd.matrixU().leftCols(params.k));
        }
        else
        {
            rsvd->setFlags(true, true, false);
            rsvd->computeUSV(params.maxp, params.tol);
            data->write_eigs_files(rsvd->S.array().square() / data->nsnps, rsvd->U, rsvd->V);
        }
    }

    delete rsvd;
}
