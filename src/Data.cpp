#include "Data.hpp"

void Data::prepare()
{
    get_matrix_dimensions();

    if (params.batch)
    {
        read_all_and_centering();
    } else {
        if (params.blocksize > 0) {
            nblocks = (unsigned int)ceil((double)nsnps / params.blocksize);
            if (nblocks == 1) {
                cerr << "Warning: only one block exists. please use --batch mode instead.\n";
                exit(EXIT_FAILURE);
            }
            start.resize(nblocks);
            stop.resize(nblocks);
            for(uint i = 0 ; i < nblocks ; i++)
            {
                start[i] = i * params.blocksize;
                stop[i] = start[i] + params.blocksize - 1;
                stop[i] = stop[i] >= nsnps ? nsnps - 1 : stop[i];
            }
        } else {
            cerr << "ERROR: --blocksize must be greater than 0.\n";
            exit(EXIT_FAILURE);
        }
    }

}


/**
  T = U'/s
  VT = (U'/s) * G = T * G
  V = G' * (U/s) // calculate V is not a good idea
 **/

void Data::calcu_vt_initial(const MatrixXf& T, MatrixXf& VT)
{
    if (nblocks == 1) {
        cerr << "Warning: only one block exists. please use --batch mode instead.\n";
        exit(EXIT_SUCCESS);
    }
    uint actual_block_size;
    open_check_file();
    for(uint i = 0 ; i < nblocks ; ++i)
    {
        actual_block_size = stop[i] - start[i] + 1;
        // G (nsamples, actual_block_size)
        read_snp_block_initial(start[i], stop[i]);
        VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
    }
    close_check_file();

    return;
}

void Data::calcu_vt_update(const MatrixXf& T, const MatrixXf& U, const VectorXf& svals, MatrixXf& VT, bool standardize)
{
    if (nblocks == 1) {
        cerr << "Warning: only one block exists. please use --batch mode instead.\n";
        exit(EXIT_SUCCESS);
    }
    uint actual_block_size;
    open_check_file();
    for(uint i = 0 ; i < nblocks ; ++i)
    {
        actual_block_size = stop[i] - start[i] + 1;
        // G (nsamples, actual_block_size)
        read_snp_block_update(start[i], stop[i], U, svals, VT, standardize);
        VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
    }
    close_check_file();

    return;
}

// MatrixXf Data::calcu_vt_from_Eb(const MatrixXf& T, const MatrixXf& U, bool standardize )
// {
//     if (nblocks == 1) {
//         cerr << "Warning: only one block exists. please use --batch mode instead.\n";
//         exit(EXIT_SUCCESS);
//     }
//     uint nrow = T.rows();
//     MatrixXf VT(nrow, nsnps);
//     uint actual_block_size;
//     open_check_file();
//     for(uint i = 0 ; i < nblocks ; ++i)
//     {
//         actual_block_size = stop[i] - start[i] + 1;
//         // G (nsamples, actual_block_size)
//         update_block_E(start[i], stop[i], U, standardize);
//         VT.block(0, start[i], nrow, actual_block_size) = T * G.leftCols(actual_block_size);
//     }
//     close_check_file();

//     return VT;
// }

void Data::write_eigs_files(const VectorXf& vals, const MatrixXf& vecs)
{
    std::ofstream out_vals(params.outfile + ".eigvals");
    std::ofstream out_vecs(params.outfile + ".eigvecs");
    if (out_vals.is_open()) {
        out_vals << vals << '\n';
    }
    if (out_vecs.is_open()) {
        out_vecs << vecs << '\n';
    }

}

void Data::update_batch_E(const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT)
{
    uint ks = svals.size();
    if (params.pcangsd)
    {
        // for gp
        #pragma omp parallel for
        for (uint j = 0; j < nsnps; ++j) {
            double p0, p1, p2;
            for (uint i = 0; i < nsamples; ++i) {
                // Rescale individual allele frequencies
                double pt = 0.0;
                for (uint k = 0; k < ks; ++k) {
                    pt += U(i, k) * svals(k) * VT(k, j);
                }
                pt = (pt + 2.0 * F(j)) / 2.0;
                pt = fmin(fmax(pt, -F(j)), 1 - F(j));
                // update E, which is G here
                p0 = P(j, 3 * i + 0) * (1.0 - pt) * (1.0 - pt);
                p1 = P(j, 3 * i + 1) * 2 * pt * (1.0 - pt);
                p2 = P(j, 3 * i + 2) * pt * pt;
                G(i, j) = (p1 + 2*p2)/(p0 + p1 + p2) - 2.0 * F(j);
            }
        }
    } else {
        // for gt
        #pragma omp parallel for
        for (uint i = 0; i < nsnps; ++i) {
            for (uint j = 0; j < nsamples; ++j) {
                if (C[i * nsamples + j] & 1) { // sites need to be predicted
                    G(j, i) = 0.0;
                    for (uint k = 0; k < ks; ++k) {
                        G(j, i) += U(j, k) * svals(k) * VT(k, i);
                    }
                    G(j, i) = fmin(fmax(G(j, i), -F(i)), 1 - F(i));
                }
            }
        }
    }
}

void Data::standardize_E()
{
    #pragma omp parallel for
    for (uint i = 0; i < nsnps; ++i) {
        for (uint j = 0; j < nsamples; ++j) {
            G(j, i) /= sqrt(F(i) * (1 - F(i)));
        }
    }
}

void Data::pcangsd_standardize_E(const MatrixXf& U, const VectorXf& svals, const MatrixXf& VT)
{
    uint ks = svals.size();
    Dc = VectorXf::Zero(nsamples);
    #pragma omp parallel
    {
        VectorXf diag_private = VectorXf::Zero(nsamples); // Thread private vector;
        #pragma omp for
        for (uint j = 0; j < nsnps; j++) {
            double p0, p1, p2, pt, pSum, tmp;
            double norm = sqrt(2.0*F(j)*(1.0 - F(j)));
            for (uint i = 0; i < nsamples; i++) {
                // Rescale individual allele frequencies
                pt = 0.0;
                for (uint k = 0; k < ks; ++k) {
                    pt += U(i, k) * svals(k) * VT(k, j);
                }
                pt = (pt + 2.0 * F(j)) / 2.0;
                pt = fmin(fmax(pt, -F(j)), 1 - F(j));
                // Update e
                p0 = P(j, 3 * i + 0) * (1.0 - pt) * (1.0 - pt);
                p1 = P(j, 3 * i + 1) * 2 * pt * (1.0 - pt);
                p2 = P(j, 3 * i + 2) * pt * pt;
                pSum = p0 + p1 + p2;
                G(i, j) = (p1 + 2*p2)/pSum - 2.0 * F(j);
                G(i, j) = G(i, j) / norm;

                // Update diag
                tmp = (0.0 - 2.0 * F(j)) * (0.0 - 2.0 * F(j)) * (p0 / pSum);
                tmp = tmp + (1.0 - 2.0 * F(j)) * (1.0 - 2.0 * F(j)) * (p1 / pSum);
                tmp = tmp + (2.0 - 2.0 * F(j)) * (2.0 - 2.0 * F(j)) * (p2 /pSum);
                diag_private[i] += tmp/(2.0*F(j)*(1.0 - F(j)));
            }
        }
        #pragma omp critical
        {
            for (uint i = 0; i < nsamples; i++) {
                Dc[i] += diag_private[i]; // Sum arrays for threads
            }
        }
    }
}

/****

MatrixXf Data::calcu_block_matmul(const MatrixXf& X, bool rightside)
{
    uint n, m;
    uint actual_block_size;
    if (rightside)
    {   // Y = X * G
        // Y.block(0, b_start_idx, n, bs) = X.block(0, 0, n, nsamples) * G.block(0, 0, nsamples, bs)
        n = X.rows();
        m = nsnps;
        MatrixXf Y = MatrixXf::Zero(n, m);

        open_check_file();
        for(uint i = 0 ; i < nblocks ; ++i)
        {
            actual_block_size = stop[i] - start[i] + 1;
            // G (nsamples, actual_block_size)
            read_snp_block_initial(start[i], stop[i]);
            Y.block(0, start[i], n, actual_block_size).noalias() = X * G.leftCols(actual_block_size);
        }
        close_check_file();

        return Y;

    } else {
        // Y = G * X
        // bs = blocksize
        // Y += G.block(0, 0, nsamples, bs) * X.block(b_start_idx, 0, bs, X.cols())
        n = nsamples;
        m = X.cols();
        MatrixXf Y = MatrixXf::Zero(n, m);

        open_check_file();
        for(uint i = 0 ; i < nblocks ; ++i)
        {
            actual_block_size = stop[i] - start[i] + 1;
            // G (nsamples, actual_block_size)
            read_snp_block_initial(start[i], stop[i]);
            Y.noalias() = Y + G.leftCols(actual_block_size) * X.block(start[i], 0, actual_block_size, m);
        }
        close_check_file();

        return Y;
    }
}

MatrixXf Data::calcu_block_matmul(const MatrixXf& X, bool rightside, const MatrixXf& U, const VectorXf& S, const MatrixXf& V, bool standardize)
{
    uint n, m;
    uint actual_block_size;
    if (rightside)
    {   // Y = X * G
        // Y.block(0, b_start_idx, n, bs) = X.block(0, 0, n, nsamples) * G.block(0, 0, nsamples, bs)
        n = X.rows();
        m = nsnps;
        MatrixXf Y = MatrixXf::Zero(n, m);

        open_check_file();
        for(uint i = 0 ; i < nblocks ; ++i)
        {
            actual_block_size = stop[i] - start[i] + 1;
            // G (nsamples, actual_block_size)
            read_snp_block_update(start[i], stop[i], U, S, V.transpose(), standardize);
            Y.block(0, start[i], n, actual_block_size).noalias() = X * G.leftCols(actual_block_size);
        }
        close_check_file();

        return Y;

    } else {
        // Y = G * X
        // bs = blocksize
        // Y += G.block(0, 0, nsamples, bs) * X.block(b_start_idx, 0, bs, X.cols())
        n = nsamples;
        m = X.cols();
        MatrixXf Y = MatrixXf::Zero(n, m);

        open_check_file();
        for(uint i = 0 ; i < nblocks ; ++i)
        {
            actual_block_size = stop[i] - start[i] + 1;
            // G (nsamples, actual_block_size)
            read_snp_block_update(start[i], stop[i], U, S, V.transpose(), standardize);
            Y.noalias() = Y + G.leftCols(actual_block_size) * X.block(start[i], 0, actual_block_size, m);
        }
        close_check_file();

        return Y;
    }
}

MatrixXf Data::calcu_block_matmul_trans(const MatrixXf& X, bool rightside)
{
    uint n, m;
    uint actual_block_size;
    if (rightside)
    {
        cerr << "shouldn't come to here.\n";
        exit(EXIT_FAILURE);
    } else {
        // Y = G' * X
        // G is nsamples x nsnps,
        //let G = {a1, a2,...}, then G' = {a1', a2',...}
        n = nsnps;
        m = X.cols();
        MatrixXf Y = MatrixXf::Zero(n, m);

        open_check_file();
        for(uint i = 0 ; i < nblocks ; ++i)
        {
            actual_block_size = stop[i] - start[i] + 1;
            // G (nsamples, actual_block_size)
            read_snp_block_initial(start[i], stop[i]);
            Y.block(start[i], 0, actual_block_size, m).noalias() = G.leftCols(actual_block_size).transpose() * X ;
        }
        close_check_file();

        return Y;
    }
}

MatrixXf Data::calcu_block_matmul_trans(const MatrixXf& X, bool rightside, const MatrixXf& U, const VectorXf& S, const MatrixXf& V, bool standardize)
{
    uint n, m, actual_block_size;
    if (rightside)
    {
        cerr << "shouldn't come to here.\n";
        exit(EXIT_FAILURE);
    } else {
        // Y = G' * X
        // G is nsamples x nsnps,
        //let G = {a1, a2,...}, then G' = {a1', a2',...}
        n = nsnps;
        m = X.cols();
        MatrixXf Y = MatrixXf::Zero(n, m);

        open_check_file();
        for(uint i = 0 ; i < nblocks ; ++i)
        {
            actual_block_size = stop[i] - start[i] + 1;
            // G (nsamples, actual_block_size)
            read_snp_block_update(start[i], stop[i], U, S, V.transpose(), standardize);
            Y.block(start[i], 0, actual_block_size, m).noalias() = G.leftCols(actual_block_size).transpose() * X ;
        }
        close_check_file();

        return Y;
    }
}

// void Data::update_block_E(uint start_idx, uint stop_idx, const MatrixXf& U, bool standardize)
// {
//     if (params.intype == "bfile")
//     {
//         uint actual_block_size = stop_idx - start_idx + 1;
//         // if G is not initial then initial it
//         // if actual_block_size is smaller than blocksize, don't resize G;
//         if (G.cols() < params.blocksize || (actual_block_size < params.blocksize))
//         {
//             G = MatrixXf::Zero(nsamples, actual_block_size);
//         }
//         Cb = MatrixXf::Zero(nsamples, actual_block_size);
//         // check where we are
//         long long offset = 3 + start_idx * bed_bytes_per_snp;
//         if (bed_ifstream.tellg() != offset)
//         {
//             cerr << "Error: something wrong with read_snp_block!\n";
//             exit(EXIT_FAILURE);
//         }
//         uint bi, ki, i, j;
//         inbed.resize(bed_bytes_per_snp * actual_block_size);
//         bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp * actual_block_size);
//         #pragma omp parallel for private(i,j,bi,ki)
//         for (i = 0; i < actual_block_size; ++i)
//         {
//             uint snp_idx = start_idx + i;
//             for (bi = 0, j = 0; bi < bed_bytes_per_snp; ++bi)
//             {
//                 uchar buf = inbed[i * bed_bytes_per_snp + bi];
//                 for (ki=0; ki<4; ++ki, ++j)
//                 {
//                     if (j < nsamples)
//                     {
//                         G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
//                         if ((buf & 3) != 1)
//                         {
//                             Cb(j, i) = 1; // not missing, weight is 1
//                         }
//                         buf = buf >> 2;  // shift packed data and throw away genotype just processed.
//                     }
//                 }
//             }
//         }

//         // get Beta matrix and update Eb

//         uint k = U.cols();
//         Bb = MatrixXf::Zero(k, actual_block_size);
//         #pragma omp parallel for
//         for (i = 0; i < actual_block_size; ++i)
//         {
//             // WXb = Wy
//             Bb.col(i) = (Cb.col(i).asDiagonal() * U).jacobiSvd(ComputeThinU | ComputeThinV).solve(Cb.col(i).asDiagonal() * G.col(i));
//             uint snp_idx = start_idx + i;
//             for(j = 0; j < nsamples; ++j)
//             {
//                 if (Cb(j, i) == 0 && G(j, i) == 0)
//                 {
//                     for(ki = 0; ki < k; ++ki)
//                     {
//                         // update
//                         G(j, i) = G(j, i) + U(j, ki) * Bb(ki, i);
//                         G(j, i) = fmin(fmax(G(j, i), -F(snp_idx)), 1 - F(snp_idx));
//                         if (standardize) G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
//                     }
//                 }
//             }
//         }

//     } else {
//         cerr << "pfile and bgen mode is coming.\n";
//     }
// }

****/

