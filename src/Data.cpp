#include "Data.hpp"

using namespace std;

Logger::Logger()
{
}
Logger::~Logger()
{
}

void Data::prepare(uint& blocksize)
{
    if (nsamples > nsnps)
        nsamples_ge_nsnps = true;

    if (params.batch)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        read_all_and_centering();
        auto t2 = std::chrono::high_resolution_clock::now();
        readtime += std::chrono::duration<double>(t2 - t1).count() * std::chrono::duration<double>::period::num / std::chrono::duration<double>::period::den;
    }
    else
    {
        // some common settings
        if (params.arnoldi)
        {
            // ram of arnoldi = n * b * 8 / 1024 kb
            blocksize = (unsigned int)ceil((double)params.memory * 134217728 / nsamples);
        }
        else
        {
            // ram of halko = (3*n*l + 2*m*l + 5*m + n*b)*8/1024 Kb
            uint64 l = params.k + params.oversamples;
            double m = (double)(3 * nsamples * l + 2 * nsnps * l + 5 * nsnps) / 134217728;
            if (params.memory > 1.1 * m)
            {
                m = 0;
            }
            else
            {
                llog << "Waring: minimum RAM required is " << m << " GB. Trying to allocate more RAM. " << endl;
            }
            blocksize = (unsigned int)ceil((double)((m + params.memory) * 134217728 - 3 * nsamples * l - 2 * nsnps * l - 5 * nsnps) / nsamples);
        }
        nblocks = (unsigned int)ceil((double)nsnps / blocksize);
        if (params.verbose)
            llog << timestamp() << "initial setting by -m/--memory: blocksize=" << blocksize << ", nblocks=" << nblocks << ", bandFactor=" << bandFactor
                 << ".\n";
        if (nblocks == 1)
        {
            throw std::invalid_argument("Warning: only one block exists. please remove -m / --memory option instead.\n");
        }
        if (params.fast)
        {
            // decrease blocksize to fit the fancy halko
            if (nblocks < params.bands)
            {
                blocksize = (unsigned int)ceil((double)nsnps / params.bands);
            }
            else
            {
                bandFactor = (unsigned int)ceil((double)nblocks / params.bands);
                blocksize = (unsigned int)ceil((double)nsnps / (params.bands * bandFactor));
            }
            nblocks = (unsigned int)ceil((double)nsnps / blocksize);
            if (params.verbose)
                llog << timestamp() << "after adjustment by -f/--fast: blocksize=" << blocksize << ", nblocks=" << nblocks << ", bandFactor=" << bandFactor
                     << ".\n";
        }
        start.resize(nblocks);
        stop.resize(nblocks);
        for (uint i = 0; i < nblocks; i++)
        {
            start[i] = i * blocksize;
            stop[i] = start[i] + blocksize - 1;
            stop[i] = stop[i] >= nsnps ? nsnps - 1 : stop[i];
        }
        // initial some variables for blockwise for specific files here.
        if (params.intype != FileType::CSV)
            F = MyVector::Zero(nsnps);
        if (params.intype == FileType::PLINK)
            centered_geno_lookup = MyArrayX::Zero(4, nsnps); // for plink input
    }
}


/**
  T = U'/s
  VT = (U'/s) * G = T * G
  V = G' * (U/s) // calculate V is not a good idea
 **/

void Data::calcu_vt_initial(const MyMatrix& T, MyMatrix& VT)
{
    if (nblocks == 1)
    {
        llog << "Warning: only one block exists. please use --batch mode instead.\n";
        exit(EXIT_SUCCESS);
    }
    uint actual_block_size;
    check_file_offset_first_var();
    for (uint i = 0; i < nblocks; ++i)
    {
        actual_block_size = stop[i] - start[i] + 1;
        // G (nsamples, actual_block_size)
        read_snp_block_initial(start[i], stop[i]);
        VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
    }

    return;
}

void Data::calcu_vt_update(const MyMatrix& T, const MyMatrix& U, const MyVector& svals, MyMatrix& VT, bool standardize)
{
    if (nblocks == 1)
    {
        llog << "Warning: only one block exists. please use --batch mode instead.\n";
        exit(EXIT_SUCCESS);
    }
    uint actual_block_size;
    check_file_offset_first_var();
    for (uint i = 0; i < nblocks; ++i)
    {
        actual_block_size = stop[i] - start[i] + 1;
        // G (nsamples, actual_block_size)
        read_snp_block_update(start[i], stop[i], U, svals, VT, standardize);
        VT.block(0, start[i], T.rows(), actual_block_size) = T * G.leftCols(actual_block_size);
    }

    return;
}

// MyMatrix Data::calcu_vt_from_Eb(const MyMatrix& T, const MyMatrix& U, bool standardize )
// {
//     if (nblocks == 1) {
//         llog << "Warning: only one block exists. please use --batch mode instead.\n";
//         exit(EXIT_SUCCESS);
//     }
//     uint nrow = T.rows();
//     MyMatrix VT(nrow, nsnps);
//     uint actual_block_size;
//     check_file_offset_first_var();
//     for(uint i = 0 ; i < nblocks ; ++i)
//     {
//         actual_block_size = stop[i] - start[i] + 1;
//         // G (nsamples, actual_block_size)
//         update_block_E(start[i], stop[i], U, standardize);
//         VT.block(0, start[i], nrow, actual_block_size) = T * G.leftCols(actual_block_size);
//     }

//     return VT;
// }

void Data::write_eigs_files(const MyVector& S, const MyMatrix& U, const MyMatrix& V)
{
    std::ofstream outs(params.outfile + ".eigvals");
    std::ofstream outu(params.outfile + ".eigvecs");
    if (outs.is_open())
    {
        outs << S << '\n';
    }
    if (outu.is_open())
    {
        outu << U << '\n';
    }
    if (params.printv)
    {
        std::ofstream outv(params.outfile + ".loadings");
        if (outv.is_open())
        {
            outv << V << '\n';
        }
    }
}

void Data::update_batch_E(const MyMatrix& U, const MyVector& svals, const MyMatrix& VT)
{
    uint ks = svals.size();
    if (params.pcangsd)
    {
// for gp
#pragma omp parallel for
        for (uint j = 0; j < nsnps; ++j)
        {
            double p0, p1, p2;
            for (uint i = 0; i < nsamples; ++i)
            {
                // Rescale individual allele frequencies
                double pt = 0.0;
                for (uint k = 0; k < ks; ++k)
                {
                    pt += U(i, k) * svals(k) * VT(k, j);
                }
                pt = (pt + 2.0 * F(j)) / 2.0;
                pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
                // update E, which is G here
                p0 = P(3 * i + 0, j) * (1.0 - pt) * (1.0 - pt);
                p1 = P(3 * i + 1, j) * 2 * pt * (1.0 - pt);
                p2 = P(3 * i + 2, j) * pt * pt;
                G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
            }
        }
    }
    else
    {
// for gt
#pragma omp parallel for
        for (uint i = 0; i < nsnps; ++i)
        {
            for (uint j = 0; j < nsamples; ++j)
            {
                if (C[i * nsamples + j] & 1)
                { // sites need to be predicted
                    G(j, i) = 0.0;
                    for (uint k = 0; k < ks; ++k)
                    {
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
    for (uint i = 0; i < nsnps; ++i)
    {
        for (uint j = 0; j < nsamples; ++j)
        {
            // in case denominator is too small.
            if (sqrt(F(i) * (1 - F(i))) > VAR_TOL)
                G(j, i) /= sqrt(F(i) * (1 - F(i)));
        }
    }
}

void Data::pcangsd_standardize_E(const MyMatrix& U, const MyVector& svals, const MyMatrix& VT)
{
    uint ks = svals.size();
    Dc = MyVector::Zero(nsamples);
#pragma omp parallel
    {
        MyVector diag_private = MyVector::Zero(nsamples); // Thread private vector;
#pragma omp for
        for (uint j = 0; j < nsnps; j++)
        {
            double p0, p1, p2, pt, pSum, tmp;
            double norm = sqrt(2.0 * F(j) * (1.0 - F(j)));
            for (uint i = 0; i < nsamples; i++)
            {
                // Rescale individual allele frequencies
                pt = 0.0;
                for (uint k = 0; k < ks; ++k)
                {
                    pt += U(i, k) * svals(k) * VT(k, j);
                }
                pt = (pt + 2.0 * F(j)) / 2.0;
                pt = fmin(fmax(pt, 1e-4), 1.0 - 1e-4);
                // Update e
                p0 = P(3 * i + 0, j) * (1.0 - pt) * (1.0 - pt);
                p1 = P(3 * i + 1, j) * 2 * pt * (1.0 - pt);
                p2 = P(3 * i + 2, j) * pt * pt;
                pSum = p0 + p1 + p2;
                G(i, j) = (p1 + 2 * p2) / pSum - 2.0 * F(j);
                if (norm > VAR_TOL)
                    G(i, j) /= norm;

                // Update diag
                tmp = (0.0 - 2.0 * F(j)) * (0.0 - 2.0 * F(j)) * (p0 / pSum);
                tmp += (1.0 - 2.0 * F(j)) * (1.0 - 2.0 * F(j)) * (p1 / pSum);
                tmp += (2.0 - 2.0 * F(j)) * (2.0 - 2.0 * F(j)) * (p2 / pSum);
                diag_private[i] += tmp / (2.0 * F(j) * (1.0 - F(j)));
            }
        }
#pragma omp critical
        {
            for (uint i = 0; i < nsamples; i++)
            {
                Dc[i] += diag_private[i]; // Sum arrays for threads
            }
        }
    }
}
