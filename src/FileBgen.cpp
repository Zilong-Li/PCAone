#include "FileBgen.hpp"

using namespace std;

void FileBgen::read_all()
{
    uint i, j, k, gc;
    double gs, af;
    if(!params.pcangsd)
    {
        F = MyVector::Zero(nsnps);
        G = MyMatrix::Zero(nsamples, nsnps);
        if(params.runem) C = ArrayXb::Zero(nsnps * nsamples);
        for(j = 0, k = 0; j < nsnps; j++)
        {
            try
            {
                var = bg->next_var();
                dosages = var.minor_allele_dosage();
                gc = 0;
                gs = 0.0;
// calculate allele frequency
#pragma omp parallel for reduction(+ : gc) reduction(+ : gs)
                for(i = 0; i < nsamples; i++)
                {
                    if(!std::isnan(dosages[i]))
                    {
                        gs += dosages[i] / 2.0; // map to [0, 1];
                        gc += 1;
                    }
                }
                if(gc == 0)
                    af = 0.0;
                else
                    af = (double)gs / gc;
                if(af > params.maf)
                    F(k) = af;
                else
                    continue;
// do centering and initialing
#pragma omp parallel for
                for(i = 0; i < nsamples; i++)
                {
                    if(std::isnan(dosages[i]))
                    {
                        if(params.runem) C[k * nsamples + i] = 1;
                        G(i, k) = 0;
                    }
                    else
                    {
                        if(params.runem) C[k * nsamples + i] = 0;
                        G(i, j) = dosages[i] / 2.0 - F(k); // map to [0, 1];
                    }
                }
                k++;
            }
            catch(const std::out_of_range & e)
            {
                throw e.what();
            }
        }
        if(k == 0)
            throw std::runtime_error("the number of SNPs after filtering should be 0!\n");
        else
            llog << timestamp() << "number of SNPs after filtering by MAF > " << params.maf << ": " << k
                 << endl;
        // resize G, F, C;
        nsnps = k; // resize nsnps;
        G.conservativeResize(Eigen::NoChange, nsnps);
        F.conservativeResize(nsnps);
        C.conservativeResize(nsnps * nsamples);
    }
    else
    {
        // read all GP data into P;
        for(j = 0; j < nsnps; j++)
        {
            try
            {
                var = bg->next_var();
                probs1d = var.probs_1d();
#pragma omp parallel for
                for(i = 0; i < nsamples; i++)
                {
                    P(i * 2 + 0, j) = probs1d[i * 3 + 0];
                    P(i * 2 + 1, j) = probs1d[i * 3 + 1];
                    // no need to parse probs1d[i * 3 + 2]
                }
            }
            catch(const std::out_of_range & e)
            {
                break;
            }
        }
        assert(j == nsnps);
        llog << timestamp() << "begin to estimate allele frequencies using GP" << endl;
        MyVector Ft(nsnps);
        F = MyVector::Constant(nsnps, 0.25);
        // run EM to estimate allele frequencies
        double diff;
        for(uint it = 0; it < params.maxiter; it++)
        {
#pragma omp parallel for
            for(uint j = 0; j < nsnps; j++)
            {
                Ft(j) = F(j);
                double p0, p1, p2, pt = 0.0;
                for(uint i = 0; i < nsamples; i++)
                {
                    p0 = P(2 * i + 0, j) * (1.0 - F(j)) * (1.0 - F(j));
                    p1 = P(2 * i + 1, j) * 2 * F(j) * (1.0 - F(j));
                    p2 = (1 - P(2 * i + 0, j) - P(2 * i + 1, j)) * F(j) * F(j);
                    pt += (p1 + 2 * p2) / (2 * (p0 + p1 + p2));
                }
                F(j) = pt / (double)nsamples;
            }
            // calculate differences between iterations
            diff = sqrt((F - Ft).array().square().sum() / nsnps);
            // Check for convergence
            if(diff < params.tolmaf)
            {
                llog << "EM (MAF) converged at iteration: " << it + 1 << endl;
                break;
            }
            else if(it == (params.maxiter - 1))
            {
                llog << "EM (MAF) did not converge.\n";
            }
        }
        filterSNPs_resizeF();
        // initial E which is G
        G = MyMatrix::Zero(nsamples, nsnps);
#pragma omp parallel for
        for(j = 0; j < nsnps; j++)
        {
            double p0, p1, p2;
            for(i = 0; i < nsamples; i++)
            {
                p0 = P(2 * i + 0, keepSNPs[j]) * (1.0 - F(j)) * (1.0 - F(j));
                p1 = P(2 * i + 1, keepSNPs[j]) * 2 * F(j) * (1.0 - F(j));
                p2 = (1 - P(2 * i + 0, keepSNPs[j]) - P(2 * i + 1, keepSNPs[j])) * F(j) * F(j);
                G(i, j) = (p1 + 2 * p2) / (p0 + p1 + p2) - 2.0 * F(j);
            }
        }
    }
}

void FileBgen::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    read_bgen_block(G, F, var, bg, dosages, probs1d, frequency_was_estimated, nsamples, nsnps,
                    params.blocksize, start_idx, stop_idx, standardize);
}

void read_bgen_block(MyMatrix & G,
                     MyVector & F,
                     bgen::Variant & var,
                     bgen::Bgen * bg,
                     float * dosages,
                     float * probs1d,
                     bool & frequency_was_estimated,
                     uint64 nsamples,
                     uint64 nsnps,
                     uint blocksize,
                     uint64 start_idx,
                     uint64 stop_idx,
                     bool standardize)
{
    uint actual_block_size = stop_idx - start_idx + 1;
    uint i, j, snp_idx;
    if(G.cols() < blocksize || (actual_block_size < blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
    }
    if(frequency_was_estimated)
    {
        for(i = 0; i < actual_block_size; ++i)
        {
            snp_idx = start_idx + i;
            var = bg->next_var();
            dosages = var.minor_allele_dosage();
#pragma omp parallel for
            for(j = 0; j < nsamples; j++)
            {
                if(std::isnan(dosages[j]))
                {
                    G(j, i) = 0;
                }
                else
                {
                    G(j, i) = dosages[j] / 2.0 - F(snp_idx);
                }
                if(standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
                    G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
            }
        }
    }
    else
    {
        uint gc;
        double gs;
        for(i = 0; i < actual_block_size; ++i)
        {
            snp_idx = start_idx + i;
            var = bg->next_var();
            dosages = var.minor_allele_dosage();
            gc = 0;
            gs = 0.0;
#pragma omp parallel for reduction(+ : gc) reduction(+ : gs)
            for(j = 0; j < nsamples; j++)
            {
                if(std::isnan(dosages[j]))
                {
                    G(j, i) = 0;
                }
                else
                {
                    G(j, i) = dosages[j] / 2.0;
                    gs += G(j, i);
                    gc += 1;
                }
            }
            if(gc == 0)
                throw std::runtime_error(
                    "Error: the allele frequency should not be 0. should do filtering first.");
            F(snp_idx) = (double)gs / gc;
// do centering
#pragma omp parallel for
            for(j = 0; j < nsamples; j++)
            {
                if(std::isnan(dosages[j]))
                {
                    G(j, i) = 0;
                }
                else
                {
                    G(j, i) -= F(snp_idx);
                }
                if(standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
                    G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
            }
        }
        if(stop_idx + 1 == nsnps) frequency_was_estimated = true;
    }
}

int shuffle_bgen_to_bin(std::string bgenfile, std::string binfile, uint gb, bool standardize)
{
    bgen::Bgen * bg = new bgen::Bgen(bgenfile, "", true);
    uint64 nsamples = bg->header.nsamples;
    uint64 nsnps = bg->header.nvariants;
    uint64 twoGB = (uint64)1073741824 * gb;
    uint64 blocksize = twoGB / (nsamples * sizeof(double));
    uint nblocks = (nsnps + blocksize - 1) / blocksize;
    std::ofstream ofs(binfile, std::ios::binary);
    ofs.write((char *)&nsamples, sizeof(nsamples));
    ofs.write((char *)&nsnps, sizeof(nsnps));
    bgen::Variant var;
    float * dosages = nullptr;
    float * probs1d = nullptr;
    bool frequency_was_estimated = false;
    std::vector<uint64> perm(nsnps);
    std::iota(perm.begin(), perm.end(), 0);
    auto rng = std::default_random_engine{};
    std::shuffle(perm.begin(), perm.end(), rng);
    MyMatrix G;
    MyVector F(nsnps);
    uint64 idx, cur = 0;
    uint64 bytes_per_snp = nsamples * sizeof(double);
    for(uint i = 0; i < nblocks; i++)
    {
        auto start_idx = i * blocksize;
        auto stop_idx = start_idx + blocksize - 1;
        stop_idx = stop_idx >= nsnps ? nsnps - 1 : stop_idx;
        read_bgen_block(G, F, var, bg, dosages, probs1d, frequency_was_estimated, nsamples, nsnps, blocksize,
                        start_idx, stop_idx, standardize);
        for(size_t p = 0; p < G.cols(); p++, cur++)
        {
            // std::cerr << cur << "," << perm[cur] << "\n";
            idx = 2 * sizeof(nsamples) + perm[cur] * bytes_per_snp;
            ofs.seekp(idx, std::ios_base::beg);
            ofs.write((char *)G.col(p).data(), bytes_per_snp);
        }
    }
    return (nsnps == cur);
}
