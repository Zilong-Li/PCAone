#include "FilePlink.hpp"

using namespace std;

// https://stackoverflow.com/questions/5166263/how-to-get-iostream-to-perform-better
// https://github.com/OpenGene/fastp/blame/master/src/fastareader.cpp#L11
void FileBed::check_file_offset_first_var()
{
    setlocale(LC_ALL,"C");
    ios_base::sync_with_stdio(false);
    long long offset = 3 + nsnps * bed_bytes_per_snp;
    if (bed_ifstream.tellg() == offset) {
        // reach the end of bed, reset the position to the first variant;
        bed_ifstream.seekg(3, std::ios_base::beg);
    } else if (bed_ifstream.tellg() == 3) {
        ;
    } else {
        cout << bed_ifstream.tellg() << ".\n";
        throw std::runtime_error("Warning: the bed_ifstream is pointing an unexpected position.\n");
    }
}

void FileBed::read_all_and_centering()
{
    F = MyVector::Zero(nsnps);
    check_file_offset_first_var();
    G = MyMatrix(nsamples, nsnps);
    // Begin to decode the plink bed
    inbed.reserve(bed_bytes_per_snp * nsnps); 
    bed_ifstream.read(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp * nsnps);
    if (params.runem) C.resize(nsnps * nsamples);
    uint64 c, i, j, b, k;
    uchar buf;
#pragma omp parallel for private(i,j,b,c,k,buf)
    for(i = 0; i < nsnps; ++i)
    {
        for (b=0, c=0, j=0; b < bed_bytes_per_snp; ++b)
        {
            buf = inbed[i * bed_bytes_per_snp + b];
            for (k=0; k<4; ++k, ++j)
            {
                if (j < nsamples)
                {
                    G(j, i) = BED2GENO[buf & 3];
                    if (G(j, i) != BED_MISSING_VALUE) {
                        // 0 indicate G(i,j) don't need to be predicted.
                        if (params.runem) C[i * nsamples + j] = 0;
                        F(i) += G(j, i);
                        c++;
                    } else {
                        // 1 indicate G(i,j) need to be predicted and updated.
                        if (params.runem) C[i * nsamples + j] = 1;
                    }
                    buf >>= 2;  // shift packed data and throw away genotype just processed.
                }
            }
        }
        if (c==0) throw std::runtime_error("Error: the allele frequency should not be 0. should do filtering first.");
        F(i) /= c;
        // do centering and initialing
        for(j=0; j<nsamples; ++j)
        {
            if (G(j, i) == BED_MISSING_VALUE) {
                G(j, i) = 0.0;
            } else {
                G(j, i) -= F(i);
            }
        }
    }
    // read bed to matrix G done and close bed_ifstream
    inbed.clear();
    inbed.shrink_to_fit();
}

void FileBed::read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    uint actual_block_size = stop_idx - start_idx + 1;
    // if G is not initial then initial it
    // if actual_block_size is smaller than blocksize, don't resize G;
    if (G.cols() < params.blocksize || (actual_block_size < params.blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
        inbed.reserve(bed_bytes_per_snp * params.blocksize);
    }
    // check where we are
    if (params.verbose) {
        long long offset = 3 + start_idx * bed_bytes_per_snp;
        if (bed_ifstream.tellg() != offset)
        {
            throw std::runtime_error("Error: something wrong with read_snp_block!\n");
        }
    }
    uint64 c, b, i, j, k, snp_idx;
    uchar buf;
    // inbed.resize(bed_bytes_per_snp * actual_block_size);
    bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp * actual_block_size);
    // bed_ifstream.rdbuf()->sgetn(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp * actual_block_size);
    if (frequency_was_estimated)
    {
#pragma omp parallel for private(i,j,b,k,snp_idx,buf)
        for (i = 0; i < actual_block_size; ++i)
        {
            snp_idx = start_idx + i;
            for (b = 0, j = 0; b < bed_bytes_per_snp; ++b)
            {
                buf = inbed[i * bed_bytes_per_snp + b];
                for (k=0; k<4; ++k, ++j)
                {
                    if (j < nsamples)
                    {
                        G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
                        if (standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL) G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
                        buf >>= 2;  // shift packed data and throw away genotype just processed.
                    }
                }
            }
        }
    } else {
        // estimate allele frequencies
#pragma omp parallel for private(c,i,j,b,k,snp_idx,buf)
        for (i = 0; i < actual_block_size; ++i)
        {
            snp_idx = start_idx + i;
            c = 0;
            for (b = 0, j = 0; b < bed_bytes_per_snp; ++b)
            {
                buf = inbed[i * bed_bytes_per_snp + b];
                for (k=0; k<4; ++k, ++j)
                {
                    if (j < nsamples)
                    {
                        if ((buf & 3) != 1) {
                            // g is {0, 0.5, 1}
                            F(snp_idx) += BED2GENO[buf & 3];
                            c++;
                        }
                        buf >>= 2;  // shift packed data and throw away genotype just processed.
                    // } else {
                    //     if (buf != 0)
                    //     {
                    //         cerr << "Row " << snp_idx + 1 << "padding is non-zero. Either the specified number of individuals is incorrect or the input file is corrupt!\n";
                    //         exit(EXIT_FAILURE);
                    //     }
                    }
                }
            }
            // calculate F and centered_geno_lookup
            if (c==0) throw std::runtime_error("Error: the allele frequency should not be 0. should do filtering first.");
            F(snp_idx) /= c;
            // do centering and initialing
            centered_geno_lookup(1, snp_idx) = 0.0; // missing
            centered_geno_lookup(0, snp_idx) = BED2GENO[0] - F(snp_idx); // minor hom
            centered_geno_lookup(2, snp_idx) = BED2GENO[2] - F(snp_idx); // het
            centered_geno_lookup(3, snp_idx) = BED2GENO[3] - F(snp_idx); // major hom
            // get centered and standardized G
            for (b = 0, j = 0; b < bed_bytes_per_snp; ++b)
            {
                buf = inbed[i * bed_bytes_per_snp + b];
                for (k=0; k<4; ++k, ++j)
                {
                    if (j < nsamples)
                    {
                        G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
                        if (standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL) G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
                        buf >>= 2;  // shift packed data and throw away genotype just processed.
                    }
                }
            }
        }
    }

    if (stop_idx + 1 == nsnps) frequency_was_estimated = true;

}

void FileBed::read_snp_block_update(uint64 start_idx, uint64 stop_idx, const MyMatrix& U, const MyVector& svals, const MyMatrix& VT, bool standardize)
{
    uint actual_block_size = stop_idx - start_idx + 1;
    if (G.cols() < params.blocksize || (actual_block_size < params.blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
        inbed.reserve(bed_bytes_per_snp * params.blocksize);
    }
    // check where we are
    if (params.verbose) {
        long long offset = 3 + start_idx * bed_bytes_per_snp;
        if (bed_ifstream.tellg() != offset)
        {
            throw std::runtime_error("Error: something wrong with read_snp_block!\n");
        }
    }
    uint64 b, i, j, snp_idx;
    uint ks = svals.rows();
    uint ki, k;
    uchar buf;
    bed_ifstream.read( reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp * actual_block_size);
#pragma omp parallel for private(i,j,b,ki,k,snp_idx,buf)
    for (i = 0; i < actual_block_size; ++i)
    {
        snp_idx = start_idx + i;
        for (b = 0, j = 0; b < bed_bytes_per_snp; ++b)
        {
            buf = inbed[i * bed_bytes_per_snp + b];
            for (ki=0; ki<4; ++ki, ++j)
            {
                if (j < nsamples)
                {
                    G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
                    if (BED2GENO[buf & 3] == BED_MISSING_VALUE)
                    {
                        G(j, i) = 0.0;
                        for (k = 0; k < ks; ++k)
                        {
                            G(j, i) += U(j, k) * svals(k) * VT(k, snp_idx);
                        }
                        // map to domain(0,1)
                        G(j, i) = fmin(fmax(G(j, i), -F(snp_idx)), 1 - F(snp_idx));
                    }
                    if (standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL) G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
                    buf >>= 2;  // shift packed data and throw away genotype just processed.
                }
            }
        }
    }

}
