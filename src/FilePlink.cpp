#include "FilePlink.hpp"

#include "Utils.hpp"

using namespace std;

void FileBed::check_file_offset_first_var()
{
    setlocale(LC_ALL, "C");
    ios_base::sync_with_stdio(false);
    long long offset = 3 + nsnps * bed_bytes_per_snp;
    if(bed_ifstream.tellg() == offset)
    {
        // reach the end of bed, reset the position to the first variant;
        bed_ifstream.seekg(3, std::ios_base::beg);
    }
    else if(bed_ifstream.tellg() == 3)
    {
        ;
    }
    else
    {
        bed_ifstream.seekg(3, std::ios_base::beg);
        if(params.verbose)
            std::cout << colwarn + "make sure you are runing PCAone algorithm2" + colend << std::endl;
    }
}

void FileBed::read_all()
{
    check_file_offset_first_var();
    // Begin to decode the plink bed
    inbed.reserve(bed_bytes_per_snp * nsnps);
    bed_ifstream.read(reinterpret_cast<char *>(&inbed[0]), bed_bytes_per_snp * nsnps);
    uint64 c, i, j, b, k;
    uchar buf;
    F = MyVector::Zero(nsnps);
    // estimate allele frequency first
#pragma omp parallel for private(i, j, b, c, k, buf)
    for(i = 0; i < nsnps; ++i)
    {
        for(b = 0, c = 0, j = 0; b < bed_bytes_per_snp; ++b)
        {
            buf = inbed[i * bed_bytes_per_snp + b];
            for(k = 0; k < 4; ++k, ++j)
            {
                if(j < nsamples)
                {
                    if(BED2GENO[buf & 3] != BED_MISSING_VALUE)
                    {
                        F(i) += BED2GENO[buf & 3];
                        c++;
                    }
                    buf >>= 2; // shift packed data and throw away genotype just processed.
                }
            }
        }
        if(c == 0)
            F(i) = 0;
        else
            F(i) /= c;
    }
    // filter and resize nsnps
    filterSNPs_resizeF();
    // fill in G with new size
    G = MyMatrix::Zero(nsamples, nsnps);
    if(params.runem) C = ArrayXb::Zero(nsnps * nsamples);
#pragma omp parallel for private(i, j, b, c, k, buf)
    for(i = 0; i < nsnps; ++i)
    {
        for(b = 0, c = 0, j = 0; b < bed_bytes_per_snp; ++b)
        {
            buf = inbed[keepSNPs[i] * bed_bytes_per_snp + b];
            for(k = 0; k < 4; ++k, ++j)
            {
                if(j < nsamples)
                {
                    G(j, i) = BED2GENO[buf & 3];
                    if(G(j, i) != BED_MISSING_VALUE)
                    {
                        // 0 indicate G(i,j) don't need to be predicted.
                        if(params.runem) C[i * nsamples + j] = 0;
                    }
                    else
                    {
                        // 1 indicate G(i,j) need to be predicted and updated.
                        if(params.runem) C[i * nsamples + j] = 1;
                    }
                    buf >>= 2; // shift packed data and throw away genotype just processed.
                }
            }
        }
        // do centering and initialing
        for(j = 0; j < nsamples; ++j)
        {
            if(G(j, i) == BED_MISSING_VALUE)
                G(j, i) = 0.0;
            else
                G(j, i) -= F(i);
        }
    }
    // read bed to matrix G done and close bed_ifstream
    inbed.clear();
    inbed.shrink_to_fit();
}

void FileBed::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    uint actual_block_size = stop_idx - start_idx + 1;
    // check where we are
    long long offset = 3 + start_idx * bed_bytes_per_snp;
    if(bed_ifstream.tellg() != offset)
        throw std::runtime_error("Error: something wrong with read_snp_block!\n");
    // if G is not initial then initial it
    // if actual_block_size is smaller than blocksize, don't resize G;
    if(G.cols() < params.blocksize || (actual_block_size < params.blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
        inbed.reserve(bed_bytes_per_snp * params.blocksize);
    }
    uint64 c, b, i, j, k, snp_idx;
    uchar buf;
    // inbed.resize(bed_bytes_per_snp * actual_block_size);
    bed_ifstream.read(reinterpret_cast<char *>(&inbed[0]), bed_bytes_per_snp * actual_block_size);
    // bed_ifstream.rdbuf()->sgetn(reinterpret_cast<char *> (&inbed[0]), bed_bytes_per_snp *
    // actual_block_size);
    if(frequency_was_estimated)
    {
#pragma omp parallel for private(i, j, b, k, snp_idx, buf)
        for(i = 0; i < actual_block_size; ++i)
        {
            snp_idx = start_idx + i;
            for(b = 0, j = 0; b < bed_bytes_per_snp; ++b)
            {
                buf = inbed[i * bed_bytes_per_snp + b];
                for(k = 0; k < 4; ++k, ++j)
                {
                    if(j < nsamples)
                    {
                        G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
                        if(standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
                            G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
                        buf >>= 2; // shift packed data and throw away genotype just processed.
                    }
                }
            }
        }
    }
    else
    {
        // estimate allele frequencies
#pragma omp parallel for private(c, i, j, b, k, snp_idx, buf)
        for(i = 0; i < actual_block_size; ++i)
        {
            snp_idx = start_idx + i;
            c = 0;
            for(b = 0, j = 0; b < bed_bytes_per_snp; ++b)
            {
                buf = inbed[i * bed_bytes_per_snp + b];
                for(k = 0; k < 4; ++k, ++j)
                {
                    if(j < nsamples)
                    {
                        if((buf & 3) != 1)
                        {
                            // g is {0, 0.5, 1}
                            F(snp_idx) += BED2GENO[buf & 3];
                            c++;
                        }
                        buf >>= 2; // shift packed data and throw away genotype just processed.
                    }
                }
            }
            // calculate F and centered_geno_lookup
            if(c == 0)
            {
                llog << colwarn << "the allele frequency should not be 0. snp index:" << snp_idx << colend;
                F(snp_idx) = 0;
            }
            else
            {
                F(snp_idx) /= c;
            }
            // do centering and initialing
            centered_geno_lookup(1, snp_idx) = 0.0; // missing
            centered_geno_lookup(0, snp_idx) = BED2GENO[0] - F(snp_idx); // minor hom
            centered_geno_lookup(2, snp_idx) = BED2GENO[2] - F(snp_idx); // het
            centered_geno_lookup(3, snp_idx) = BED2GENO[3] - F(snp_idx); // major hom
            // get centered and standardized G
            for(b = 0, j = 0; b < bed_bytes_per_snp; ++b)
            {
                buf = inbed[i * bed_bytes_per_snp + b];
                for(k = 0; k < 4; ++k, ++j)
                {
                    if(j < nsamples)
                    {
                        G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
                        if(standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
                            G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
                        buf >>= 2; // shift packed data and throw away genotype just processed.
                    }
                }
            }
        }
    }

    if(stop_idx + 1 == nsnps) frequency_was_estimated = true;
}

void FileBed::read_block_update(uint64 start_idx,
                                uint64 stop_idx,
                                const MyMatrix & U,
                                const MyVector & svals,
                                const MyMatrix & VT,
                                bool standardize)
{
    uint actual_block_size = stop_idx - start_idx + 1;
    if(G.cols() < params.blocksize || (actual_block_size < params.blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
        inbed.reserve(bed_bytes_per_snp * params.blocksize);
    }
    // check where we are
    if(params.verbose)
    {
        long long offset = 3 + start_idx * bed_bytes_per_snp;
        if(bed_ifstream.tellg() != offset)
        {
            throw std::runtime_error("Error: something wrong with read_snp_block!\n");
        }
    }
    uint64 b, i, j, snp_idx;
    uint ks = svals.rows();
    uint ki, k;
    uchar buf;
    bed_ifstream.read(reinterpret_cast<char *>(&inbed[0]), bed_bytes_per_snp * actual_block_size);
#pragma omp parallel for private(i, j, b, ki, k, snp_idx, buf)
    for(i = 0; i < actual_block_size; ++i)
    {
        snp_idx = start_idx + i;
        for(b = 0, j = 0; b < bed_bytes_per_snp; ++b)
        {
            buf = inbed[i * bed_bytes_per_snp + b];
            for(ki = 0; ki < 4; ++ki, ++j)
            {
                if(j < nsamples)
                {
                    G(j, i) = centered_geno_lookup(buf & 3, snp_idx);
                    if(BED2GENO[buf & 3] == BED_MISSING_VALUE)
                    {
                        G(j, i) = 0.0;
                        for(k = 0; k < ks; ++k)
                        {
                            G(j, i) += U(j, k) * svals(k) * VT(k, snp_idx);
                        }
                        // map to domain(0,1)
                        G(j, i) = fmin(fmax(G(j, i), -F(snp_idx)), 1 - F(snp_idx));
                    }
                    if(standardize && sqrt(F(snp_idx) * (1 - F(snp_idx))) > VAR_TOL)
                        G(j, i) /= sqrt(F(snp_idx) * (1 - F(snp_idx)));
                    buf >>= 2; // shift packed data and throw away genotype just processed.
                }
            }
        }
    }
}

// structured permutation with cached buffer
void permute_plink(std::string & fin, const std::string & fout, uint gb, uint nbands)
{
    uint nsnps = count_lines(fin + ".bim");
    uint nsamples = count_lines(fin + ".fam");
    uint bed_bytes_per_snp = (nsamples + 3) >> 2;
    cout << timestamp() << "permute plink files. nsnps:" << nsnps << ", nsamples:" << nsamples << endl;

    // calculate the readin number of snps of certain big buffer like 2GB.
    // must be a multiple of nbands.
    uint twoGB_snps = (uint)floor((double)1073741824 * gb / bed_bytes_per_snp);
    if(twoGB_snps > nsnps) twoGB_snps = nsnps;
    uint bufsize = (uint)floor((double)twoGB_snps / nbands);
    twoGB_snps = bufsize * nbands; // initially twoGB_snps is a multiple of nbands
    assert(nsnps >= twoGB_snps);
    uint nblocks = (nsnps + twoGB_snps - 1) / twoGB_snps;
    uint modr2 = nsnps % twoGB_snps;
    uint64 bed_bytes_per_block = bed_bytes_per_snp * twoGB_snps;
    vector<uchar> inbed; // keep the input buffer
    inbed.resize(bed_bytes_per_block);
    vector<uchar> outbed; // keep the output buffer
    uint64 out_bytes_per_block = bed_bytes_per_snp * bufsize;
    outbed.resize(out_bytes_per_block);

    // get index of first snp of each band
    vector<uint64> bandidx;
    bandidx.resize(nbands);
    uint modr = nsnps % nbands;
    uint bandsize = (nsnps + nbands - 1) / nbands;
    if(modr == 0)
    {
        for(uint i = 0; i < nbands; ++i)
        {
            bandidx[i] = i * bandsize;
        }
    }
    else
    {
        for(uint i = 0; i < nbands; ++i)
        {
            if(i < modr)
            {
                bandidx[i] = i * bandsize;
            }
            else
            {
                bandidx[i] = modr * bandsize + (bandsize - 1) * (i - modr);
            }
        }
    }

    ios_base::sync_with_stdio(false);
    std::ifstream in(fin + ".bed", std::ios::binary);
    std::ofstream out(fout + ".perm.bed", std::ios::binary);
    if(!in.is_open())
    {
        throw std::invalid_argument(colerror + "Cannot open bed file.\n");
    }
    uchar header[3];
    in.read(reinterpret_cast<char *>(&header[0]), 3);
    if((header[0] != 0x6c) || (header[1] != 0x1b) || (header[2] != 0x01))
    {
        throw std::invalid_argument(colerror + "Incorrect magic number in bed file.\n");
    }
    out.write(reinterpret_cast<char *>(&header[0]), 3);
    std::ifstream in_bim(fin + ".bim", std::ios::in);
    std::ofstream out_bim(fout + ".perm.bim", std::ios::out);
    vector<std::string> bims(std::istream_iterator<Line>{in_bim}, std::istream_iterator<Line>{});
    vector<std::string> bims2;
    bims2.resize(nsnps);
    uint64 b, i, j, twoGB_snps2, idx, bufidx = bufsize;
    for(i = 0; i < nblocks; i++)
    {
        if(i == nblocks - 1 && modr2 != 0)
        {
            twoGB_snps2 = nsnps - (nblocks - 1) * twoGB_snps;
            bed_bytes_per_block = bed_bytes_per_snp * twoGB_snps2;
            inbed.resize(bed_bytes_per_block);
            // in last block, twoGB_snps is not neccessary a multiple of nbands and smaller than the previous
            bufsize = (uint64)(twoGB_snps2 + nbands - 1) / nbands;
            modr2 = twoGB_snps2 % nbands;
            out_bytes_per_block = bed_bytes_per_snp * bufsize;
            outbed.resize(out_bytes_per_block);
        }
        in.read(reinterpret_cast<char *>(&inbed[0]), bed_bytes_per_block);
        for(b = 0; b < nbands; b++)
        {
            idx = 3 + (i * bufidx + bandidx[b]) * bed_bytes_per_snp;
            for(j = 0; j < bufsize - 1; j++)
            {
                std::copy(inbed.begin() + (j * nbands + b) * bed_bytes_per_snp,
                          inbed.begin() + (j * nbands + b + 1) * bed_bytes_per_snp,
                          outbed.begin() + j * bed_bytes_per_snp);
                // cout << i * twoGB_snps + j * nbands + b << endl;
                bims2[i * bufidx + bandidx[b] + j] = bims[i * twoGB_snps + j * nbands + b];
            }
            if(i != nblocks - 1 || (i == nblocks - 1 && b < modr2) || modr2 == 0)
            {
                std::copy(inbed.begin() + (j * nbands + b) * bed_bytes_per_snp,
                          inbed.begin() + (j * nbands + b + 1) * bed_bytes_per_snp,
                          outbed.begin() + j * bed_bytes_per_snp);
                bims2[i * bufidx + bandidx[b] + j] = bims[i * twoGB_snps + j * nbands + b];
            }
            else
            {
                out_bytes_per_block = bed_bytes_per_snp * (bufsize - 1);
            }
            out.seekp(idx, std::ios_base::beg);
            out.write(reinterpret_cast<char *>(&outbed[0]), out_bytes_per_block);
        }
    }
    for(auto b : bims2)
    {
        out_bim << b << "\n";
    }
    in.close();
    out.close();
    out_bim.close();

    std::ifstream in_fam(fin + ".fam");
    std::ofstream out_fam(fout + ".perm.fam");
    out_fam << in_fam.rdbuf();
    fin = fout + ".perm";
}
