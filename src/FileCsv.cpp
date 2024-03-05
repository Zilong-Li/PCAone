#include "FileCsv.hpp"

using namespace std;

void FileCsv::read_all()
{
    auto buffIn = const_cast<void *>(static_cast<const void *>(zbuf.buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(zbuf.buffOutTmp.c_str()));
    size_t read, i, j, e, lastSNP = 0;
    zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
    zbuf.buffCur = "";
    G = MyMatrix::Zero(nsamples, nsnps);
    while((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin)))
    {
        ZSTD_inBuffer input = {buffIn, read, 0};
        while(input.pos < input.size)
        {
            ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
            zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
            zbuf.buffCur += std::string((char *)buffOut, output.pos);
            while((e = zbuf.buffCur.find("\n")) != std::string::npos)
            {
                zbuf.buffLine = zbuf.buffCur.substr(0, e);
                zbuf.buffCur.erase(0, e + 1);
                for(i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                {
                    if(zbuf.buffLine[i] == ',')
                    {
                        tidx[j] = i + 1;
                        j++;
                    }
                }
                tidx[nsamples] = zbuf.buffLine.size() + 1;

#pragma omp parallel for
                for(size_t i = 0; i < nsamples; i++)
                {
                    G(i, lastSNP) = std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                    if(params.scale == 1)
                        G(i, lastSNP) = log10(G(i, lastSNP) + 0.01);
                    else if(params.scale == 2)
                        G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
                }

                // G.col(lastSNP).array() -= G.col(lastSNP).mean(); // only do centering
                lastSNP++;
            }
        }
    }

    G.rowwise() -= G.colwise().mean(); // do centering

    if(zbuf.lastRet != 0) throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");

    // deal with the case there is no "\n" for the last line of file
    if(lastSNP != nsnps) throw std::runtime_error("error when parsing csv file\n");
}

void FileCsv::check_file_offset_first_var()
{
    if(zbuf.fin == nullptr)
    {
        zbuf.fin = fopenOrDie(params.filein.c_str(), "rb");
    }
    else if(feof(zbuf.fin) || zbuf.lastRet == 0)
    {
        rewind(zbuf.fin);
    }
    else
    {
        rewind(zbuf.fin);
        if(params.verbose) cao.warning("make sure you are runing PCAone algorithm2");
    }
    zbuf.lastRet = 1;
    zbuf.buffCur = "";
}

void FileCsv::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    read_csvzstd_block(zbuf, params.blocksize, start_idx, stop_idx, G, nsamples, libsize, tidx,
                       median_libsize, params.scale);
}

void parse_csvzstd(ZstdBuffer & zbuf,
                   uint & nsamples,
                   uint & nsnps,
                   uint scale,
                   std::vector<double> & libsize,
                   std::vector<size_t> & tidx,
                   double & median_libsize)
{
    auto buffIn = const_cast<void *>(static_cast<const void *>(zbuf.buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(zbuf.buffOutTmp.c_str()));
    size_t read, i, j, p, ncol = 0, lastCol = 0;
    int isEmpty = 1;
    nsnps = 0;
    while((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin)))
    {
        isEmpty = 0;
        ZSTD_inBuffer input = {buffIn, read, 0};
        while(input.pos < input.size)
        {
            ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
            zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
            zbuf.buffCur += std::string((char *)buffOut, output.pos);
            while((p = zbuf.buffCur.find("\n")) != std::string::npos)
            {
                nsnps++;
                zbuf.buffLine = zbuf.buffCur.substr(0, p);
                zbuf.buffCur.erase(0, p + 1);
                lastCol = ncol;
                ncol = 1;
                for(i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                {
                    if(zbuf.buffLine[i] == ',')
                    {
                        ncol++;
                        if(nsnps > 1)
                        {
                            tidx[j++] = i + 1;
                        }
                    }
                }
                // get ncol from the first line or header
                if(nsnps == 1)
                {
                    libsize.resize(ncol);
                    tidx.resize(ncol + 1);
                    for(i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                    {
                        if(zbuf.buffLine[i] == ',')
                        {
                            tidx[j++] = i + 1;
                        }
                    }
                }

                tidx[ncol] = zbuf.buffLine.size() + 1;
                if(scale == 2) // cpmed
                {
#pragma omp parallel for
                    for(size_t i = 0; i < ncol; i++)
                    {
                        libsize[i] += std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                    }
                }
                if(nsnps > 2 && (lastCol != ncol)) cao.error("the csv file has unaligned columns");
            }
        }
    }

    if(isEmpty) cao.error("input file is empty.");
    if(zbuf.lastRet != 0) throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");

    nsamples = ncol;
    zbuf.lastRet = 1;
    if(scale == 2) median_libsize = get_median(libsize);
}

void read_csvzstd_block(ZstdBuffer & zbuf,
                        int blocksize,
                        uint64 start_idx,
                        uint64 stop_idx,
                        MyMatrix & G,
                        uint nsamples,
                        std::vector<double> & libsize,
                        std::vector<size_t> & tidx,
                        double median_libsize,
                        uint scale)
{
    const uint actual_block_size = stop_idx - start_idx + 1;

    if(G.cols() < blocksize || (actual_block_size < blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
    }
    auto buffIn = const_cast<void *>(static_cast<const void *>(zbuf.buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(zbuf.buffOutTmp.c_str()));
    size_t read, i, j, e, lastSNP = 0;
    if(zbuf.buffCur != "")
    {
        while(lastSNP < actual_block_size && ((e = zbuf.buffCur.find("\n")) != std::string::npos))
        {
            zbuf.buffLine = zbuf.buffCur.substr(0, e);
            zbuf.buffCur.erase(0, e + 1);
            for(i = 0, j = 1; i < zbuf.buffLine.size(); i++)
            {
                if(zbuf.buffLine[i] == ',')
                {
                    tidx[j++] = i + 1;
                }
            }
            tidx[nsamples] = zbuf.buffLine.size() + 1;

#pragma omp parallel for
            for(size_t i = 0; i < nsamples; i++)
            {
                G(i, lastSNP) = std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                if(scale == 1)
                    G(i, lastSNP) = log10(G(i, lastSNP) + 0.01);
                else if(scale == 2)
                    G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
            }

            lastSNP++;
            G.col(lastSNP).array() -= G.col(lastSNP).mean(); // only do centering
        }
    }

    if(zbuf.lastRet != 0 && lastSNP < actual_block_size)
    {
        while((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin)))
        {
            ZSTD_inBuffer input = {buffIn, read, 0};
            while(input.pos < input.size)
            {
                ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
                zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
                zbuf.buffCur += std::string((char *)buffOut, output.pos);
                while(lastSNP < actual_block_size && ((e = zbuf.buffCur.find("\n")) != std::string::npos))
                {
                    zbuf.buffLine = zbuf.buffCur.substr(0, e);
                    zbuf.buffCur.erase(0, e + 1);
                    for(i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                    {
                        if(zbuf.buffLine[i] == ',')
                        {
                            tidx[j] = i + 1;
                            j++;
                        }
                    }
                    tidx[nsamples] = zbuf.buffLine.size() + 1;

#pragma omp parallel for
                    for(size_t i = 0; i < nsamples; i++)
                    {
                        G(i, lastSNP) = std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                        if(scale == 1)
                            G(i, lastSNP) = log10(G(i, lastSNP) + 0.01);
                        else if(scale == 2)
                            G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
                    }

                    G.col(lastSNP).array() -= G.col(lastSNP).mean(); // only do centering
                    lastSNP++;
                }
            }
            if(lastSNP >= actual_block_size) break;
        }
    }

    if(lastSNP != actual_block_size) cao.error("something wrong when read_block_initial");
}

PermMat shuffle_csvzstd_to_bin(std::string & fin, std::string fout, uint gb, uint scale)
{
    std::vector<size_t> tidx;
    std::vector<double> libsize;
    double median_libsize;
    uint nsnps, nsamples;
    const uint ibyte = 4;
    ZstdBuffer zbuf;
    {
        zbuf.fin = fopenOrDie(fin.c_str(), "rb");
        parse_csvzstd(zbuf, nsamples, nsnps, scale, libsize, tidx, median_libsize);
        fcloseOrDie(zbuf.fin);
    }
    uint64 bytes_per_snp = nsamples * ibyte;
    uint blocksize = 1073741824 * gb / bytes_per_snp;
    uint nblocks = (nsnps + blocksize - 1) / blocksize;
    std::ofstream ofs(fout + ".perm.bin", std::ios::binary);
    std::ofstream ofs2(fout + ".perm.txt");
    ofs.write((char *)&nsnps, ibyte);
    ofs.write((char *)&nsamples, ibyte);
    uint magic = ibyte * 2;
    zbuf.fin = fopenOrDie(fin.c_str(), "rb");
    zbuf.lastRet = 1;
    zbuf.buffCur = "";
    MyMatrix G;
    std::vector<int> perm(nsnps);
    std::iota(perm.begin(), perm.end(), 0);
    auto rng = std::default_random_engine{};
    std::shuffle(perm.begin(), perm.end(), rng);
    Eigen::VectorXf fg;
    uint64 start_idx, stop_idx, idx;
    int ia{0}, ib{0};
    Eigen::VectorXi indices(nsnps);
    for(uint i = 0; i < nblocks; i++)
    {
        start_idx = i * blocksize;
        stop_idx = start_idx + blocksize - 1;
        stop_idx = stop_idx >= nsnps ? nsnps - 1 : stop_idx;
        read_csvzstd_block(zbuf, blocksize, start_idx, stop_idx, G, nsamples, libsize, tidx, median_libsize,
                           scale);
        for(Eigen::Index p = 0; p < G.cols(); p++, ia++)
        {
            ib = perm[ia];
            indices(ib) = ia;
            idx = magic + ib * bytes_per_snp;
            ofs.seekp(idx, std::ios_base::beg);
            fg = G.col(p).cast<float>();
            ofs.write((char *)fg.data(), bytes_per_snp);
        }
    }
    fin = fout + ".perm.bin";
    ofs2 << indices << "\n";
    return PermMat(indices);
}
