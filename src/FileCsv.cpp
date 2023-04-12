#include "FileCsv.hpp"

using namespace std;

void FileCsv::read_all()
{
    auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
    auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
    size_t read, i, j, e, lastSNP = 0;
    zbuf.fin = fopenOrDie(params.csvfile.c_str(), "rb");
    zbuf.buffCur = "";
    G = MyMatrix::Zero(nsamples, nsnps);
    while ((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin)))
    {
        ZSTD_inBuffer input = {buffIn, read, 0};
        while (input.pos < input.size)
        {
            ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
            zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
            zbuf.buffCur += std::string((char*)buffOut, output.pos);
            while ((e = zbuf.buffCur.find("\n")) != std::string::npos)
            {
                zbuf.buffLine = zbuf.buffCur.substr(0, e);
                zbuf.buffCur.erase(0, e + 1);
                for (i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                {
                    if (zbuf.buffLine[i] == ',')
                    {
                        tidx[j] = i + 1;
                        j++;
                    }
                }
                tidx[nsamples] = zbuf.buffLine.size() + 1;

#pragma omp parallel for
                for (size_t i = 0; i < nsamples; i++)
                {
                    G(i, lastSNP) = std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                    if (params.cpmed)
                        G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
                }

                // G.col(lastSNP).array() -= G.col(lastSNP).mean(); // only do centering
                lastSNP++;
            }
        }
    }

    if (params.center)
        G.rowwise() -= G.colwise().mean();

    if (zbuf.lastRet != 0)
        throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");

    // deal with the case there is no "\n" for the last line of file
    if (lastSNP != nsnps)
        throw std::runtime_error("error when parsing csv file\n");
}

void FileCsv::check_file_offset_first_var()
{
    if (zbuf.fin == nullptr)
    {
        zbuf.fin = fopenOrDie(params.csvfile.c_str(), "rb");
    }
    else if (feof(zbuf.fin) || zbuf.lastRet == 0)
    {
        rewind(zbuf.fin);
    }
    else
    {
        rewind(zbuf.fin);
        if (params.verbose)
            std::cout << colwarn + "make sure you are runing PCAone algorithm2" + colend << std::endl;
    }
    zbuf.lastRet = 1;
    zbuf.buffCur = "";
}

void FileCsv::read_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{
    read_csvzstd_block(zbuf, start_idx, stop_idx, standardize, G, params.blocksize, nsamples, params.cpmed, libsize, tidx, median_libsize);
}

void parse_csvzstd(ZstdBuffer& zbuf, uint64& nsamples, uint64& nsnps, bool cpmed, std::vector<double>& libsize, std::vector<size_t>& tidx,
                   double& median_libsize)
{
    auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
    auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
    size_t read, i, j, p, ncol = 0, lastCol = 0;
    int isEmpty = 1;
    nsnps = 0;
    while (zbuf.lastRet != 0 && (read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin)))
    {
        isEmpty = 0;
        ZSTD_inBuffer input = {buffIn, read, 0};
        while (input.pos < input.size)
        {
            ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
            zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
            zbuf.buffCur += std::string((char*)buffOut, output.pos);
            while ((p = zbuf.buffCur.find("\n")) != std::string::npos)
            {
                nsnps++;
                zbuf.buffLine = zbuf.buffCur.substr(0, p);
                zbuf.buffCur.erase(0, p + 1);
                lastCol = ncol;
                ncol = 1;
                for (i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                {
                    if (zbuf.buffLine[i] == ',')
                    {
                        ncol++;
                        if (nsnps > 1 && cpmed)
                        {
                            tidx[j++] = i + 1;
                        }
                    }
                }
                // get ncol from the first line or header
                if (nsnps == 1 && cpmed)
                {
                    libsize.resize(ncol);
                    tidx.resize(ncol + 1);
                    for (i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                    {
                        if (zbuf.buffLine[i] == ',')
                        {
                            tidx[j++] = i + 1;
                        }
                    }
                }

                if (cpmed)
                {
                    tidx[ncol] = zbuf.buffLine.size() + 1;
#pragma omp parallel for
                    for (size_t i = 0; i < ncol; i++)
                    {
                        libsize[i] += std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                    }
                }

                if (nsnps > 2 && (lastCol != ncol))
                {
                    throw std::invalid_argument(colerror + "the csv file has unaligned columns\n");
                }
            }
        }
    }

    if (isEmpty)
        throw std::invalid_argument(colerror + "input file is empty.\n");
    if (zbuf.lastRet != 0)
        throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");

    nsamples = ncol;
    zbuf.lastRet = 1;
    if (cpmed)
        median_libsize = get_median(libsize);
}

void read_csvzstd_block(ZstdBuffer& zbuf, uint64 start_idx, uint64 stop_idx, bool standardize, MyMatrix& G, int blocksize, uint64 nsamples, bool cpmed,
                        std::vector<double>& libsize, std::vector<size_t>& tidx, double median_libsize)
{

    const uint actual_block_size = stop_idx - start_idx + 1;

    if (G.cols() < blocksize || (actual_block_size < blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
    }
    auto buffIn = const_cast<void*>(static_cast<const void*>(zbuf.buffInTmp.c_str()));
    auto buffOut = const_cast<void*>(static_cast<const void*>(zbuf.buffOutTmp.c_str()));
    size_t read, i, j, e, lastSNP = 0;
    if (zbuf.buffCur != "")
    {
        while (lastSNP < actual_block_size && ((e = zbuf.buffCur.find("\n")) != std::string::npos))
        {
            zbuf.buffLine = zbuf.buffCur.substr(0, e);
            zbuf.buffCur.erase(0, e + 1);
            for (i = 0, j = 1; i < zbuf.buffLine.size(); i++)
            {
                if (zbuf.buffLine[i] == ',')
                {
                    tidx[j++] = i + 1;
                }
            }
            tidx[nsamples] = zbuf.buffLine.size() + 1;

#pragma omp parallel for
            for (size_t i = 0; i < nsamples; i++)
            {
                G(i, lastSNP) = std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                if (cpmed)
                    G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
            }

            G.col(lastSNP).array() -= G.col(lastSNP).mean(); // only do centering
            lastSNP++;
        }
    }

    if (zbuf.lastRet != 0 && lastSNP < actual_block_size)
    {
        while ((read = freadOrDie(buffIn, zbuf.buffInSize, zbuf.fin)))
        {
            ZSTD_inBuffer input = {buffIn, read, 0};
            while (input.pos < input.size)
            {
                ZSTD_outBuffer output = {buffOut, zbuf.buffOutSize, 0};
                zbuf.lastRet = ZSTD_decompressStream(zbuf.dctx, &output, &input);
                zbuf.buffCur += std::string((char*)buffOut, output.pos);
                while (lastSNP < actual_block_size && ((e = zbuf.buffCur.find("\n")) != std::string::npos))
                {
                    zbuf.buffLine = zbuf.buffCur.substr(0, e);
                    zbuf.buffCur.erase(0, e + 1);
                    for (i = 0, j = 1; i < zbuf.buffLine.size(); i++)
                    {
                        if (zbuf.buffLine[i] == ',')
                        {
                            tidx[j] = i + 1;
                            j++;
                        }
                    }
                    tidx[nsamples] = zbuf.buffLine.size() + 1;

#pragma omp parallel for
                    for (size_t i = 0; i < nsamples; i++)
                    {
                        G(i, lastSNP) = std::stod(zbuf.buffLine.substr(tidx[i], tidx[i + 1] - tidx[i] - 1));
                        if (cpmed)
                            G(i, lastSNP) = log10(G(i, lastSNP) * median_libsize / libsize[i] + 1);
                    }

                    G.col(lastSNP).array() -= G.col(lastSNP).mean(); // only do centering
                    lastSNP++;
                }
            }
            if (lastSNP >= actual_block_size)
                break;
        }
    }

    if (lastSNP != actual_block_size)
    {
        throw std::runtime_error("something wrong when read_block_initial");
    }
}

void shuffle_csvzstd_to_bin(std::string csvfile, std::string binfile, int blocksize, bool standardize, bool cpmed)
{
    std::vector<size_t> tidx;
    std::vector<double> libsize;
    double median_libsize;
    uint64 nsnps, nsamples, idx, cur = 0;
    ZstdBuffer zbuf;
    {
        zbuf.fin = fopenOrDie(csvfile.c_str(), "rb");
        parse_csvzstd(zbuf, nsamples, nsnps, cpmed, libsize, tidx, median_libsize);
        fcloseOrDie(zbuf.fin);
    }
    uint nblocks = (nsnps + blocksize - 1) / blocksize;
    std::ofstream ofs(binfile, std::ios::binary);
    ofs.write((char*)&nsamples, sizeof(nsamples));
    ofs.write((char*)&nsnps, sizeof(nsnps));
    zbuf.fin = fopenOrDie(csvfile.c_str(), "rb");
    zbuf.lastRet = 1;
    zbuf.buffCur = "";
    MyMatrix G;
    std::vector<uint64> perm(nsnps);
    auto rng = std::default_random_engine{};
    std::uniform_int_distribution<> dis(0, nsnps - 1);
    for (auto& p : perm)
        p = dis(rng);
    for (uint i = 0; i < nblocks; i++)
    {
        auto start_idx = i * blocksize;
        auto stop_idx = start_idx + blocksize - 1;
        stop_idx = stop_idx >= nsnps ? nsnps - 1 : stop_idx;
        read_csvzstd_block(zbuf, start_idx, stop_idx, standardize, G, blocksize, nsamples, cpmed, libsize, tidx, median_libsize);
        for (size_t p = 0; p < G.cols(); p++, cur++)
        {
            idx = 2 * sizeof(uint64) + perm[cur] * nsamples * sizeof(double);
            ofs.seekp(idx, std::ios_base::beg);
            ofs.write((char*)G.col(p).data(), nsamples * sizeof(double));
        }
    }
    // std::cout << nsnps << "," << cur << "\n";
}
