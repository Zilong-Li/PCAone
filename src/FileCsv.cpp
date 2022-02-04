#include "FileCsv.hpp"

using namespace std;

void FileCsv::read_all_and_centering()
{
    G = MyMatrix(nsamples, nsnps);

    auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));
    size_t read, i, j, e, lastSNP = 0;
    fin = fopenOrDie(params.csvfile.c_str(), "rb");
    buffCur = "";

    while ((read = freadOrDie(buffIn, buffInSize, fin))) {
        ZSTD_inBuffer input = {buffIn, read, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
            lastRet = ZSTD_decompressStream(dctx, &output, &input);
            buffCur += std::string((char *)buffOut, output.pos);
            while ( (e = buffCur.find("\n")) != std::string::npos ) {
                buffLine = buffCur.substr(0, e);
                buffCur.erase(0, e + 1);
                for (i = 0, j = 1; i < buffLine.size(); i ++) {
                    if (buffLine[i] == ',') {
                        tidx[j] = i + 1;
                        j++;
                    }
                }
                tidx[nsamples] = buffLine.size() + 1;

                #pragma omp parallel for
                for (size_t i = 0; i < nsamples; i++) {
                    G(i, lastSNP) = std::stod(buffLine.substr(tidx[i], tidx[i+1] - tidx[i] - 1));
                }

                G.col(lastSNP).array() -= G.col(lastSNP).mean();  // only do centering
                lastSNP++;
            }
        }
    }

    if (lastRet != 0) {
        throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");
    }

    // deal with the case there is no "\n" for the last line of file

    if (lastSNP != nsnps) {
        throw std::runtime_error("error when parsing csv file\n");
    }

}

void FileCsv::check_file_offset_first_var()
{
    if ( fin == nullptr) {
        fin = fopenOrDie(params.csvfile.c_str(), "rb");
    } else if (feof(fin) || lastRet == 0) {
        rewind(fin);
    } else {
        throw std::runtime_error("no eof detected. something wrong.\n");
    }
    lastRet = 1;
    buffCur = "";

}

void FileCsv::read_snp_block_initial(uint64 start_tidx, uint64 stop_tidx, bool standardize)
{

    const uint actual_block_size = stop_tidx - start_tidx + 1;

    if (G.cols() < params.blocksize || (actual_block_size < params.blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
    }
    auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));
    size_t read, i, j, e, lastSNP = 0;

    if (buffCur != "") {
        while (lastSNP < actual_block_size && ( (e = buffCur.find("\n")) != std::string::npos )) {
            buffLine = buffCur.substr(0, e);
            buffCur.erase(0, e + 1);
            for (i = 0, j = 1; i < buffLine.size(); i ++) {
                if (buffLine[i] == ',') {
                    tidx[j++] = i + 1;
                }
            }
            tidx[nsamples] = buffLine.size() + 1;

            #pragma omp parallel for
            for (size_t i = 0; i < nsamples; i++) {
                G(i, lastSNP) = std::stod(buffLine.substr(tidx[i], tidx[i+1] - tidx[i] - 1));
            }

            G.col(lastSNP).array() -= G.col(lastSNP).mean();  // only do centering
            lastSNP++;
        }
    }


    if (lastRet != 0 && lastSNP < actual_block_size) {
        while ((read = freadOrDie(buffIn, buffInSize, fin)) ) {
            ZSTD_inBuffer input = {buffIn, read, 0};
            while (input.pos < input.size) {
                ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
                lastRet = ZSTD_decompressStream(dctx, &output, &input);
                buffCur += std::string((char *)buffOut, output.pos);
                while (lastSNP < actual_block_size && ( (e = buffCur.find("\n")) != std::string::npos )) {
                    buffLine = buffCur.substr(0, e);
                    buffCur.erase(0, e + 1);
                    for (i = 0, j = 1; i < buffLine.size(); i ++) {
                        if (buffLine[i] == ',') {
                            tidx[j] = i + 1;
                            j++;
                        }
                    }
                    tidx[nsamples] = buffLine.size() + 1;

                    #pragma omp parallel for
                    for (size_t i = 0; i < nsamples; i++) {
                        G(i, lastSNP) = std::stod(buffLine.substr(tidx[i], tidx[i+1] - tidx[i] - 1));
                    }

                    G.col(lastSNP).array() -= G.col(lastSNP).mean();  // only do centering
                    lastSNP++;
                }
            }
            if (lastSNP >= actual_block_size)
                break;
        }
    }

    if (lastSNP != actual_block_size) {
        throw std::runtime_error("something wrong when read_snp_block_initial");
    }
}
