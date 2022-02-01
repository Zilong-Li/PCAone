#include "FileCsv.hpp"

using namespace std;

void FileCsv::read_all_and_centering()
{
    auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));
    size_t read, p;
    fin = fopenOrDie(params.csvfile.c_str(), "rb");
    buffCur = "";
    while ((read = freadOrDie(buffIn, buffInSize, fin))) {
        ZSTD_inBuffer input = {buffIn, read, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
            lastRet = ZSTD_decompressStream(dctx, &output, &input);
            buffCur += std::string((char *)buffOut, output.pos);
            while (( p = buffCur.find("\n") ) != std::string::npos) {
                buffLine = buffCur.substr(0, p);
                buffVecs.push_back(buffLine);
                buffCur.erase(0, p + 1);
            }
        }
    }

    if (lastRet != 0) {
        throw std::runtime_error("EOF before end of ZSTD_decompressStream.\n");
    }

    // in case there is no "\n" for the last line of file
    if (buffCur != "") {
        buffVecs.push_back(buffCur);
    }

    if (buffVecs.size() != nsnps) {
        throw std::runtime_error("error when parsing csv file\n");
    }

    G = MyMatrix(nsamples, nsnps);
    // start parsing all inputs in parallel
    #pragma omp parallel for
    for(size_t i = 0; i < nsnps; i++) {
        size_t j = 0, p;
        std::string line = buffVecs[i];
        while (( p = line.find(",") ) != std::string::npos) {
            G(j, i) = std::stod(line.substr(0, p));
            j++;
            line.erase(0, p + 1);
        }
        G(j, i) = std::stod(line);
        G.col(i).array() -= G.col(i).mean();  // centeringg
    }
    buffVecs.clear();
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

}

void FileCsv::read_snp_block_initial(uint64 start_idx, uint64 stop_idx, bool standardize)
{

    const uint actual_block_size = stop_idx - start_idx + 1;

    if (G.cols() < params.blocksize || (actual_block_size < params.blocksize))
    {
        G = MyMatrix::Zero(nsamples, actual_block_size);
    }
    auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));
    auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));
    size_t read, p, lastRead = 0;
    if (buffCur != "") {
        while (( p = buffCur.find("\n") ) != std::string::npos && lastRead < actual_block_size) {
            lastRead++;
            buffLine = buffCur.substr(0, p);
            buffVecs.push_back(buffLine);
            buffCur.erase(0, p + 1);
        }
    }

    if (lastRet != 0 && lastRead < actual_block_size) {
        while ((read = freadOrDie(buffIn, buffInSize, fin)) ) {
            ZSTD_inBuffer input = {buffIn, read, 0};
            while (input.pos < input.size) {
                ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
                lastRet = ZSTD_decompressStream(dctx, &output, &input);
                buffCur += std::string((char *)buffOut, output.pos);
                while (( p = buffCur.find("\n") ) != std::string::npos && lastRead < actual_block_size) {
                    lastRead++;
                    buffLine = buffCur.substr(0, p);
                    buffVecs.push_back(buffLine);
                    buffCur.erase(0, p + 1);
                }
            }
            if (buffVecs.size() >= actual_block_size)
                break;
        }
    }

    if (buffVecs.size() != actual_block_size) {
        cout << buffVecs.size() << "\t" << actual_block_size << "\t" << lastRead << endl;
        throw std::runtime_error("something wrong when read_snp_block_initial");
    }

    #pragma omp parallel for
    for (size_t i = 0; i < actual_block_size; ++i) {
        size_t j = 0, p;
        std::string line = buffVecs[i];
        while (( p = line.find(",") ) != std::string::npos) {
            G(j, i) = std::stod(line.substr(0, p));
            j++;
            line.erase(0, p + 1);
        }
        G(j, i) = std::stod(line);
        G.col(i).array() -= G.col(i).mean();  // centeringg
    }

    buffVecs.clear();
}
