#include "FileCsv.hpp"

using namespace std;

void FileCsv::read_all_and_centering()
{
    cout << timestamp() << "begin to read whole data" << endl;
    std::vector<std::string> buffVecs;
    size_t const buffInSize = ZSTD_DStreamInSize();
    buffInTmp.reserve(buffInSize);
    auto buffIn = const_cast<void *>(static_cast<const void *>(buffInTmp.c_str()));
    auto buffOutSize = ZSTD_DStreamOutSize();
    buffOutTmp.reserve(buffOutSize);
    auto buffOut = const_cast<void *>(static_cast<const void *>(buffOutTmp.c_str()));
    ZSTD_DCtx *const dctx = ZSTD_createDCtx();
    size_t const toRead = buffInSize;
    size_t lastRet = 0;
    size_t read, p;
    fin = fopen_orDie(params.csvfile.c_str(), "rb");
    buffCur = "";
    while ((read = fread_orDie(buffIn, toRead, fin))) {
        ZSTD_inBuffer input = {buffIn, read, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
            lastRet = ZSTD_decompressStream(dctx, &output, &input);
            buffCur += std::string((char *)buffOut, output.pos);
            while (( p = buffCur.find("\n") ) != std::string::npos) {
                nsamples++;
                buffLine = buffCur.substr(0, p);
                buffCur.erase(0, p + 1);
                buffVecs.push_back(buffLine);
            }
        }
    }
    buffVecs.push_back(buffCur);
    ZSTD_freeDCtx(dctx);
    fclose_orDie(fin);

    assert(buffVecs.size() == nsamples);

    G = MatrixXd(nsamples, nsnps);
    // start parsing all inputs in parallel
    #pragma omp parallel for
    for(size_t i = 0; i < nsamples; i++) {
        std::string line = buffVecs[i];
        size_t j = 0;
        while (( p = line.find(",") ) != std::string::npos) {
            G(i, j) = std::stod(line.substr(0, p));
            j++;
            line.erase(0, p + 1);
        }
        G(i, j) = std::stod(line);
    }
}
