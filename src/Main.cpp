#include "Arnoldi.hpp"
#include "Cmd.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileBinary.hpp"
#include "FileCsv.hpp"
#include "FilePlink.hpp"
#include "Halko.hpp"
#include <omp.h>

#ifdef WITH_OPENBLAS
#    include "lapacke.h"
#elif defined WITH_MKL
#    include "mkl_lapacke.h"
#endif

using namespace std;

int main(int argc, char * argv[])
{
    auto t1 = std::chrono::steady_clock::now();
    Param params(argc, argv);
    string commandargs = params.ss.str();
    // set number of threads
    // openblas_set_num_threads(params.threads);
    omp_set_num_threads(params.threads);
    Data * data;
    if(params.svd_t == SvdType::PCAoneAlg2 && !params.noshuffle && params.out_of_core)
    {
        auto ts = std::chrono::steady_clock::now();
        if(params.file_t == FileType::PLINK)
        {
            permute_plink(params.filein, params.fileout, params.buffer, params.bands);
            data = new FileBed(params);
        }
        else if(params.file_t == FileType::BGEN)
        {
            permute_bgen(params.filein, params.fileout);
            data = new FileBgen(params);
        }
        else if(params.file_t == FileType::CSV)
        {
            shuffle_csvzstd_to_bin(params.filein, params.fileout, params.buffer, params.scale);
            data = new FileBin(params);
        }
        else
        {
            throw runtime_error("wrong file type used!\n");
        }
        auto te = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(te - ts);
        data->llog << timestamp() << "total elapsed time of permuting data: " << duration.count()
                   << " seconds" << endl;
    }
    else
    {
        if(params.file_t == FileType::PLINK)
        {
            data = new FileBed(params);
        }
        else if(params.file_t == FileType::BGEN)
        {
            data = new FileBgen(params);
        }
        else if(params.file_t == FileType::BEAGLE)
        {
            data = new FileBeagle(params);
        }
        else if(params.file_t == FileType::BINARY)
        {
            data = new FileBin(params);
        }
        else if(params.file_t == FileType::CSV)
        {
            data = new FileCsv(params);
        }
        else
        {
            throw std::invalid_argument(colerror + "invalid input files" + colend);
        }
    }
    // ready for run
    data->prepare(params.blocksize);
    // begin to run
    if(params.svd_t == SvdType::IRAM)
        run_pca_with_arnoldi(data, params);
    else if(params.svd_t == SvdType::PCAoneAlg1 || params.svd_t == SvdType::PCAoneAlg2)
        run_pca_with_halko(data, params);
    else if(params.svd_t == SvdType::FULL)
    {
        data->llog << timestamp() << "running the Full SVD with in-core mode." << endl;
        Eigen::JacobiSVD<MyMatrix> svd(data->G, Eigen::ComputeThinU | Eigen::ComputeThinV);
        data->write_eigs_files(svd.singularValues().head(params.k).array().square() / data->nsnps,
                               svd.matrixU().leftCols(params.k), svd.matrixU().leftCols(params.k));
    }
    else
        throw invalid_argument("unsupported PCA method to apply");
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(t2 - t1).count()
                    * std::chrono::duration<double>::period::num / std::chrono::duration<double>::period::den;
    data->llog << timestamp() << "total elapsed reading time: " << data->readtime << " seconds" << endl;
    data->llog << timestamp() << "total elapsed wall time: " << duration << " seconds" << endl;
    data->llog << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";
    data->llog << commandargs << endl;

    delete data;

    return 0;
}
