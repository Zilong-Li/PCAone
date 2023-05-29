#define _DECLARE_TOOLBOX_HERE

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
    Param params(argc, argv);
    cao.cao.open(params.fileout + ".log");
    string commandargs = params.ss.str();
    cao << commandargs << endl;
    // set number of threads
    // openblas_set_num_threads(params.threads);
    omp_set_num_threads(params.threads);
    Data * data;
    if(params.svd_t == SvdType::PCAoneAlg2 && !params.noshuffle && params.out_of_core)
    {
        tick.clock();
        if(params.file_t == FileType::PLINK)
        {
            permute_plink(params.filein, params.fileout, params.buffer, params.bands);
            data = new FileBed(params);
        }
        else if(params.file_t == FileType::BGEN)
        {
            shuffle_bgen_to_bin(params.filein, params.fileout, params.buffer, params.scale);
            params.file_t = FileType::BINARY;
            data = new FileBin(params);
        }
        else if(params.file_t == FileType::CSV)
        {
            shuffle_csvzstd_to_bin(params.filein, params.fileout, params.buffer, params.scale);
            params.file_t = FileType::BINARY;
            data = new FileBin(params);
        }
        else if(params.file_t == FileType::BINARY)
        {
            data = new FileBin(params);
        }
        else
        {
            throw runtime_error("wrong file type used!\n");
        }
        cao << tick.date() << "total elapsed time of permuting data: " << tick.reltime() << " seconds" << endl;
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
    data->prepare();
    // begin to run
    if(params.svd_t == SvdType::IRAM)
        run_pca_with_arnoldi(data, params);
    else if(params.svd_t == SvdType::PCAoneAlg1 || params.svd_t == SvdType::PCAoneAlg2)
        run_pca_with_halko(data, params);
    else if(params.svd_t == SvdType::FULL)
    {
        cao << tick.date() << "running the Full SVD with in-core mode." << endl;
        if(params.file_t == FileType::PLINK || params.file_t == FileType::BGEN) data->standardize_E();
        Eigen::JacobiSVD<MyMatrix> svd(data->G, Eigen::ComputeThinU | Eigen::ComputeThinV);
        data->write_eigs_files(svd.singularValues().array().square() / data->nsnps, svd.matrixU(),
                               svd.matrixV());
    }
    else
        throw invalid_argument("unsupported PCA method was applied");

    cao << tick.date() << "total elapsed reading time: " << data->readtime << " seconds" << endl;
    cao << tick.date() << "total elapsed wall time: " << tick.abstime() << " seconds" << endl;
    cao << tick.date() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";

    delete data;

    return 0;
}
