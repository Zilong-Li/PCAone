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
    if(!params.noshuffle && !params.batch)
    {
        params.binfile = params.outfile + ".perm";
        if(params.tmpfile != "") params.binfile = params.tmpfile;
        if(params.intype == FileType::PLINK)
        {
            auto ts = std::chrono::steady_clock::now();
            permute_plink(params.bed_prefix, params.binfile, params.buffer, params.bands);
            auto te = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(te - ts);
            data = new FileBed(params);
            data->llog << timestamp() << "total elapsed time of permuting data: " << duration.count()
                       << " seconds" << endl;
        }
        else if(params.intype == FileType::CSV)
        {
            shuffle_csvzstd_to_bin(params.csvfile, params.binfile, params.buffer, params.cpmed);
        }
        data = new FileBin(params);
    }
    else
    {
        if(params.intype == FileType::PLINK)
        {
            data = new FileBed(params);
            if(!params.batch && params.fast)
                data->llog << timestamp()
                           << colwarn
                                  + "running PCAone (algorithm2) without shuffuling the input data. make "
                                    "sure it's permuted."
                                  + colend
                           << endl;
        }
        else if(params.intype == FileType::BGEN)
        {
            data = new FileBgen(params);
        }
        else if(params.intype == FileType::BEAGLE)
        {
            data = new FileBeagle(params);
        }
        else if(params.intype == FileType::CSV)
        {
            data = new FileCsv(params);
        }
        else
        {
            throw std::invalid_argument(colerror
                                        + "\nplease specify the input file using one of --bfile, --bgen, "
                                          "--beagle, --csv option!\nusing --help to show all options.\n");
        }
    }
    // ready for run
    data->prepare(params.blocksize);
    // begin to run
    if(params.arnoldi)
    {
        run_pca_with_arnoldi(data, params);
    }
    else
    {
        run_pca_with_halko(data, params);
    }
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
