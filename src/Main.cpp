#include "Arnoldi.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileCsv.hpp"
#include "FilePlink.hpp"
#include "Halko.hpp"
#include "Cmd.hpp"
#include <omp.h>

#ifdef WITH_OPENBLAS
#include "lapacke.h"
#elif defined WITH_MKL
#include "mkl_lapacke.h"
#endif

using namespace std;
using namespace popl;

int main(int argc, char* argv[])
{
    auto t1 = std::chrono::steady_clock::now();
    Param params(argc, argv);
    string commandargs = params.ss.str();
    // parse params and check before run
    // string commandargs = parse_params(argc, argv, &params);
    // set number of threads
    // openblas_set_num_threads(params.threads);
    omp_set_num_threads(params.threads);
    Data* data;
    if (params.intype == FileType::PLINK)
    {
        if (!params.batch && params.fast)
        {
            if (params.noshuffle)
            {
                cout << timestamp() << "warning: running fast fancy RSVD without shuffling the data!" << endl;
            }
            else
            {
                auto ts = std::chrono::steady_clock::now();
                string fout = params.outfile + ".perm";
                if (params.tmpfile != "")
                    fout = params.tmpfile;
                permute_plink(params.bed_prefix, fout, params.buffer);
                auto te = std::chrono::steady_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::seconds>(te - ts);
                cout << timestamp() << "total elapsed time of permuting data: " << duration.count() << " seconds" << endl;
            }
        }
        data = new FileBed(params);
    }
    else if (params.intype == FileType::BGEN)
    {
        data = new FileBgen(params);
    }
    else if (params.intype == FileType::BEAGLE)
    {
        data = new FileBeagle(params);
    }
    else if (params.intype == FileType::CSV)
    {
        data = new FileCsv(params);
    }
    else
    {
        throw std::invalid_argument(
            "\nplease specify the input file using one of --bfile, --bgen, --beagle, --csv option!\nusing --help to show all options.\n");
    }
    // ready for run
    data->prepare(params.blocksize);
    // begin to run
    if (params.arnoldi)
    {
        run_pca_with_arnoldi(data, params);
    }
    else
    {
        run_pca_with_halko(data, params);
    }
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(t2 - t1).count() * std::chrono::duration<double>::period::num / std::chrono::duration<double>::period::den;
    data->llog << timestamp() << "total elapsed reading time: " << data->readtime << " seconds" << endl;
    data->llog << timestamp() << "total elapsed wall time: " << duration << " seconds" << endl;
    data->llog << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";
    data->llog << commandargs << endl;

    delete data;

    return 0;
}

