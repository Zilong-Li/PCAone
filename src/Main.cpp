#include "Arnoldi.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileCsv.hpp"
#include "FilePlink.hpp"
#include "Halko.hpp"
#include "popl/popl.hpp"
#include <omp.h>

#ifdef WITH_OPENBLAS
#include "lapacke.h"
#elif defined WITH_MKL
#include "mkl_lapacke.h"
#endif

using namespace std;
using namespace popl;

string parse_params(int argc, char* argv[], struct Param* params);

int main(int argc, char* argv[])
{
    auto t1 = std::chrono::steady_clock::now();
    Param params;
    // parse params and check before run
    string commandargs = parse_params(argc, argv, &params);
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

string parse_params(int argc, char* argv[], struct Param* params)
{
    string copyr{"PCA All In One (v" + (string)VERSION +
                 ")        https://github.com/Zilong-Li/PCAone\n(C) 2021-2022 Zilong Li        GNU General Public License v3\n\nUsage: PCAone [OPTION]\n\n"};
    OptionParser opts(copyr + "Main options");
    auto help_opt = opts.add<Switch>("", "help", "print list of all options\n");
    opts.add<Switch>("a", "arnoldi", "use IRAM algorithm instead", &params->arnoldi);
    opts.add<Value<string>>("b", "bfile", "prefix of PLINK .bed/.bim/.fam files", "", &params->bed_prefix);
    opts.add<Value<string>>("B", "bgen", "path of BGEN file", "", &params->bgen);
    opts.add<Value<string>>("g", "beagle", "path of beagle file", "", &params->beagle);
    opts.add<Value<string>>("c", "csv", "path of CSV file compressed by zstd", "", &params->csvfile);
    opts.add<Switch>("C", "cpmed", "normalize values by count per median (CPMED) for scRNAs", &params->cpmed);
    opts.add<Switch>("e", "emu", "use EMU algorithm for data with lots of missingness", &params->emu);
    // opts.add<Switch>("f", "fast", "use fast RSVD algorithm with super power iterations", &params->fast);
    opts.add<Switch>("h", "halko", "use normal RSVD algorithm instead", &params->halko);
    opts.add<Value<uint>>("k", "eigs", "top k components to be calculated", params->k, &params->k);
    opts.add<Value<double>>("", "maf", "remove variants with minor allele frequency below maf", params->maf, &params->maf);
    opts.add<Value<double>>("m", "memory", "specify the RAM usage in GB unit", params->memory, &params->memory);
    opts.add<Value<uint>>("n", "threads", "number of threads to use", params->threads, &params->threads);
    opts.add<Value<string>>("o", "out", "prefix of output files", params->outfile, &params->outfile);
    opts.add<Switch>("p", "pcangsd", "use PCAngsd algorithm for genotype likelihood input", &params->pcangsd);
    opts.add<Value<uint>>("P", "maxp", "maximum number of power iteration for RSVD", params->maxp, &params->maxp);
    opts.add<Switch>("V", "printv", "output another eigen vectors with suffix .loadings", &params->printv);
    opts.add<Value<string>>("T", "tmp", "prefix of temporary permuted data", "", &params->tmpfile);
    opts.add<Switch>("v", "verbose", "verbose message output", &params->verbose);
    opts.add<Switch>("", "no-shuffle", "do not shuffle the data if it is already permuted", &params->noshuffle);
    opts.add<Value<uint64>, Attribute::advanced>("M", "", "number of features, eg. SNPs", 0, &params->nsnps);
    opts.add<Value<uint64>, Attribute::advanced>("N", "", "number of samples", 0, &params->nsamples);
    opts.add<Value<uint>, Attribute::advanced>("", "bands", "number of bands to use for fast RSVD", params->bands, &params->bands);
    opts.add<Value<uint>, Attribute::advanced>("", "buffer", "buffer in GB uint used for permuting the data", params->buffer, &params->buffer);
    opts.add<Value<uint>, Attribute::advanced>("", "imaxiter", "maximum number of IRAM interations", params->imaxiter, &params->imaxiter);
    opts.add<Value<double>, Attribute::advanced>("", "itol", "tolerance for IRAM algorithm", params->itol, &params->itol);
    opts.add<Value<uint>, Attribute::advanced>("", "ncv", "number of Lanzcos basis vectors for IRAM", params->ncv, &params->ncv);
    opts.add<Value<uint>, Attribute::advanced>("", "oversamples", "number of oversampling columns for RSVD", params->oversamples, &params->oversamples);
    opts.add<Value<double>, Attribute::advanced>("", "tol", "tolerance for RSVD algorithm", params->tol, &params->tol);
    opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance for EMU/PCAngsd algorithm", params->tolem, &params->tolem);
    opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance for minor allele frequencies estimation update by EM ", params->tolmaf,
                                                 &params->tolmaf);

    std::ostringstream ss;
    // print command line options
    ss << (string) "PCAone (v" + VERSION + ")    https://github.com/Zilong-Li/PCAone\n";
    ss << "Options in effect:\n";
    std::copy(argv, argv + argc, std::ostream_iterator<char*>(ss, " "));
    try
    {
        opts.parse(argc, argv);
        if (params->bed_prefix != "")
        {
            params->intype = FileType::PLINK;
        }
        else if (params->bgen != "")
        {
            params->intype = FileType::BGEN;
        }
        else if (params->beagle != "")
        {
            params->intype = FileType::BEAGLE;
        }
        else if (params->csvfile != "")
        {
            params->intype = FileType::CSV;
        }
        else if (help_opt->count() == 1)
        {
            cout << opts.help(Attribute::advanced) << "\n";
            exit(EXIT_SUCCESS);
        }
        else if (argc == 1)
        {
            cout << opts << "\n";
            exit(EXIT_FAILURE);
        }
        params->ncv = fmax(20, 2 * params->k + 1);
        params->oversamples = fmax(10, params->k);
        // beagle only represents genotype likelihood for pcangsd algorithm now
        if (params->intype == FileType::BEAGLE)
            params->pcangsd = true;
        if (params->emu || params->pcangsd)
        {
            params->runem = true;
        }
        else
        {
            params->maxiter = 0;
        }
        if (params->halko || params->arnoldi)
            params->fast = false;
        if (params->memory > 0)
        {
            params->batch = false;
            if (params->pcangsd)
                throw std::invalid_argument("not support -m option for PCAngsd algorithm yet, but the feature is on the way!");
        }
    }
    catch (const popl::invalid_option& e)
    {
        cerr << "Invalid Option Exception: " << e.what() << "\n";
        cerr << "error:  ";
        if (e.error() == invalid_option::Error::missing_argument)
            cerr << "missing_argument\n";
        else if (e.error() == invalid_option::Error::invalid_argument)
            cerr << "invalid_argument\n";
        else if (e.error() == invalid_option::Error::too_many_arguments)
            cerr << "too_many_arguments\n";
        else if (e.error() == invalid_option::Error::missing_option)
            cerr << "missing_option\n";

        if (e.error() == invalid_option::Error::missing_option)
        {
            string option_name(e.option()->name(OptionName::short_name, true));
            if (option_name.empty())
                option_name = e.option()->name(OptionName::long_name, true);
            cerr << "option: " << option_name << "\n";
        }
        else
        {
            cerr << "option: " << e.option()->name(e.what_name()) << "\n";
            cerr << "value:  " << e.value() << "\n";
        }
        exit(EXIT_FAILURE);
    }
    catch (const std::exception& e)
    {
        cerr << "Exception: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    return ss.str();
}
