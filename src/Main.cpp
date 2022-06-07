#include "Data.hpp"
#include "FilePlink.hpp"
#include "FileBeagle.hpp"
#include "FileBgen.hpp"
#include "FileCsv.hpp"
#include "Halko.hpp"
#include "Arnoldi.hpp"
#include "cxxopts.hpp"
#include <omp.h>

#ifdef WITH_OPENBLAS
#include "lapacke.h"
#elif defined WITH_MKL
#include "mkl_lapacke.h"
#endif

using namespace std;

string parse_params(int argc, char* argv[], struct Param* params);

int main(int argc, char *argv[])
{
    auto t1 = std::chrono::steady_clock::now();
    Param params;
    // parse params and check before run
    string commandargs = parse_params(argc, argv, &params);
    // set number of threads
    // openblas_set_num_threads(params.threads);
    omp_set_num_threads(params.threads);
    Data *data;
    if (params.intype == "bfile") {
        if (!params.batch && params.fast) {
            if (params.shuffle) {
                auto ts = std::chrono::steady_clock::now();
                string fout = params.outfile + ".perm";
                if (params.tmpfile != "")
                  fout = params.tmpfile;
                permute_plink(params.bed_prefix, fout, params.buffer);
                auto te = std::chrono::steady_clock::now();
                auto duration =
                    std::chrono::duration_cast<std::chrono::seconds>(te - ts);
                cout << timestamp() << "total elapsed time of permuting data: "
                     << duration.count() << " seconds" << endl;
            } else {
                cout << timestamp() << "warning: running fast fancy RSVD without shuffling the "
                      "data!" << endl;
            }
        }
        data = new FileBed(params);
    } else if( params.intype == "bgen" ) {
        data = new FileBgen(params);
    } else if( params.intype == "beagle" ) {
        data = new FileBeagle(params);
    } else if( params.intype == "csv" ) {
        data = new FileCsv(params);
    } else {
        exit(EXIT_FAILURE);
    }
    // start logging
    data->llog.clog.open(string(params.outfile + ".log").c_str(), ios::out | ios::trunc);
    // ready for run
    data->prepare(params.blocksize);
    // begin to run
    if (params.arnoldi)
    {
        run_pca_with_arnoldi(data, params);
    } else {
        run_pca_with_halko(data, params);
    }
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    data->llog << timestamp() << "total elapsed reading time: " << data->readtime << " seconds" << endl;
    data->llog << timestamp() << "total elapsed wall time: " << duration.count() << " seconds" << endl;
    data->llog << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";
    data->llog << commandargs << endl;

    delete data;

    return 0;
}

string parse_params(int argc, char* argv[], struct Param* params)
{
    cxxopts::Options opts(argv[0], (string)"PCA All In One (v" + VERSION + ")        https://github.com/Zilong-Li/PCAone\n(C) 2021-2022 Zilong Li        GNU General Public License v3");
    opts.add_options()
        ("help", "Print list of all options")
        ;

    opts.add_options("Main")
        ("a, arnoldi", "use IRAM algorithm", cxxopts::value<bool>()->default_value("false"))
        ("beagle", "path of beagle file", cxxopts::value<std::string>())
        ("bfile", "prefix of PLINK .bed/.bim/.fam files", cxxopts::value<std::string>())
        // ("pfile", "prefix to PLINK2 .pgen/.pvar/.psam files.", cxxopts::value<std::string>(), "PREFIX")
        ("bgen", "path of BGEN file", cxxopts::value<std::string>())
        ("csv", "path of zstd compressed csv file", cxxopts::value<std::string>())
        ("cpmed", "normalize values by count per median (CPMED) for scRNAs", cxxopts::value<bool>()->default_value("false"))
        ("e,emu", "use EMU algorithm for data with lots of missingness", cxxopts::value<bool>()->default_value("false"))
        ("f, fast", "use fast RSVD algorithm with super power iterations", cxxopts::value<bool>()->default_value("true"))
        ("h, halko", "use normal RSVD algorithm", cxxopts::value<bool>()->default_value("false"))
        ("k,eigs", "top k components to be calculated", cxxopts::value<int>()->default_value("10"))
        ("maxp", "maximum number of power iteration for RSVD", cxxopts::value<int>()->default_value("20"))
        ("m,memory", "specify the RAM usage in GB unit", cxxopts::value<double>())
        ("n,threads", "number of threads", cxxopts::value<int>()->default_value("10"))
        ("no-shuffle", "permute data by features for fast RSVD.", cxxopts::value<bool>()->default_value("false"))
        ("o,out", "prefix of output files", cxxopts::value<string>())
        ("p,pcangsd", "use PCAngsd algorithm for genotype likelihood input", cxxopts::value<bool>()->default_value("false"))
        ("printv", "output another eigen vectors with suffix .loadings", cxxopts::value<bool>()->default_value("false"))
        ("tmp", "prefix of permuted data", cxxopts::value<string>())
        ("v,verbose", "verbose message output", cxxopts::value<bool>()->default_value("false"))
        ("M", "number of features, eg. SNPs", cxxopts::value<int>())
        ("N", "number of samples", cxxopts::value<int>())
        ;
    opts.add_options("More")
        ("bands", "number of bands to use for fast RSVD", cxxopts::value<int>()->default_value("64"))
        ("buffer", "buffer in GB uint used for permuting the data", cxxopts::value<int>()->default_value("2"))
        ("imaxiter", "maximum number of IRAM interations", cxxopts::value<int>()->default_value("500"))
        ("itol", "tolerance for IRAM algorithm", cxxopts::value<double>()->default_value("1e-6"))
        ("maxiter", "maximum number of EMU/PCAngsd interations", cxxopts::value<int>()->default_value("100"))
        ("ncv", "number of Lanzcos basis vectors", cxxopts::value<string>()->default_value("max(20, 2k+1)"))
        ("oversamples", "number of oversampling columns for RSVD", cxxopts::value<string>()->default_value("max(k, 10)"))
        ("tol-em", "tolerance for EMU/PCAngsd algorithm", cxxopts::value<double>()->default_value("1e-4"))
        ("tol-rsvd", "tolerance for RSVD algorithm", cxxopts::value<double>()->default_value("1e-4"))
        ("tol-maf", "MAF tolerance for PCAngsd algorithm", cxxopts::value<double>()->default_value("1e-4"))
        ;

    std::ostringstream ss;
    // print command line options
    ss << (string)"PCAone (v" + VERSION + ")    https://github.com/Zilong-Li/PCAone\n";
    ss << "Options in effect:\n";
    try {
        auto vm = opts.parse(argc, argv);
        auto args = vm.arguments();
        // help menu
        if (vm.count("help")){
            cout << opts.help({"", "Main", "More"}) << "\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("eigs") ) {
            params->k = vm["eigs"].as<int>();
            params->ncv = fmax(20, 2 * params->k + 1);
            params->oversamples = fmax(10, params->k);
        }
        if( vm.count("threads") ) params->threads = vm["threads"].as<int>();
        if( vm.count("N") ) params->nsamples = vm["N"].as<int>();
        if( vm.count("M") ) params->nsnps = vm["M"].as<int>();
        if( vm.count("tol-rsvd") ) params->tol_halko = vm["tol-rsvd"].as<double>();
        if( vm.count("imaxiter") ) params->imaxiter = vm["imaxiter"].as<int>();
        if( vm.count("maxp") ) params->maxp = vm["maxp"].as<int>();
        if( vm.count("oversamples") ) params->oversamples = vm["oversamples"].as<int>();
        if( vm.count("bands") ) params->bands = vm["bands"].as<int>();
        if( vm.count("buffer") ) params->buffer = vm["buffer"].as<int>();
        if( vm.count("ncv") ) params->ncv = vm["ncv"].as<int>();
        if( vm.count("itol") ) params->itol = vm["itol"].as<double>();
        if( vm.count("verbose") ) params->verbose = vm["verbose"].as<bool>();
        if( vm.count("printv") ) params->printv = vm["printv"].as<bool>();
        if( vm.count("pcangsd") ) params->pcangsd = vm["pcangsd"].as<bool>();
        if( vm.count("emu") ) params->emu = vm["emu"].as<bool>();
        if( vm.count("no-shuffle") ) params->shuffle = (vm["no-shuffle"].as<bool>() == 0);
        if( vm.count("cpmed") ) params->cpmed = vm["cpmed"].as<bool>();
        if (params->emu || params->pcangsd) {
            params->runem = true;
            if( vm.count("maxiter") ) params->maxiter = vm["maxiter"].as<int>();
            if( vm.count("tol-maf") ) params->tolmaf = vm["tol-maf"].as<double>();
            if( vm.count("tol-em") ) params->tol = vm["tol-em"].as<double>();
        } else {
            params->maxiter = 0;
        }
        if( vm.count("memory") ) {
            params->memory = vm["memory"].as<double>();
            params->batch = false;
        }
        if( vm.count("bfile") ) {
            params->intype = "bfile";
            params->bed_prefix = vm["bfile"].as<string>();
            params->pcangsd = false;
        } else if( vm.count("pfile") ) {
            params->intype = "pfile";
            params->pgen_prefix = vm["pfile"].as<string>();
        } else if( vm.count("bgen") ) {
            params->intype = "bgen";
            params->bgen = vm["bgen"].as<string>();
        } else if( vm.count("beagle") ) {
            params->intype = "beagle";
            params->beagle = vm["beagle"].as<string>();
            params->pcangsd = true;
        } else if( vm.count("csv") ) {
            params->intype = "csv";
            params->csvfile = vm["csv"].as<string>();
        } else {
            cout << opts.help({"", "Main"}) << "\n\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("arnoldi") ) {
            params->arnoldi = vm["arnoldi"].as<bool>();
        } else if( vm.count("halko") ) {
            params->halko = vm["halko"].as<bool>();
        } else {
            params->fast = true;
        }
        if( vm.count("tmp") ) params->tmpfile = vm["tmp"].as<string>();
        if( vm.count("out") ) params->outfile = vm["out"].as<string>();

        for (size_t counter = 1; counter < args.size(); counter++) {
            ss << "  --" << args[counter - 1].key() << " ";
            if (args[counter - 1].value() == "true") {
                ss << "\\\n";
                continue;
            }
            ss << args[counter - 1].value() << " \\" << endl;
        }
        // last option (skip \ at the end)
        ss << "  --" << args.back().key() << " ";
        if (args.back().value() != "true") ss << args.back().value();

    } catch (const cxxopts::OptionException& e) {
        throw e.what();
    }

    return ss.str();
}
