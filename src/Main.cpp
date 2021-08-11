#include "cxxopts.hpp"
#include "Data.hpp"
#include "FilePlink.hpp"
#include "FileBeagle.hpp"
#include "Halko.hpp"
#include "Arnoldi.hpp"
// use bgen lib for parsing bgen file
#ifdef WITH_BGEN
#include "FileBgen.hpp"
#endif
// if using external BLAS LAPACK routines
#ifdef WITH_OPENBLAS
#include "lapacke.h"
#elif defined WITH_MKL
#include "mkl_lapacke.h"
#endif
#include <omp.h>

using namespace std;

void parse_params(int argc, char* argv[], struct Param* params);

int main(int argc, char *argv[])
{
    Param params;
    // parse params and check before run
    parse_params(argc, argv, &params);
    // set number of threads
    // openblas_set_num_threads(params.threads);
    omp_set_num_threads(params.threads);
    Data *data;
    if (params.intype == "bfile") {
        data = new FileBed(params);
    #ifdef WITH_BGEN
    } else if( params.intype == "bgen" ) {
        data = new FileBgen(params);
    #endif
    } else if( params.intype == "beagle" ) {
        data = new FileBeagle(params);
    } else {
        exit(EXIT_FAILURE);
    }
    // ready for run
    data->prepare(params.blocksize);
    // begin to run
    if (params.arnoldi)
    {
        run_pca_with_arnoldi(data, params);
    } else {
        run_pca_with_halko(data, params);
    }
    delete data;

    return 0;
}

void parse_params(int argc, char* argv[], struct Param* params)
{
    cxxopts::Options opts(argv[0], "PCA All In One.");
    opts.add_options()
        ("h,help", "print list of main options.")
        ("H,helpall", "print list of all options.")
        ;

    opts.add_options("Main")
        ("k,eigs", "top k components to be calculated.[10]", cxxopts::value<int>(),"INT")
        ("bfile", "prefix to PLINK .bed/.bim/.fam files.", cxxopts::value<std::string>(), "PREFIX")
        // ("pfile", "prefix to PLINK2 .pgen/.pvar/.psam files.", cxxopts::value<std::string>(), "PREFIX")
        #ifdef WITH_BGEN
        ("bgen", "BGEN file.", cxxopts::value<std::string>(), "FILE")
        #endif
        ("beagle", "beagle file.", cxxopts::value<std::string>(), "FILE")
        ("a, arnoldi", "use implicit restarted Arnoldi method instead of default Randomized SVD (Halko).", cxxopts::value<bool>()->default_value("false"))
        ("f, fast", "force to use fast super power iterations for Halko.", cxxopts::value<bool>()->default_value("false"))
        ("emu", "use EMU algorithm for data with large proportion of missingness.", cxxopts::value<bool>()->default_value("false"))
        ("pcangsd", "use PCAngsd algorithm for data with genotype probability.", cxxopts::value<bool>()->default_value("false"))
        ("m,memory", "specify the RAM usage in GB unit instead of exploiting the RAM of the server.", cxxopts::value<double>(),"DOUBLE")
        ("n,threads", "number of threads. [1]", cxxopts::value<int>(),"INT")
        ("o,out", "prefix for output files.", cxxopts::value<string>(),"PREFIX")
        ("v,verbose", "verbose message output.", cxxopts::value<bool>()->default_value("false"))
        ;
    opts.add_options("More")
        ("bands", "number of bands to use for fast Halko.[128]", cxxopts::value<int>(),"INT")
        ("maxp", "maximum number of power iteration for Halko.[20]", cxxopts::value<int>(),"INT")
        ("tol_halko", "tolerance for Halko algorithm. [1e-4]", cxxopts::value<double>(),"DOUBLE")
        ("tol_emu", "tolerance for EMU algorithm. [5e-7]", cxxopts::value<double>(),"DOUBLE")
        ("tol_pcangsd", "tolerance for PCAngsd algorithm. [1e-4]", cxxopts::value<double>(),"DOUBLE")
        ("maxiter", "maximum number of EMU/PCAngsd interations. [100]", cxxopts::value<int>(),"INT")
        ("imaxiter", "maximum number of Arnoldi interations. [500]", cxxopts::value<int>(),"INT")
        ("itol", "tolerance for Arnoldi algorithm. [1e-6]", cxxopts::value<double>(),"DOUBLE")
        ("ncv", "number of Lanzcos basis vectors to use.[max(20, 2*k+1)]", cxxopts::value<int>(),"INT")
        // ("tolmaf", "MAF tolerance for PCAngsd algorithm. [1e-5]", cxxopts::value<double>(),"DOUBLE")
        ;

    try {
        auto vm = opts.parse(argc, argv);
        auto args = vm.arguments();
        // help menu
        if (vm.count("help")){
            cout << opts.help({"", "Main"}) << "\n\n";
            exit(EXIT_SUCCESS);
        }
        if (vm.count("helpall")){
            cout << opts.help({"", "Main", "More"}) << "\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("out") ) params->outfile = vm["out"].as<string>();
        if( vm.count("eigs") ) {params->k = vm["eigs"].as<int>(); params->ncv = fmax(20, 2 * params->k + 1);}
        if( vm.count("threads") ) params->threads = vm["threads"].as<int>();
        if( vm.count("tol_halko") ) params->tol_halko = vm["tol_halko"].as<double>();
        if( vm.count("imaxiter") ) params->imaxiter = vm["imaxiter"].as<int>();
        if( vm.count("maxp") ) params->p = vm["maxp"].as<int>();
        if( vm.count("bands") ) params->bands = vm["bands"].as<int>();
        if( vm.count("ncv") ) params->ncv = vm["ncv"].as<int>();
        if( vm.count("itol") ) params->itol = vm["itol"].as<double>();
        if( vm.count("arnoldi") ) params->arnoldi = vm["arnoldi"].as<bool>();
        if( vm.count("verbose") ) params->verbose = vm["verbose"].as<bool>();
        if( vm.count("pcangsd") ) params->pcangsd = vm["pcangsd"].as<bool>();
        if( vm.count("emu") ) params->emu = vm["emu"].as<bool>();
        if( vm.count("fast") ) params->fast = vm["fast"].as<bool>();
        if (params->emu || params->pcangsd) {
            if( vm.count("maxiter") ) params->maxiter = vm["maxiter"].as<int>();
            if( vm.count("tol_pcangsd") ) params->tol_pcangsd = vm["tol_pcangsd"].as<double>();
            if( vm.count("tolmaf") ) params->tolmaf = vm["tolmaf"].as<double>();
            if( vm.count("tol_emu") ) params->tol_emu = vm["tol_emu"].as<double>();
            if( params->pcangsd ) params->tol = params->tol_pcangsd;
            if( params->emu ) params->tol = params->tol_emu;
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
        #ifdef WITH_BGEN
        } else if( vm.count("bgen") ) {
            params->intype = "bgen";
            params->bgen = vm["bgen"].as<string>();
        #endif
        } else if( vm.count("beagle") ) {
            params->intype = "beagle";
            params->beagle = vm["beagle"].as<string>();
            params->pcangsd = true;
        } else {
            cerr << "ERROR: please check the options using -h.\n";
            exit(EXIT_SUCCESS);
        }

        if (vm.count("out") != 1) {
            cerr << "ERROR: You must specify the output prefix.\n";
            exit(EXIT_SUCCESS);
        }
        if (vm.count("eigs") != 1) {
            cerr << "ERROR: You must specify the number of top eigs to be calculated.\n";
            exit(EXIT_SUCCESS);
        }

    } catch (const cxxopts::OptionException& e) {
        cerr << "ERROR: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    return;
}
