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
    data->prepare();
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
        ("helpall", "print list of all options.")
        ;

    opts.add_options("Main")
        ("k,eigs", "top k components to be calculated.", cxxopts::value<int>(),"INT")
        ("bfile", "prefix to PLINK .bed/.bim/.fam files.", cxxopts::value<std::string>(), "PREFIX")
        ("pfile", "prefix to PLINK2 .pgen/.pvar/.psam files.", cxxopts::value<std::string>(), "PREFIX")
        #ifdef WITH_BGEN
        ("bgen", "BGEN file.", cxxopts::value<std::string>(), "FILE")
        #endif
        ("beagle", "beagle file.", cxxopts::value<std::string>(), "FILE")
        ("b,blocksize", "size of block and in number of SNPs.[0]", cxxopts::value<int>(),"INT")
        ("emu", "using EMU algorithm for data with missingness.", cxxopts::value<bool>()->default_value("false"))
        ("pcangsd", "using PCAngsd algorithm for data with genotype probability.", cxxopts::value<bool>()->default_value("false"))
        ("n,threads", "number of threads. [1]", cxxopts::value<int>(),"INT")
        ("p,power", "number of power iterations for Halko method.[2]", cxxopts::value<int>(),"INT")
        ("A, arnoldi", "using implicit restarted Arnoldi method instead of default Randomized SVD (Halko)", cxxopts::value<bool>()->default_value("false"))
        ("o,out", "prefix for output files.", cxxopts::value<string>(),"PREFIX")
        ("v,verbose", "verbose message output.", cxxopts::value<bool>()->default_value("false"))
        ;
    opts.add_options("More")
        ("ncv", "number of Lanzcos basis vectors to use.[max(20, 2*k+1)]", cxxopts::value<int>(),"INT")
        ("imaxiter", "maximum number of Arnoldi interations. [1000]", cxxopts::value<int>(),"INT")
        ("itol", "tolerance for Arnoldi algorithm. [1e-6]", cxxopts::value<double>(),"DOUBLE")
        ("true", "true file with eigen vectors.", cxxopts::value<std::string>(), "FILE")
        ("maxiter", "maximum number of EMU or PCAngsd interations. [100]", cxxopts::value<int>(),"INT")
        ("tol", "tolerance for EMU or PCAngsd algorithm. [5e-7/1e-4]", cxxopts::value<double>(),"DOUBLE")
        ("tolmaf", "MAF tolerance for PCAngsd algorithm. [1e-5]", cxxopts::value<double>(),"DOUBLE")
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
        if( vm.count("true") ) params->truefile = vm["true"].as<string>();
        if( vm.count("out") ) params->outfile = vm["out"].as<string>();
        if( vm.count("eigs") ) {params->k = vm["eigs"].as<int>(); params->ncv = fmax(20, 2 * params->k + 1);}
        if( vm.count("power") ) params->p = vm["power"].as<int>();
        if( vm.count("threads") ) params->threads = vm["threads"].as<int>();
        if( vm.count("imaxiter") ) params->imaxiter = vm["imaxiter"].as<int>();
        if( vm.count("ncv") ) params->ncv = vm["ncv"].as<int>();
        if( vm.count("itol") ) params->itol = vm["itol"].as<double>();
        if( vm.count("arnoldi") ) params->arnoldi = vm["arnoldi"].as<bool>();
        if( vm.count("verbose") ) params->verbose = vm["verbose"].as<bool>();
        if( vm.count("pcangsd") ) params->pcangsd = vm["pcangsd"].as<bool>();
        if( vm.count("emu") ) params->emu = vm["emu"].as<bool>();
        if (params->emu || params->pcangsd) {
            if( vm.count("maxiter") ) params->maxiter = vm["maxiter"].as<int>();
            if( params->pcangsd ) params->tol = params->tol2;
            if( vm.count("tolmaf") ) params->tolmaf = vm["tolmaf"].as<double>();
            if( vm.count("tol") ) params->tol = vm["tol"].as<double>();
        } else {
            params->maxiter = 0;
        }
        if( vm.count("blocksize") ) {
            params->blocksize = vm["blocksize"].as<int>();
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
            cerr << "ERROR: You must use either --beagle, --bfile, --pfile or --bgen.\n";
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
