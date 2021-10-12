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
        if (!params.batch && params.fast && !params.noshuffle)
            permute_plink(params.bed_prefix, params.buffer);
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
    cout << timestamp() << "eigenvecs and eigenvals are saved. have a nice day. bye!\n";

    delete data;

    return 0;
}

void parse_params(int argc, char* argv[], struct Param* params)
{
    cxxopts::Options opts(argv[0], (string)"PCA All In One (v" + VERSION + ")");
    opts.add_options()
        ("h,help", "print list of main options.")
        ("H", "print list of all options.")
        ;

    opts.add_options("Main")
        ("a, arnoldi", "use Implicit Restarted Arnoldi method instead of default Randomized SVD(Halko).", cxxopts::value<bool>()->default_value("false"))
        ("beagle", "beagle file.", cxxopts::value<std::string>(), "FILE")
        ("bfile", "prefix of PLINK .bed/.bim/.fam files.", cxxopts::value<std::string>(), "PREFIX")
        // ("pfile", "prefix to PLINK2 .pgen/.pvar/.psam files.", cxxopts::value<std::string>(), "PREFIX")
        #ifdef WITH_BGEN
        ("bgen", "BGEN file.", cxxopts::value<std::string>(), "FILE")
        #endif
        ("e,emu", "use EMU algorithm for data with large proportion of missingness.", cxxopts::value<bool>()->default_value("false"))
        ("f, fast", "force to use fast super power iterations for Halko.", cxxopts::value<bool>()->default_value("false"))
        ("k,eigs", "top k components to be calculated.[10]", cxxopts::value<int>(),"INT")
        ("m,memory", "specify the RAM usage in GB unit instead of exploiting the RAM of the server.", cxxopts::value<double>(),"DOUBLE")
        ("n,threads", "number of threads.[1]", cxxopts::value<int>(),"INT")
        ("o,out", "prefix of output files.", cxxopts::value<string>(),"PREFIX")
        ("p,pcangsd", "use PCAngsd algorithm for data with genotype probability.", cxxopts::value<bool>()->default_value("false"))
        ("v,verbose", "verbose message output.", cxxopts::value<bool>()->default_value("false"))
        ;
    opts.add_options("More")
        ("bands", "number of bands to use for fast Halko.[64]", cxxopts::value<int>(),"INT")
        ("buffer", "number of snps as a buffer for permuting the data.[1]", cxxopts::value<int>(),"INT")
        ("imaxiter", "maximum number of Arnoldi interations.[500]", cxxopts::value<int>(),"INT")
        ("itol", "tolerance for Arnoldi algorithm.[1e-6]", cxxopts::value<double>(),"DOUBLE")
        ("maxp", "maximum number of power iteration for Halko.[20]", cxxopts::value<int>(),"INT")
        ("maxiter", "maximum number of EMU/PCAngsd interations.[100]", cxxopts::value<int>(),"INT")
        ("ncv", "number of Lanzcos basis vectors.[max(20, 2*k+1)]", cxxopts::value<int>(),"INT")
        ("no-shuffle", "do not shuffle the matrix for fast Halko blocksize mode.", cxxopts::value<bool>()->default_value("false"))
        ("oversamples", "the number of oversampling columns for Halko.[10]", cxxopts::value<int>(),"INT")
        ("tol-em", "tolerance for EMU/PCAngsd algorithm.[1e-5]", cxxopts::value<double>(),"DOUBLE")
        ("tol-halko", "tolerance for Halko algorithm.[1e-4]", cxxopts::value<double>(),"DOUBLE")
        ("tol-maf", "MAF tolerance for PCAngsd algorithm.[1e-4]", cxxopts::value<double>(),"DOUBLE")
        ;

    try {
        auto vm = opts.parse(argc, argv);
        auto args = vm.arguments();
        // help menu
        if (vm.count("help")){
            cout << opts.help({"", "Main"}) << "\n\n";
            exit(EXIT_SUCCESS);
        }
        if (vm.count("H")){
            cout << opts.help({"", "Main", "More"}) << "\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("out") ) params->outfile = vm["out"].as<string>();
        if( vm.count("eigs") ) {params->k = vm["eigs"].as<int>(); params->ncv = fmax(20, 2 * params->k + 1);}
        if( vm.count("threads") ) params->threads = vm["threads"].as<int>();
        if( vm.count("tol-halko") ) params->tol_halko = vm["tol-halko"].as<double>();
        if( vm.count("imaxiter") ) params->imaxiter = vm["imaxiter"].as<int>();
        if( vm.count("maxp") ) params->p = vm["maxp"].as<int>();
        if( vm.count("oversamples") ) params->oversamples = vm["oversamples"].as<int>();
        if( vm.count("bands") ) params->bands = vm["bands"].as<int>();
        if( vm.count("buffer") ) params->buffer = vm["buffer"].as<int>();
        if( vm.count("ncv") ) params->ncv = vm["ncv"].as<int>();
        if( vm.count("itol") ) params->itol = vm["itol"].as<double>();
        if( vm.count("arnoldi") ) params->arnoldi = vm["arnoldi"].as<bool>();
        if( vm.count("verbose") ) params->verbose = vm["verbose"].as<bool>();
        if( vm.count("pcangsd") ) params->pcangsd = vm["pcangsd"].as<bool>();
        if( vm.count("emu") ) params->emu = vm["emu"].as<bool>();
        if( vm.count("fast") ) params->fast = vm["fast"].as<bool>();
        if( vm.count("no-shuffle") ) params->noshuffle = vm["no-shuffle"].as<bool>();
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
            cout << opts.help({"", "Main"}) << "\n\n";
            exit(EXIT_SUCCESS);
        }

        if (vm.count("out") != 1) {
            throw std::invalid_argument("ERROR: You must specify the output prefix.\n");
        }

    } catch (const cxxopts::OptionException& e) {
        throw e.what();
    }

    return;
}
