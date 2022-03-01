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
        if (!params.batch && params.fast && params.shuffle)
            permute_plink2(params.bed_prefix, params.buffer);
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
    cxxopts::Options opts(argv[0], (string)"PCA All In One (v" + VERSION + ")        https://github.com/Zilong-Li/PCAone\n(C) 2021-2022 Zilong Li        GNU General Public License v3");
    opts.add_options()
        ("help", "Print list of main options.")
        ("helpall", "Print list of all options.")
        ;

    opts.add_options("Main")
        ("beagle", "path of beagle file.", cxxopts::value<std::string>(), "FILE")
        ("bfile", "prefix of PLINK .bed/.bim/.fam files.", cxxopts::value<std::string>(), "PREFIX")
        // ("pfile", "prefix to PLINK2 .pgen/.pvar/.psam files.", cxxopts::value<std::string>(), "PREFIX")
        ("bgen", "path of BGEN file.", cxxopts::value<std::string>(), "FILE")
        ("csv", "path of zstd compressed csv file.", cxxopts::value<std::string>(), "FILE")
        ("cpmed", "normalize values by count per median (CPMED)for scRNAs", cxxopts::value<bool>()->default_value("false"))
        ("maxp", "maximum number of power iteration for Halko.[20]", cxxopts::value<int>(),"INT")
        ("printv", "print out another eigen vectors with suffix .loadings.", cxxopts::value<bool>()->default_value("false"))
        ("e,emu", "use EMU algorithm for data with lots of missingness.", cxxopts::value<bool>()->default_value("false"))
        ("f, fast", "force to use fast super power iterations for Halko.", cxxopts::value<bool>()->default_value("false"))
        ("h, halko", "use Halko method instead of default Arnoldi method.", cxxopts::value<bool>()->default_value("false"))
        ("k,eigs", "top k components to be calculated.[10]", cxxopts::value<int>(),"INT")
        ("m,memory", "specify the RAM usage in GB unit instead of exploiting the RAM of the server.", cxxopts::value<double>(),"DOUBLE")
        ("n,threads", "number of threads.[1]", cxxopts::value<int>(),"INT")
        ("o,out", "prefix of output files.", cxxopts::value<string>(),"PREFIX")
        ("p,pcangsd", "use PCAngsd algorithm for genotype likelihood input.", cxxopts::value<bool>()->default_value("false"))
        ("v,verbose", "verbose message output.", cxxopts::value<bool>()->default_value("false"))
        ("M", "number of features. eg. SNPs.", cxxopts::value<int>(),"INT")
        ("N", "number of samples.", cxxopts::value<int>(),"INT")
        ;
    opts.add_options("More")
        ("bands", "number of bands to use for fast Halko.[32]", cxxopts::value<int>(),"INT")
        ("buffer", "buffer in GB uint used for permuting the data.[2]", cxxopts::value<int>(),"INT")
        ("imaxiter", "maximum number of Arnoldi interations.[500]", cxxopts::value<int>(),"INT")
        ("itol", "tolerance for Arnoldi algorithm.[1e-6]", cxxopts::value<double>(),"DOUBLE")
        ("maxiter", "maximum number of EMU/PCAngsd interations.[100]", cxxopts::value<int>(),"INT")
        ("ncv", "number of Lanzcos basis vectors.[max(20, 2*k+1)]", cxxopts::value<int>(),"INT")
        ("shuffle", "permute data by features for fast Halko.", cxxopts::value<bool>()->default_value("false"))
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
        if (vm.count("helpall")){
            cout << opts.help({"", "Main", "More"}) << "\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("out") ) params->outfile = vm["out"].as<string>();
        if( vm.count("eigs") ) {params->k = vm["eigs"].as<int>(); params->ncv = fmax(20, 2 * params->k + 1);}
        if( vm.count("threads") ) params->threads = vm["threads"].as<int>();
        if( vm.count("N") ) params->nsamples = vm["N"].as<int>();
        if( vm.count("M") ) params->nsnps = vm["M"].as<int>();
        if( vm.count("tol-halko") ) params->tol_halko = vm["tol-halko"].as<double>();
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
        if( vm.count("shuffle") ) params->shuffle = vm["shuffle"].as<bool>();
        if( vm.count("cpmed") ) params->cpmed = vm["cpmed"].as<bool>();
        if( vm.count("halko") ) {
            params->halko = vm["halko"].as<bool>();
            params->arnoldi = false;
        } else if( vm.count("fast") ) {
            params->fast = vm["fast"].as<bool>();
            params->arnoldi = false;
        }
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

        if (vm.count("out") != 1) {
            throw std::invalid_argument("ERROR: You must specify the output prefix.\n");
        }

    } catch (const cxxopts::OptionException& e) {
        throw e.what();
    }

    return;
}
