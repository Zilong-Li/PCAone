#include "EMU.hpp"
#include "Data.hpp"
// #include "Utils.hpp"
#include "cxxopts.hpp"
// #include "SvdBlas.hpp"

#include <cblas.h>


int main(int argc, char *argv[])
{
    Param params;
    // parse params and check before run
    parse_params(argc, argv, &params);
    // ready for run
    Data data(params);
    // begin to run
    data.run();

    // data.read_bed_whole_genomat();
    // cerr << timestamp() << "begin to init E\n";
    // data.init_E();

    // SvdBlasWide op(data.G);
    // SymEigsSolver< double, LARGEST_ALGE, SvdBlasWide > eigs(&op, k, k * 2 + 1);
    // SvdOp op(data.G);
    // SymEigsSolver< double, LARGEST_ALGE, SvdOp > eigs(&op, params.k, params.k * 2 + 1);
    // eigs.init();
    // cerr << timestamp() << "begin to do eigs\n";
    // eigs.compute(params.maxiter, params.tol);
    // VectorXd evalues;
    // if(eigs.info() == Spectra::SUCCESSFUL)
    //     evalues = eigs.eigenvalues();

    // cout << timestamp() << "Eigenvalues found:\n" << evalues << endl;

    return 0;
}

void parse_params(int argc, char* argv[], struct Param* params)
{
    cxxopts::Options opts(argv[0], "Blockwise Version of EMU");
    opts.add_options("Main")
        ("h,help", "print list of available options")
        ("bfile", "prefix to PLINK .bed/.bim/.fam files", cxxopts::value<std::string>(), "PREFIX")
        ("pfile", "prefix to PLINK2 .pgen/.pvar/.psam files", cxxopts::value<std::string>(), "PREFIX")
        ("bgen", "BGEN file", cxxopts::value<std::string>(), "FILE")
        ("k,eigs", "top k eigs to be calculated", cxxopts::value<int>(),"INT")
        ("n,threads", "number of threads[1]", cxxopts::value<int>(),"INT")
        ("maxiter", "maximum number of arnoldi interations[500]", cxxopts::value<int>(),"INT")
        ("batch", "load all genotypes into RAM", cxxopts::value<bool>()->default_value("false"))
        ;

    try {
        auto vm = opts.parse(argc, argv);
        auto args = vm.arguments();
        // help menu
        if (vm.count("help")){
            cout << opts.help({"Main"}) << "\n\n";
            exit(EXIT_SUCCESS);
        }
        if (vm.count("eigs") != 1) {
            cerr << "ERROR: You must specify the number of top eigs to be calculated.\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("bfile") ) {
            params->intype = "bfile";
            params->bed_prefix = vm["bfile"].as<string>();
        } else if( vm.count("pfile") ) {
            params->intype = "pfile";
            params->pgen_prefix = vm["pfile"].as<string>();
        } else if( vm.count("bgen") ) {
            params->intype = "bgen";
            params->bgen = vm["bgen"].as<string>();
        }
        if (params->intype == "") {
            cerr << "ERROR: You must use either --bfile, --pfile or --bgen.\n";
            exit(EXIT_SUCCESS);
        }
        if( vm.count("eigs") ) params->k = vm["eigs"].as<int>();
        if( vm.count("threads") ) params->threads = vm["threads"].as<int>();
        if( vm.count("maxiter") ) params->maxiter = vm["maxiter"].as<int>();
        if( vm.count("batch") ) params->batch = vm["batch"].as<bool>();
        openblas_set_num_threads(params->threads);

    } catch (const cxxopts::OptionException& e) {
        cerr << "ERROR: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    return;
}
