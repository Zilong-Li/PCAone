#include "Cmd.hpp"

#include "popl/popl.hpp"

using namespace popl;

Param::Param(int argc, char ** argv)
{
    // clang-format off
    std::string copyr{"PCA All In One (v" + (std::string)VERSION + ")        https://github.com/Zilong-Li/PCAone\n" +
                      "(C) 2021-2022 Zilong Li        GNU General Public License v3\n\n" +
    "\x1B[32m" +
                      "Usage: use plink files as input and apply the window-based RSVD method\n" +
                      "       PCAone --bfile plink -m 2 -n 20 -k 10\n\n" +
                      "       use csv file as input and apply the Implicitly Restarted Arnoldi Method\n" +
                      "       PCAone --csv csv.zst --svd 0 -m 2 -n 20 -k 10\n\n" +
    "\033[0m"};
    OptionParser opts(copyr + "Main options");
    auto help_opt = opts.add<Switch>("h", "help", "print list of all options including hidden advanced options");
    auto svd_opt = opts.add<Value<uint>>("d", "svd", "svd method to be applied. default 2 is recommended for big data.\n"
                                         "0: the Implicitly Restarted Arnoldi Method\n"
                                         "1: the Yu's single-pass Randomized SVD with power iterations\n"
                                         "2: the proposed window-based Randomized SVD  method\n"
                                         "3: the full Singular Value Decomposition.", 2);
    auto plinkfile = opts.add<Value<std::string>>("b", "bfile", "prefix to PLINK .bed/.bim/.fam files", "", &filein);
    auto binfile = opts.add<Value<std::string>>("B", "binary", "path of binary file", "", &filein);
    auto csvfile = opts.add<Value<std::string>>("c", "csv", "path of comma seperated CSV file compressed by zstd", "", &filein);
    auto bgenfile = opts.add<Value<std::string>>("g", "bgen", "path of BGEN file", "", &filein);
    auto beaglefile = opts.add<Value<std::string>>("G", "beagle", "path of BEAGLE file", "", &filein);
    opts.add<Value<uint>>("k", "pc", "top k components to be calculated", k, &k);
    opts.add<Value<double>>("m", "memory", "specify the RAM usage in GB unit. default [0] uses all RAM", memory, &memory);
    opts.add<Value<uint>>("n", "threads", "number of threads for multithreading", threads, &threads);
    opts.add<Value<std::string>>("o", "out", "prefix to output files. default [pcaone]", fileout, &fileout);
    opts.add<Value<uint>>("p", "maxp", "maximum number of power iterations for RSVD algorithm", maxp, &maxp);
    opts.add<Switch>("S", "no-shuffle", "do not shuffle the data if it is already permuted", &noshuffle);
    opts.add<Switch>("v", "verbose", "verbose message output", &verbose);
    opts.add<Value<uint>>("w", "batches", "number of mini-batches to be used by PCAone (algorithm2)", bands, &bands);
    opts.add<Value<uint>>("C", "scale", "do scaling for input file.\n"
                          "0: do just centering\n"
                          "1: do log transformation eg. log(x+0.01) for RNA-seq data\n"
                          "2: do count per median log transformation (CPMED) for scRNAs",
                          scale,  &scale);
    opts.add<Switch>("", "emu", "use EMU algorithm for genotype data with missingness", &emu);
    opts.add<Switch>("", "pcangsd", "use PCAngsd algorithm for genotype likelihood input", &pcangsd);
    opts.add<Switch>("", "ld", "estimate ld for admixed population", &ld);
    opts.add<Value<uint>>("", "ld-window", "ld window size in base units instead of number of sites", ld_window_bp, &ld_window_bp);
    opts.add<Value<double>>("", "maf", "skip variants with minor allele frequency below maf", maf, &maf);
    opts.add<Switch>("U", "printu", "output eigen vector of each epoch (for tests)", &printu);
    opts.add<Switch>("V", "printv", "output the right eigen vectors with suffix .loadings", &printv);
    opts.add<Value<uint>, Attribute::advanced>("", "M", "number of features (eg. SNPs) if already known", 0, &nsnps);
    opts.add<Value<uint>, Attribute::advanced>("", "N", "number of samples if already known", 0, &nsamples);
    opts.add<Switch, Attribute::advanced>("", "haploid", "the plink format represents haploid data", &haploid);
    opts.add<Value<uint>, Attribute::advanced>("", "buffer", "buffer in GB uint used for permuting the data", buffer, &buffer);
    opts.add<Value<uint>, Attribute::advanced>("", "imaxiter", "maximum number of IRAM interations", imaxiter, &imaxiter);
    opts.add<Value<double>, Attribute::advanced>("", "itol", "tolerance for IRAM algorithm", itol, &itol);
    // opts.add<Switch, Attribute::advanced>("", "mev", "use mev measurement instead of default minSSE", &mev);
    opts.add<Value<uint>, Attribute::advanced>("", "ncv", "number of Lanzcos basis vectors for IRAM", ncv, &ncv);
    opts.add<Value<uint>, Attribute::advanced>("", "oversamples", "number of oversampling columns for RSVD", oversamples, &oversamples);
    opts.add<Value<uint>, Attribute::advanced>("", "rand", "the random matrix type. 0: uniform, 1: guassian", rand, &rand);
    opts.add<Value<double>, Attribute::advanced>("", "tol-rsvd", "tolerance for RSVD algorithm", tol, &tol);
    opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance for EMU/PCAngsd algorithm", tolem, &tolem);
    opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance for MAF estimation updated by EM", tolmaf, &tolmaf);
    opts.add<Switch, Attribute::hidden>("", "groff", "print groff formatted help message", &groff);
    // collect command line options acutal in effect
    ss << (std::string) "PCAone (v" + VERSION + ")    https://github.com/Zilong-Li/PCAone\n";
    ss << "Options in effect:\n";
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(ss, " "));
    try
    {
        opts.parse(argc, argv);
        if(groff)
        {
            GroffOptionPrinter groff_printer(&opts);
            std::cout << groff_printer.print(Attribute::advanced);
            exit(EXIT_SUCCESS);
        }
        if(plinkfile->is_set())
            file_t = FileType::PLINK;
        else if(binfile->is_set())
            file_t = FileType::BINARY;
        else if(bgenfile->is_set())
            file_t = FileType::BGEN;
        else if(beaglefile->is_set())
            file_t = FileType::BEAGLE;
        else if(csvfile->is_set())
            file_t = FileType::CSV;
        else if(help_opt->count() == 1)
        {
            std::cout << opts.help(Attribute::advanced) << "\n";
            exit(EXIT_SUCCESS);
        }
        else if(argc == 1)
        {
            std::cout << opts << "\n";
            exit(EXIT_FAILURE);
        }
        if(svd_opt->value()==0)
            svd_t = SvdType::IRAM;
        else if(svd_opt->value()==1)
            svd_t = SvdType::PCAoneAlg1;
        else if(svd_opt->value()==2)
            svd_t = SvdType::PCAoneAlg2;
        else if(svd_opt->value()==3)
            svd_t = SvdType::FULL, out_of_core = false;
        else
            svd_t = SvdType::PCAoneAlg2;

        ncv = 20 > (2 * k + 1) ? 20 : (2 * k + 1);
        oversamples = 10 > k ? 10 : k;
        // beagle only represents genotype likelihood for pcangsd algorithm now
        if(file_t == FileType::BEAGLE) pcangsd = true;
        if((file_t == FileType::PLINK || file_t == FileType::BGEN) && !haploid) diploid = true;
        if(emu || pcangsd)
            runem = true;
        else
            maxiter = 0;
        if(memory > 0)
        {
            out_of_core = true;
            if(pcangsd)
                throw std::invalid_argument(
                    "not support -m option for PCAngsd algorithm yet, but the feature is on the way!");
        }
        if(bands < 4 || bands % 2 != 0)
            throw std::invalid_argument("the --batches must be a power of 2 and the minimun is 4. the recommended is 64\n");
    }
    catch(const popl::invalid_option & e)
    {
        std::cerr << "Invalid Option Exception: " << e.what() << "\n";
        std::cerr << "error:  ";
        if(e.error() == invalid_option::Error::missing_argument)
            std::cerr << "missing_argument\n";
        else if(e.error() == invalid_option::Error::invalid_argument)
            std::cerr << "invalid_argument\n";
        else if(e.error() == invalid_option::Error::too_many_arguments)
            std::cerr << "too_many_arguments\n";
        else if(e.error() == invalid_option::Error::missing_option)
            std::cerr << "missing_option\n";

        if(e.error() == invalid_option::Error::missing_option)
        {
            std::string option_name(e.option()->name(OptionName::short_name, true));
            if(option_name.empty()) option_name = e.option()->name(OptionName::long_name, true);
            std::cerr << "option: " << option_name << "\n";
        }
        else
        {
            std::cerr << "option: " << e.option()->name(e.what_name()) << "\n";
            std::cerr << "value:  " << e.value() << "\n";
        }
        exit(EXIT_FAILURE);
    }
    catch(const std::exception & e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}

Param::~Param() {}
