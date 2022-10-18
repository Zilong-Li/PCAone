#include "Cmd.hpp"

Param::Param(int argc, char** argv)
{
    std::string copyr{"PCA All In One (v" + (std::string)VERSION +
                 ")        https://github.com/Zilong-Li/PCAone\n(C) 2021-2022 Zilong Li        GNU General Public License v3\n\n\x1B[32mUsage: PCAone [OPTION]. Use --help to see all hidden options\033[0m\n\n"};
    OptionParser opts(copyr + "Main options");
    auto help_opt = opts.add<Switch>("h", "help", "print list of all options including hidden advanced options\n");
    opts.add<Switch>("a", "arnoldi", "use Implicitly Restarted Arnoldi Method instead", &arnoldi);
    opts.add<Value<std::string>>("b", "bfile", "prefix of PLINK .bed/.bim/.fam files", "", &bed_prefix);
    opts.add<Value<std::string>>("c", "csv", "path of CSV file compressed by zstd", "", &csvfile);
    opts.add<Switch>("C", "cpmed", "normalize values by count per median (CPMED) for scRNAs", &cpmed);
    opts.add<Switch>("e", "emu", "use EMU algorithm for genotype data with lots of missingness", &emu);
    opts.add<Value<double>>("f", "maf", "remove variants with minor allele frequency below maf", maf, &maf);
    opts.add<Value<std::string>>("g", "bgen", "path of BGEN file", "", &bgen);
    opts.add<Value<uint>>("k", "eigs", "top k components to be calculated", k, &k);
    opts.add<Value<std::string>>("l", "beagle", "path of beagle file", "", &beagle);
    opts.add<Value<double>>("m", "memory", "specify the RAM usage in GB unit", memory, &memory);
    opts.add<Value<uint>>("n", "threads", "number of threads for multithreading", threads, &threads);
    opts.add<Switch>("y", "halko", "use Yu RSVD + Halko power iteration algorithm instead", &halko);
    opts.add<Switch>("S", "no-shuffle", "do not shuffle the data if it is already permuted", &noshuffle);
    opts.add<Value<uint64>, Attribute::advanced>("M", "", "number of features (eg. SNPs) if already known", 0, &nsnps);
    opts.add<Value<uint64>, Attribute::advanced>("N", "", "number of samples if already known", 0, &nsamples);
    opts.add<Value<std::string>>("o", "out", "prefix of output files", outfile, &outfile);
    opts.add<Switch>("p", "pcangsd", "use PCAngsd algorithm for genotype likelihood input", &pcangsd);
    opts.add<Value<uint>>("P", "maxp", "maximum number of power iterations for RSVD algorithm", maxp, &maxp);
    opts.add<Value<std::string>>("T", "tmp", "prefix of temporary permuted data", "", &tmpfile);
    opts.add<Switch>("U", "printu", "output eigen vector of each epoch (for tests)", &printu);
    opts.add<Switch>("v", "verbose", "verbose message output", &verbose);
    opts.add<Switch>("V", "printv", "output the right eigen vectors with suffix .loadings", &printv);
    opts.add<Value<uint>, Attribute::advanced>("", "buffer", "buffer in GB uint used for permuting the data", buffer, &buffer);
    opts.add<Value<uint>, Attribute::advanced>("", "imaxiter", "maximum number of IRAM interations", imaxiter, &imaxiter);
    opts.add<Value<double>, Attribute::advanced>("", "itol", "tolerance for IRAM algorithm", itol, &itol);
    // opts.add<Switch, Attribute::advanced>("", "mev", "use mev measurement instead of default minSSE", &mev);
    opts.add<Value<uint>, Attribute::advanced>("", "ncv", "number of Lanzcos basis vectors for IRAM", ncv, &ncv);
    opts.add<Value<uint>, Attribute::advanced>("", "oversamples", "number of oversampling columns for RSVD", oversamples, &oversamples);
    opts.add<Value<uint>, Attribute::advanced>("", "rand", "the random matrix type. 0: uniform, 1: guassian", rand, &rand);
    opts.add<Value<double>, Attribute::advanced>("", "tol", "tolerance for RSVD algorithm", tol, &tol);
    opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance for EMU/PCAngsd algorithm", tolem, &tolem);
    opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance for MAF estimation updated by EM", tolmaf,
                                                 &tolmaf);
    opts.add<Value<int>, Attribute::advanced>("", "windows", "number of windows to use for PCAone (algorithm2)", bands, &bands);
    opts.add<Switch, Attribute::hidden>("", "groff", "print groff formatted help message", &groff);
    // collect command line options acutal in effect
    ss << (std::string) "PCAone (v" + VERSION + ")    https://github.com/Zilong-Li/PCAone\n";
    ss << "Options in effect:\n";
    std::copy(argv, argv + argc, std::ostream_iterator<char*>(ss, " "));
    try
    {
        opts.parse(argc, argv);
        if (groff)
        {
            GroffOptionPrinter groff_printer(&opts);
            std::cout << groff_printer.print(Attribute::advanced);
            exit(EXIT_SUCCESS);
        }
        if (bed_prefix != "")
        {
            intype = FileType::PLINK;
        }
        else if (bgen != "")
        {
            intype = FileType::BGEN;
        }
        else if (beagle != "")
        {
            intype = FileType::BEAGLE;
        }
        else if (csvfile != "")
        {
            intype = FileType::CSV;
        }
        else if (help_opt->count() == 1)
        {
            std::cout << opts.help(Attribute::advanced) << "\n";
            exit(EXIT_SUCCESS);
        }
        else if (argc == 1)
        {
            std::cout << opts << "\n";
            exit(EXIT_FAILURE);
        }
        ncv = 20 > (2 * k + 1) ? 20 : (2 * k + 1);
        oversamples = 10 > k ? 10 : k;
        // beagle only represents genotype likelihood for pcangsd algorithm now
        if (intype == FileType::BEAGLE)
            pcangsd = true;
        if (emu || pcangsd)
        {
            runem = true;
        }
        else
        {
            maxiter = 0;
        }
        if (halko || arnoldi)
            fast = false;
        if (memory > 0)
        {
            batch = false;
            if (pcangsd)
                throw std::invalid_argument("not support -m option for PCAngsd algorithm yet, but the feature is on the way!");
        }
    }
    catch (const popl::invalid_option& e)
    {
        std::cerr << "Invalid Option Exception: " << e.what() << "\n";
        std::cerr << "error:  ";
        if (e.error() == invalid_option::Error::missing_argument)
            std::cerr << "missing_argument\n";
        else if (e.error() == invalid_option::Error::invalid_argument)
            std::cerr << "invalid_argument\n";
        else if (e.error() == invalid_option::Error::too_many_arguments)
            std::cerr << "too_many_arguments\n";
        else if (e.error() == invalid_option::Error::missing_option)
            std::cerr << "missing_option\n";

        if (e.error() == invalid_option::Error::missing_option)
        {
            std::string option_name(e.option()->name(OptionName::short_name, true));
            if (option_name.empty())
                option_name = e.option()->name(OptionName::long_name, true);
            std::cerr << "option: " << option_name << "\n";
        }
        else
        {
            std::cerr << "option: " << e.option()->name(e.what_name()) << "\n";
            std::cerr << "value:  " << e.value() << "\n";
        }
        exit(EXIT_FAILURE);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

}

Param::~Param()
{
}
