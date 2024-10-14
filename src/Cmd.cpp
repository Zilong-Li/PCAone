/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Cmd.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Cmd.hpp"

#include <cstdlib>
#include <iostream>

#include "popl/popl.hpp"

using namespace popl;

Param::Param(int argc, char **argv) {
  // clang-format off
  bool haploid = false;
  std::string copyr{"PCA All In One (v" + (std::string)VERSION + ")        https://github.com/Zilong-Li/PCAone\n" +
                    "(C) 2021-2024 Zilong Li        GNU General Public License v3\n" +
  "\x1B[32m\n" +
                    "Usage: use plink files as input and apply default window-based RSVD method\n" +
                    "       PCAone --bfile plink -n 20 \n\n" +
                    "       use csv file as input and apply the Implicitly Restarted Arnoldi Method\n" +
                    "       PCAone --csv csv.zst --svd 0 \n" +
  "\033[0m\n"};
  OptionParser opts(copyr + "Main options");
  auto help_opt = opts.add<Switch>("h", "help", "print all options including hidden advanced options");
  auto svd_opt = opts.add<Value<uint>>("d", "svd", "SVD method to be applied. default 2 is recommended for big data.\n"
                                       "0: the Implicitly Restarted Arnoldi Method (IRAM)\n"
                                       "1: the Yu's single-pass Randomized SVD with power iterations\n"
                                       "2: the accurate window-based Randomized SVD method (PCAone)\n"
                                       "3: the full Singular Value Decomposition.", 2);
  auto plinkfile = opts.add<Value<std::string>>("b", "bfile", "prefix to PLINK .bed/.bim/.fam files", "", &filein);
  auto binfile = opts.add<Value<std::string>>("B", "binary", "path of binary file", "", &filein);
  auto csvfile = opts.add<Value<std::string>>("c", "csv", "path of comma seperated CSV file compressed by zstd", "", &filein);
  auto bgenfile = opts.add<Value<std::string>>("g", "bgen", "path of BGEN file compressed by gzip/zstd", "", &filein);
  auto beaglefile = opts.add<Value<std::string>>("G", "beagle", "path of BEAGLE file compressed by gzip", "", &filein);
  opts.add<Value<std::string>>("", "read-U", "path of file with left singular vectors (.eigvecs)", "", &fileU);
  opts.add<Value<std::string>>("", "read-V", "path of file with right singular vectors (.loadings)", "", &fileV);
  opts.add<Value<std::string>>("", "read-S", "path of file with eigen values (.eigvals)", "", &fileS);
  opts.add<Value<uint>>("k", "pc", "top k principal components (PCs) to be calculated", k, &k);
  opts.add<Value<double>>("m", "memory", "RAM usage in GB unit for out-of-core mode. default is in-core mode", memory, &memory);
  opts.add<Value<uint>>("n", "threads", "the number of threads to be used", threads, &threads);
  opts.add<Value<std::string>>("o", "out", "prefix to output files. default [pcaone]", fileout, &fileout);
  opts.add<Value<uint>>("p", "maxp", "maximum number of power iterations for RSVD algorithm", maxp, &maxp);
  opts.add<Switch>("S", "no-shuffle", "do not shuffle columns of data for --svd 2 (if not locally correlated)", &noshuffle);
  opts.add<Switch>("v", "verbose", "verbose message output", &verbose);
  opts.add<Switch>("V", "printv", "output the right eigenvectors with suffix .loadings", &printv);
  opts.add<Value<uint>, Attribute::advanced>("w", "batches", "the number of mini-batches used by --svd 2", bands, &bands);
  opts.add<Value<uint>>("C", "scale", "do scaling for input file.\n"
                        "0: do just centering\n"
                        "1: do log transformation eg. log(x+0.01) for RNA-seq data\n"
                        "2: do count per median log transformation (CPMED) for scRNAs",
                        scale,  &scale);
  opts.add<Switch>("", "emu", "use EMU algorithm for genotype input with missingness", &emu);
  opts.add<Switch>("", "pcangsd", "use PCAngsd algorithm for genotype likelihood input", &pcangsd);
  opts.add<Value<double>>("", "maf", "exclude variants with MAF lower than this value", maf, &maf);
  opts.add<Value<int>>("", "project", "project the new samples onto the existing PCs.\n"
                                      "0: disabled\n"
                                      "1: by multiplying the loadings with mean imputation for missing genotypes\n"
                                      "2: by solving the least squares system Vx=g. skip sites with missingness\n"
                                      "3: by Augmentation, Decomposition and Procrusters transformation\n",
                       project, &project);
  opts.add<Value<int>>("", "inbreed", "compute the inbreeding coefficient accounting for population structure.\n"
                                      "0: disabled\n"
                                      "1: compute per-site inbreeding coefficient and HWE test\n",
                       inbreed, &inbreed);
  opts.add<Value<std::string>>("", "match-bim", "the .mbim file to be matched, where the 7th column is allele frequency", "", &filebim);
  opts.add<Switch>("", "ld", "output a binary matrix for downstream LD related analysis", &ld);
  opts.add<Value<double>>("", "ld-r2", "r2 cutoff for LD-based pruning. (usually 0.2)", ld_r2, &ld_r2);
  opts.add<Value<uint>>("", "ld-bp", "physical distance threshold in bases for LD. (usually 1000000)", ld_bp, &ld_bp);
  opts.add<Value<int>>("", "ld-stats", "statistics to calculate LD r2 for pairwise SNPs.\n"
                                       "0: the ancestry adjusted, i.e. correlation between residuals\n"
                                       "1: the standard, i.e. correlation between two alleles\n",
                       ld_stats, &ld_stats);
  opts.add<Switch>("", "print-r2", "print LD r2 to *.ld.gz file for pairwise SNPs within a window", &print_r2);
  auto clumpfile = opts.add<Value<std::string>>("", "clump", "assoc-like file with target variants and pvalues for clumping", "", &clump);
  auto assocnames = opts.add<Value<std::string>>("", "clump-names", "column names in assoc-like file for locating chr, pos and pvalue", "CHR,BP,P", &assoc_colnames);
  opts.add<Value<double>>("", "clump-p1", "significance threshold for index SNPs", clump_p1, &clump_p1);
  opts.add<Value<double>>("", "clump-p2", "secondary significance threshold for clumped SNPs", clump_p2, &clump_p2);
  opts.add<Value<double>>("", "clump-r2", "r2 cutoff for LD-based clumping", clump_r2, &clump_r2);
  opts.add<Value<uint>>("", "clump-bp", "physical distance threshold in bases for clumping", clump_bp, &clump_bp);
  opts.add<Switch, Attribute::advanced>("", "printu", "output eigen vector of each epoch (for tests)", &printu);
  opts.add<Value<uint>, Attribute::advanced>("", "M", "the number of features (eg. SNPs) if already known", 0, &nsnps);
  opts.add<Value<uint>, Attribute::advanced>("", "N", "the number of samples if already known", 0, &nsamples);
  // opts.add<Switch, Attribute::advanced>("", "debug", "turn on debugging mode", &debug);
  opts.add<Switch, Attribute::advanced>("", "haploid", "the plink format represents haploid data", &haploid);
  opts.add<Value<uint>, Attribute::advanced>("", "buffer", "memory buffer in GB unit for permuting the data", buffer, &buffer);
  opts.add<Value<uint>, Attribute::advanced>("", "imaxiter", "maximum number of IRAM iterations", imaxiter, &imaxiter);
  opts.add<Value<double>, Attribute::advanced>("", "itol", "stopping tolerance for IRAM algorithm", itol, &itol);
  opts.add<Value<uint>, Attribute::advanced>("", "ncv", "the number of Lanzcos basis vectors for IRAM", ncv, &ncv);
  opts.add<Value<uint>, Attribute::advanced>("", "oversamples", "the number of oversampling columns for RSVD", oversamples, &oversamples);
  opts.add<Value<uint>, Attribute::advanced>("", "rand", "the random matrix type. 0: uniform, 1: guassian", rand, &rand);
  opts.add<Value<double>, Attribute::advanced>("", "tol-rsvd", "tolerance for RSVD algorithm", tol, &tol);
  opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance for EMU/PCAngsd algorithm", tolem, &tolem);
  opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance for MAF estimation by EM", tolmaf, &tolmaf);
  opts.add<Switch, Attribute::hidden>("", "groff", "print groff formatted help message", &groff);
  // collect command line options acutal in effect
  ss << (std::string) "PCAone (v" + VERSION + ")    https://github.com/Zilong-Li/PCAone\n";
  ss << "Options in effect:\n";
  std::copy(argv, argv + argc, std::ostream_iterator<char *>(ss, " "));
  // clang-format on
  try {
    opts.parse(argc, argv);
    if (groff) {
      GroffOptionPrinter groff_printer(&opts);
      std::cout << groff_printer.print(Attribute::advanced);
      exit(EXIT_SUCCESS);
    }
    if (!opts.unknown_options().empty()) {
      for (const auto &uo : opts.unknown_options()) std::cerr << "unknown option: " << uo << "\n";
      exit(EXIT_FAILURE);
    }
    if (plinkfile->is_set())
      file_t = FileType::PLINK;
    else if (binfile->is_set())
      file_t = FileType::BINARY;
    else if (bgenfile->is_set())
      file_t = FileType::BGEN;
    else if (beaglefile->is_set())
      file_t = FileType::BEAGLE;
    else if (csvfile->is_set())
      file_t = FileType::CSV;
    else if (help_opt->count() == 1) {
      std::cout << opts.help(Attribute::advanced) << "\n";
      exit(EXIT_SUCCESS);
    } else if (argc == 1) {
      std::cout << opts << "\n";
      exit(EXIT_SUCCESS);
    }
    if (svd_opt->value() == 0)
      svd_t = SvdType::IRAM;
    else if (svd_opt->value() == 1)
      svd_t = SvdType::PCAoneAlg1;
    else if (svd_opt->value() == 2)
      svd_t = SvdType::PCAoneAlg2;
    else if (svd_opt->value() == 3)
      svd_t = SvdType::FULL;
    else
      svd_t = SvdType::PCAoneAlg2;

    ncv = 20 > (2 * k + 1) ? 20 : (2 * k + 1);
    oversamples = 10 > k ? 10 : k;
    // beagle.gz only represents genotype likelihood for pcangsd algorithm now
    if (file_t == FileType::BEAGLE && svd_t == SvdType::PCAoneAlg2)
      throw std::invalid_argument(
          "not supporting PCAone --svd 2 for PCAngsd algorithm yet! "
          "please use --svd 1 or 0 option");
    if (file_t == FileType::BEAGLE) pcangsd = true;
    if (haploid && (file_t == FileType::PLINK || file_t == FileType::BGEN)) ploidy = 1;
    if (emu || pcangsd) {
      impute = true;
      if (svd_t == SvdType::PCAoneAlg2)
        throw std::invalid_argument(
            "not supporting PCAone --svd 2 for PCAngsd/EMU algorithm yet! "
            "please specify --svd 1 or 0");
    } else {
      maxiter = 0;
    }
    if (memory > 0) {
      out_of_core = true;
      if (pcangsd) throw std::invalid_argument("not supporting -m option for PCAngsd algorithm yet!");
    }
    if (bands < 4 || bands % 2 != 0)
      throw std::invalid_argument(
          "the -w/--batches must be a power of 2 and the minimun is 4. "
          "the recommended is 64.");
    if (maf > 0.5) {
      std::cerr << "warning: you specify '--maf' a value greater than 0.5. "
                   "will do 1-maf for you!\n";
      maf = 1 - maf;
    }
    keepsnp = maf > 0 ? true : false;
    if (maf && out_of_core)
      throw std::invalid_argument(
          "does not support --maf filters for out-of-core mode yet! "
          "please remove --maf option");
    if (svd_t == SvdType::PCAoneAlg2 && !noshuffle) perm = true;
    if (svd_t == SvdType::FULL) out_of_core = false;
    if (print_r2 || ld_bp > 0 || ld_r2 > 0 || !clump.empty()) {
      pca = false;
      memory /= 2.0;  // adjust memory estimator
    }

    // handle projection
    if (project > 0) {
      if (fileV.empty() || fileS.empty())
        throw std::invalid_argument("please use --read-S and --read-V together with --project");
      if (project > 2) throw std::invalid_argument("more projection methods are coming. stay tuned!");
      impute = true;
      memory = 0;
      out_of_core = false;
    }

    // handle inbreeding
    if (inbreed > 0) {
      if (fileU.empty() || fileV.empty() || fileS.empty())
        throw std::invalid_argument("please use --read-U, --read-S and --read-V together with --inbreed");
    }
  } catch (const popl::invalid_option &e) {
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

    if (e.error() == invalid_option::Error::missing_option) {
      std::string option_name(e.option()->name(OptionName::short_name, true));
      if (option_name.empty()) option_name = e.option()->name(OptionName::long_name, true);
      std::cerr << "option: " << option_name << "\n";
    } else {
      std::cerr << "option: " << e.option()->name(e.what_name()) << "\n";
      std::cerr << "value:  " << e.value() << "\n";
    }
    exit(EXIT_FAILURE);
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
}

Param::~Param() {}
