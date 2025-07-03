/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Cmd.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Cmd.hpp"

#include "popl/popl.hpp"

using namespace popl;

Param::Param(int argc, char **argv) {
  // clang-format off
  bool haploid = false;
  std::string copyr{"PCA All In One (v" + (std::string)VERSION + ")        https://github.com/Zilong-Li/PCAone\n" +
                    "(C) 2021-2024 Zilong Li        GNU General Public License v3\n" +
  "\n" +
                    "Usage: 1) use PLINK files as input and apply default window-based RSVD method\n" +
                    "       $ PCAone -b plink \n\n" +
                    "       2) use CSV file as input and apply the Implicitly Restarted Arnoldi Method\n" +
                    "       $ PCAone -c csv.zst -d 0 \n\n" +
                    "       3) compute ancestry adjusted LD matrix and R2\n" +
                    "       $ PCAone -b plink -k 2 -D -o adj \n" +
                    "       $ PCAone -B adj.residuals -f adj.mbim -R --ld-bp 1000" +
  "\n"};
  OptionParser opts(copyr);
  opts.add<Value<std::string>, Attribute::headline>("","PCAone","General options:");
  auto help_opt = opts.add<Switch>("h", "help", "print all options including hidden advanced options");
  opts.add<Value<double>>("m", "memory", "RAM usage in GB unit for out-of-core mode. default is in-core mode", memory, &memory);
  opts.add<Value<uint>>("n", "threads", "the number of threads to be used", threads, &threads);
  opts.add<Value<uint>>("v", "verbose", "verbosity level for logs. any level x includes messages for all levels (1...x). Options are\n"
                                        "0: silent, no messages on screen;\n"
                                        "1: concise messages to screen;\n"
                                        "2: more verbose information;\n"
                                        "3: enable debug information."
                        , verbose, &verbose);
  opts.add<Value<std::string>, Attribute::headline>("","PCA","PCA algorithms:");
  auto svd_opt = opts.add<Value<uint>>("d", "svd", "SVD method to be applied. default 2 is recommended for big data. Options are\n"
                                                   "0: the Implicitly Restarted Arnoldi Method (IRAM);\n"
                                                   "1: the Yu's single-pass Randomized SVD with power iterations;\n"
                                                   "2: the accurate window-based Randomized SVD method (PCAone);\n"
                                                   "3: the full Singular Value Decomposition.", 2);
  opts.add<Value<uint>>("k", "pc", "top k principal components (PCs) to be calculated", k, &k);
  opts.add<Value<uint>>("C", "scale", "do scaling for input file. Options are\n"
                                      "0: do nothing and proceed to SVD;\n"
                                      "1: do only standardization;\n"
                                      "2: do count per median log transformation (CPMED);\n"
                                      "3: do log1p transformation;\n"
                                      "4: do relative counts.", scale,  &scale);
  opts.add<Value<uint>>("p", "maxp", "maximum number of power iterations for RSVD algorithm.", maxp, &maxp);
  opts.add<Switch>("S", "no-shuffle", "do not shuffle columns of data for --svd 2 (if not locally correlated).", &noshuffle);
  opts.add<Value<uint>, Attribute::advanced>("w", "batches", "the number of mini-batches used by --svd 2.", bands, &bands);
  opts.add<Switch>("", "emu", "use EMU algorithm for genotype input with missingness.", &emu);
  opts.add<Switch>("", "pcangsd", "use PCAngsd algorithm for genotype likelihood input.", &pcangsd);
  opts.add<Value<uint>, Attribute::advanced>("", "M", "the number of features (eg. SNPs) if already known.", 0, &nsnps);
  opts.add<Value<uint>, Attribute::advanced>("", "N", "the number of samples if already known.", 0, &nsamples);
  opts.add<Value<double>, Attribute::advanced>("", "scale-factor", "feature counts for each sample are normalized and multiplied by this value", 1.0, &scaleFactor);
  opts.add<Value<uint>, Attribute::advanced>("", "buffer", "memory buffer in GB unit for permuting the data.", buffer, &buffer);
  opts.add<Value<uint>, Attribute::advanced>("", "imaxiter", "maximum number of IRAM iterations.", imaxiter, &imaxiter);
  opts.add<Value<double>, Attribute::advanced>("", "itol", "stopping tolerance for IRAM algorithm.", itol, &itol);
  opts.add<Value<uint>, Attribute::advanced>("", "ncv", "the number of Lanzcos basis vectors for IRAM.", ncv, &ncv);
  opts.add<Value<uint>, Attribute::advanced>("", "oversamples", "the number of oversampling columns for RSVD.", oversamples, &oversamples);
  opts.add<Value<uint>, Attribute::advanced>("", "rand", "the random matrix type. 0: uniform; 1: guassian.", rand, &rand);
  opts.add<Value<uint>, Attribute::advanced>("", "maxiter", "maximum number of EM iterations.", maxiter, &maxiter);
  opts.add<Value<double>, Attribute::advanced>("", "tol-rsvd", "tolerance for RSVD algorithm.", tol, &tol);
  opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance for EMU/PCAngsd algorithm.", tolem, &tolem);
  opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance for MAF estimation by EM.", tolmaf, &tolmaf);
  
  opts.add<Value<std::string>, Attribute::headline>("","INPUT","Input options:");
  auto plinkfile = opts.add<Value<std::string>>("b", "bfile", "prefix of PLINK .bed/.bim/.fam files.", "", &filein);
  opts.add<Switch, Attribute::advanced>("", "haploid", "the plink format represents haploid data.", &haploid);
  auto binfile = opts.add<Value<std::string>>("B", "binary", "path of binary file.", "", &filein);
  auto csvfile = opts.add<Value<std::string>>("c", "csv", "path of comma seperated CSV file compressed by zstd.", "", &filein);
  auto bgenfile = opts.add<Value<std::string>>("g", "bgen", "path of BGEN file compressed by gzip/zstd.", "", &filein);
  auto beaglefile = opts.add<Value<std::string>>("G", "beagle", "path of BEAGLE file compressed by gzip.", "", &filein);
  opts.add<Value<std::string>>("f", "match-bim", "the .mbim file to be matched, where the 7th column is allele frequency.", "", &filebim);
  auto usvprefix = opts.add<Value<std::string>>("", "USV", "prefix of PCAone .eigvecs/.eigvals/.loadings/.mbim.");
  opts.add<Value<std::string>, Attribute::hidden>("", "read-U", "path of file with left singular vectors (.eigvecs).", "", &fileU);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-V", "path of file with right singular vectors (.loadings).", "", &fileV);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-S", "path of file with sigular values (.sigvals).", "", &fileS);
  
  opts.add<Value<std::string>, Attribute::headline>("","OUTPUT","Output options:");
  opts.add<Value<std::string>>("o", "out", "prefix of output files. default [pcaone].", fileout, &fileout);
  opts.add<Switch>("V", "printv", "output the right eigenvectors with suffix .loadings.", &printv);
  opts.add<Switch>("D", "ld", "output a binary matrix for downstream LD related analysis.", &ld);
  opts.add<Switch>("R", "print-r2", "print LD R2 to *.ld.gz file for pairwise SNPs within a window controlled by --ld-bp.", &print_r2);
  
  opts.add<Value<std::string>, Attribute::headline>("","MISC","Misc options:");
  opts.add<Value<double>>("", "maf", "exclude variants with MAF lower than this value", maf, &maf);
  opts.add<Value<int>>("", "project", "project the new samples onto the existing PCs. Options are\n"
                                      "0: disabled;\n"
                                      "1: by multiplying the loadings with mean imputation for missing genotypes;\n"
                                      "2: by solving the least squares system Vx=g. skip sites with missingness;\n"
                                      "3: by Augmentation, Decomposition and Procrusters transformation.\n", project, &project);
  opts.add<Value<int>>("", "inbreed", "compute the inbreeding coefficient accounting for population structure. Options are\n"
                                      "0: disabled;\n"
                                      "1: compute per-site inbreeding coefficient and HWE test.\n", inbreed, &inbreed);
  opts.add<Value<double>>("", "ld-r2", "R2 cutoff for LD-based pruning (usually 0.2).", ld_r2, &ld_r2);
  opts.add<Value<uint>>("", "ld-bp", "physical distance threshold in bases for LD window (usually 1000000).", ld_bp, &ld_bp);
  opts.add<Value<int>>("", "ld-stats", "statistics to compute LD R2 for pairwise SNPs. Options are\n"
                                       "0: the ancestry adjusted, i.e. correlation between residuals;\n"
                                       "1: the standard, i.e. correlation between two alleles.\n", ld_stats, &ld_stats);
  auto clumpfile = opts.add<Value<std::string>>("", "clump", "assoc-like file with target variants and pvalues for clumping.", "", &clump);
  auto assocnames = opts.add<Value<std::string>>("", "clump-names", "column names in assoc-like file for locating chr, pos and pvalue.", "CHR,BP,P", &assoc_colnames);
  opts.add<Value<double>>("", "clump-p1", "significance threshold for index SNPs.", clump_p1, &clump_p1);
  opts.add<Value<double>>("", "clump-p2", "secondary significance threshold for clumped SNPs.", clump_p2, &clump_p2);
  opts.add<Value<double>>("", "clump-r2", "r2 cutoff for LD-based clumping.", clump_r2, &clump_r2);
  opts.add<Value<uint>>("", "clump-bp", "physical distance threshold in bases for clumping.", clump_bp, &clump_bp);
  opts.add<Switch, Attribute::hidden>("", "groff", "PCAone 1 \"24 December 2024\" \"PCAone-v"+ std::string(VERSION)+"\"  \"Bioinformatics tools\"", &groff);
  
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
    // handle PI, i.e U,S,V
    if (usvprefix->is_set()) {
      if (fileU.empty()) fileU = usvprefix->value() + ".eigvecs";
      if (fileS.empty()) fileS = usvprefix->value() + ".sigvals";
      if (fileV.empty()) fileV = usvprefix->value() + ".loadings";
      if (filebim.empty()) filebim = usvprefix->value() + ".mbim";
    }

    // handle LD
    if (print_r2 || ld_bp > 0 || ld_r2 > 0 || !clump.empty()) {
      pca = false;
      memory /= 2.0;  // adjust memory estimator
    }

    // handle projection
    if (project > 0) {
      if (fileV.empty() || fileS.empty())
        throw std::invalid_argument(
            "please use --read-S and --read-V together with --project, or simply --USV");
      if (project > 2) throw std::invalid_argument("more projection methods are coming. stay tuned!");
      estaf = false, impute = true, out_of_core = false, pca = false;
      memory = 0;
    }

    // handle inbreeding
    if (inbreed > 0) {
      estaf = false, center = false, pca = false;
      if (fileU.empty() || fileV.empty() || fileS.empty())
        throw std::invalid_argument("please use --USV together with --inbreed");
    }

    // handle memory and misc options
    ncv = 20 > (2 * k + 1) ? 20 : (2 * k + 1);
    oversamples = oversamples > k ? oversamples : k;
    if (haploid && (file_t == FileType::PLINK || file_t == FileType::BGEN)) ploidy = 1;
    if (memory > 0 && svd_t != SvdType::FULL) out_of_core = true;
    if (maf > 0.5) {
      std::cerr << "warning: you specify '--maf' a value greater than 0.5.\n";
      maf = 1 - maf;
    }
    keepsnp = maf > 0 ? true : false;
    if (maf && out_of_core)
      throw std::invalid_argument("does not support --maf filters for out-of-core mode yet! ");

    // handle EM-PCA
    if (pca && file_t == FileType::BEAGLE) pcangsd = true;
    if (emu || pcangsd) {
      impute = true;
      if (svd_t == SvdType::PCAoneAlg2)
        throw std::invalid_argument("fancy EM-PCA with --svd 2 is on the way!");
    } else if (pca) {
      maxiter = 0;
    }
    if (out_of_core && pca && pcangsd)
      throw std::invalid_argument("not supporting -m option for PCAngsd with BEAGLE file yet!");
    if (bands < 4 || bands % 2 != 0)
      throw std::invalid_argument("the -w/--batches must be a power of 2 and the minimun is 4.");

    if (svd_t == SvdType::PCAoneAlg2 && !noshuffle) perm = true;

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
