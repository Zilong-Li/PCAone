/*******************************************************************************
 * @file        https://github.com/Zilong-Li/PCAone/src/Cmd.cpp
 * @author      Zilong Li
 * Copyright (C) 2022-2024. Use of this code is governed by the LICENSE file.
 ******************************************************************************/

#include "Cmd.hpp"

#include <cstring>
#include <iterator>

#include "popl/popl.hpp"

using namespace popl;

namespace {

// popl reads "-1" into an unsigned option as 4294967295 (istream wraps it), so
// a negative value would pass every range check below. Refuse the sign instead.
class Unsigned : public Value<uint> {
 public:
  using Value<uint>::Value;

 protected:
  void parse(OptionName what_name, const char* value) override {
    if (value != nullptr && std::strchr(value, '-') != nullptr)
      throw invalid_option(
          this, invalid_option::Error::invalid_argument, what_name, value,
          "invalid argument for " + name(what_name, true) + ": '" + value + "' (must be a non-negative integer)");
    Value<uint>::parse(what_name, value);
  }
};

void require(bool ok, const std::string& msg) {
  if (!ok) throw std::invalid_argument(msg);
}

// true if s is exactly three non-empty comma-separated names, e.g. CHR,BP,P
bool three_names(const std::string& s) {
  std::istringstream in(s);
  std::string f;
  int n = 0;
  while (std::getline(in, f, ',')) {
    if (f.empty()) return false;
    ++n;
  }
  return n == 3 && s.back() != ',';
}

}  // namespace

Param::Param(int argc, char** argv) {
  // clang-format off
  bool haploid = false;
  std::string copyr{"PCA All In One (v" + (std::string)VERSION + ")        https://github.com/Zilong-Li/PCAone\n" +
                    "(C) 2021-2024 Zilong Li        GNU General Public License v3\n" +
  "\n" +
                    "Usage: 1) use PLINK files as input and apply default window-based RSVD method\n" +
                    "       $ PCAone -b plink \n\n" +
                    "       2) use CSV file as input and apply the Implicitly Restarted Arnoldi Method\n" +
                    "       $ PCAone -c csv.zst -d 0 \n\n" +
                    "       3) compute the ancestry adjusted LD R2, removing the PCs of a previous run\n" +
                    "       $ PCAone -b plink -k 2 -o pcs \n" +
                    "       $ PCAone -b plink -P pcs -R --ld-bp 1000 -o adj" +
  "\n"};
  OptionParser opts(copyr);
  opts.add<Value<std::string>, Attribute::headline>("","PCAone","General options:");
  auto help_opt = opts.add<Switch>("h", "help", "print all options including hidden advanced options");
  opts.add<Value<double>>("m", "memory", "RAM usage in GB unit for out-of-core mode. default is in-core mode.\n"
                                        "with --svd 3, it sets the blocks the GRM is streamed in", memory, &memory);
  opts.add<Unsigned>("n", "threads", "the number of threads to be used", threads, &threads);
  opts.add<Unsigned>("v", "verbose", "verbosity level for logs. Options are\n"
                                     "0: silent, no messages on screen;\n"
                                     "1: concise messages to screen;\n"
                                     "2: more verbose information;\n"
                                     "3: enable debug information."
                     , verbose, &verbose);
  opts.add<Value<std::string>, Attribute::headline>("","PCA","PCA algorithms:");
  auto svd_opt = opts.add<Unsigned>("d", "svd", "SVD method to be applied. default 2 is recommended for big data. Options are\n"
                                                "0: the Implicitly Restarted Arnoldi Method (IRAM);\n"
                                                "1: the Yu's single-pass Randomized SVD with power iterations;\n"
                                                "2: the accurate window-based Randomized SVD method (PCAone);\n"
                                                "3: exact PCA by eigendecomposition of the sample GRM, streamed block by block when N <= M,\n"
                                                "   in N x N memory (no EM-PCA support).", 2);
  opts.add<Unsigned>("k", "pc", "top k principal components (PCs) to be calculated", k, &k);
  opts.add<Value<int>>("C", "scale", "do normalization or scaling for input file. Options are\n"
                                     "-9: standardize genetic data by sqrt(ploidy*f*(1-f));\n"
                                     " 0: do nothing and proceed to SVD;\n"
                                     " 1: do direct standardization, as the scale(x, center=TRUE, scale=TRUE) function in R;\n"
                                     " 2: do first count per median log transformation (CPMED), then standardization;\n"
                                     " 3: do first log1p transformation, then standardization;\n"
                                     " 4: do first relative counts, then standardization.", scale,  &scale);
  opts.add<Unsigned>("", "maxp", "maximum number of power iterations for RSVD algorithm.", maxp, &maxp);
  opts.add<Switch>("S", "no-shuffle", "do not shuffle columns of data for --svd 2 (if not locally correlated).", &noshuffle);
  opts.add<Unsigned, Attribute::advanced>("w", "batches", "the number of mini-batches used by --svd 2.", bands, &bands);
  opts.add<Value<int>>("", "seed", "seeds for reproducing results.\n", seed, &seed);
  opts.add<Switch>("", "emu", "use EMU algorithm for genotype input with missingness. not with --svd 3.", &emu);
  opts.add<Switch>("", "pcangsd", "use PCAngsd algorithm for genotype likelihood input. not with --svd 3.", &pcangsd);
  opts.add<Unsigned, Attribute::advanced>("", "M", "the number of features (eg. SNPs) if already known.", 0, &nsnps);
  opts.add<Unsigned, Attribute::advanced>("", "N", "the number of samples if already known.", 0, &nsamples);
  opts.add<Value<double>, Attribute::advanced>("", "scale-factor", "feature counts for each sample are normalized and multiplied by this value", 1.0, &scaleFactor);
  opts.add<Unsigned, Attribute::advanced>("", "buffer", "memory buffer in GB unit for permuting the data.", buffer, &buffer);
  opts.add<Unsigned, Attribute::advanced>("", "imaxiter", "maximum number of IRAM iterations.", imaxiter, &imaxiter);
  opts.add<Value<double>, Attribute::advanced>("", "itol", "stopping tolerance for IRAM algorithm.", itol, &itol);
  auto ncv_opt = opts.add<Unsigned, Attribute::advanced>("", "ncv", "the number of Lanzcos basis vectors for IRAM.", ncv, &ncv);
  opts.add<Unsigned, Attribute::advanced>("", "oversamples", "the number of oversampling columns for RSVD.", oversamples, &oversamples);
  opts.add<Unsigned, Attribute::advanced>("", "rand", "the random matrix type. 0: uniform; 1: guassian.", rand, &rand);
  opts.add<Unsigned, Attribute::advanced>("", "maxiter", "maximum number of EM iterations.", maxiter, &maxiter);
  opts.add<Value<double>, Attribute::advanced>("", "tol-rsvd", "tolerance for RSVD algorithm.", tol, &tol);
  opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance for EMU/PCAngsd algorithm.", tolem, &tolem);
  opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance for MAF estimation by EM.", tolmaf, &tolmaf);
  
  opts.add<Value<std::string>, Attribute::headline>("","INPUT","Input options:");
  auto plinkfile = opts.add<Value<std::string>>("b", "bfile", "prefix of PLINK .bed/.bim/.fam files.", "", &filein);
  opts.add<Switch, Attribute::advanced>("", "haploid", "the plink format represents haploid data.", &haploid);
  auto pgenfile = opts.add<Value<std::string>>("p", "pgen", "prefix of PLINK2 .pgen/.pvar/.psam files.", "", &filein);
  opts.add<Switch, Attribute::advanced>("", "hardcall", "use hardcall genotype instead of dosages.", &hardcall);
  // removed in v0.8.0; kept hidden only to say what replaces them
  auto binfile = opts.add<Value<std::string>, Attribute::hidden>("B", "binary", "removed. LD now reads the genotypes directly.");
  auto csvfile = opts.add<Value<std::string>>("c", "csv", "path of comma seperated CSV file compressed by zstd.", "", &filein);
  auto bgenfile = opts.add<Value<std::string>>("g", "bgen", "path of BGEN file compressed by gzip/zstd.", "", &filein);
  auto beaglefile = opts.add<Value<std::string>>("G", "beagle", "path of BEAGLE file compressed by gzip.", "", &filein);
  opts.add<Value<std::string>>("F", "match-bim", "the .mbim file to be matched, where the 7th column is allele frequency.", "", &filebim);
  auto usvprefix = opts.add<Value<std::string>>("P", "USV", "prefix of PCAone .eigvecs/.sigvals/.loadings/.mbim.");
  opts.add<Value<std::string>, Attribute::hidden>("", "read-U", "path of file with left singular vectors (.eigvecs).", "", &fileU);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-V", "path of file with right singular vectors (.loadings).", "", &fileV);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-S", "path of file with sigular values (.sigvals).", "", &fileS);
  
  opts.add<Value<std::string>, Attribute::headline>("","OUTPUT","Output options:");
  opts.add<Value<std::string>>("o", "out", "prefix of output files. default [pcaone].", fileout, &fileout);
  opts.add<Switch>("V", "printv", "output the right eigenvectors with suffix .loadings.", &printv);
  auto ld_opt = opts.add<Switch, Attribute::hidden>("D", "ld", "removed. LD no longer needs a residual matrix.");
  opts.add<Switch>("R", "print-r2", "print LD R2 to *.ld.gz file for pairwise SNPs within a window controlled by --ld-bp.", &print_r2);
  
  opts.add<Value<std::string>, Attribute::headline>("","MISC","Misc options:");
  auto maf_opt = opts.add<Value<double>>("", "maf", "exclude variants with MAF lower than this value. default is 0.05 for\n"
                                         "BEAGLE input, as in PCAngsd, and 0 (no filter) otherwise", maf, &maf);
  opts.add<Value<int>>("", "project", "project the new samples onto the existing PCs. Options are\n"
                                      "0: disabled;\n"
                                      "1: by multiplying the loadings with mean imputation for missing genotypes;\n"
                                      "2: by solving the least squares system Vx=g. skip sites with missingness;\n"
                                      "3: by EM to account for genotype uncertainty (BEAGLE input);\n"
                                      "4: by Augmentation, Decomposition and Procrusters transformation.\n", project, &project);
  opts.add<Unsigned>("", "project-bootstrap", "run SNP bootstrap diagnostics for --project 2 using this many replicates.", project_bootstrap, &project_bootstrap);
  opts.add<Switch>("", "project-bootstrap-save", "save raw bootstrap projection coordinates to *.proj.bootstrap.eigvecs.", &project_bootstrap_save);
  opts.add<Value<int>>("", "inbreed", "compute the inbreeding coefficient accounting for population structure. Options are\n"
                                      "0: disabled;\n"
                                      "1: compute per-site inbreeding coefficient and HWE test.\n", inbreed, &inbreed);
  opts.add<Switch>("", "evaladmix", "compute the correlation of residuals (evalAdmix) given the top PCs.", &evaladmix);
  opts.add<Value<int>>("", "evaladmix-k", "number of PCs used by --evaladmix. default is all computed PCs (use K-1 for an admixture model with K populations).", evaladmix_k, &evaladmix_k);
  opts.add<Value<int>>("", "selection", "compute selection statistics. Options are\n"
                                      "0: disabled;\n"
                                      "1: perform selection scan using Galinsky et al method;\n"
                                      "2: perform selection scan using PCAdapt method.\n", selection, &selection);
  opts.add<Value<double>>("", "ld-r2", "R2 cutoff for LD-based pruning (usually 0.2).", ld_r2, &ld_r2);
  opts.add<Unsigned>("", "ld-bp", "physical distance threshold in bases for LD window.", ld_bp, &ld_bp);
  opts.add<Value<int>>("", "ld-stats", "statistics to compute LD R2 for pairwise SNPs. Options are\n"
                                       "0: the ancestry adjusted, i.e. correlation between the residuals after\n"
                                       "   removing the PCs given by -P/--USV (reads .eigvecs of the same samples);\n"
                                       "1: the standard, i.e. correlation between two alleles.\n", ld_stats, &ld_stats);
  auto clumpfile = opts.add<Value<std::string>>("", "clump", "assoc-like file with target variants and pvalues for clumping.", "", &clump);
  auto assocnames = opts.add<Value<std::string>>("", "clump-names", "column names in assoc-like file for locating chr, pos and pvalue.", "CHR,BP,P", &assoc_colnames);
  opts.add<Value<double>>("", "clump-p1", "significance threshold for index SNPs.", clump_p1, &clump_p1);
  opts.add<Value<double>>("", "clump-p2", "secondary significance threshold for clumped SNPs.", clump_p2, &clump_p2);
  opts.add<Value<double>>("", "clump-r2", "r2 cutoff for LD-based clumping.", clump_r2, &clump_r2);
  opts.add<Unsigned>("", "clump-bp", "physical distance threshold in bases for clumping.", clump_bp, &clump_bp);
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
      for (const auto& uo : opts.unknown_options()) std::cerr << "unknown option: " << uo << "\n";
      exit(EXIT_FAILURE);
    }
    if (!opts.non_option_args().empty()) {
      for (const auto& a : opts.non_option_args()) std::cerr << "unexpected argument: " << a << "\n";
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
      throw std::invalid_argument("-d/--svd supports only 0, 1, 2 or 3");

    if (binfile->is_set() || ld_opt->is_set())
      throw std::invalid_argument(
          "-B/--binary, -D/--ld and the .residuals file were removed in v0.8.0. run the LD analysis on the genotypes "
          "directly, removing the PCs of a previous run:\n"
          "  PCAone -b plink -k 2 -o pcs\n"
          "  PCAone -b plink -P pcs --ld-r2 0.8 --ld-bp 1000000 -o adj\n"
          "-D/--ld computed its PCs without standardizing the sites; add --scale 0 to the first run to do the same");

    // all the inputs write to filein, so a second one would silently replace the first
    const int ninputs =
        plinkfile->is_set() + pgenfile->is_set() + csvfile->is_set() + bgenfile->is_set() + beaglefile->is_set();
    require(ninputs <= 1, "please give only one of -b/--bfile, -p/--pgen, -c/--csv, -g/--bgen and -G/--beagle");

    if (plinkfile->is_set())
      file_t = FileType::PLINK;
    else if (bgenfile->is_set())
      file_t = FileType::BGEN;
    else if (beaglefile->is_set())
      file_t = FileType::BEAGLE;
    else if (csvfile->is_set())
      file_t = FileType::CSV;
    else if (pgenfile->is_set())
      file_t = FileType::PGEN;
    else if (help_opt->is_set()) {
      std::cout << opts.help(Attribute::advanced) << "\n";
      exit(EXIT_SUCCESS);
    } else if (argc == 1) {
      std::cout << opts << "\n";
      exit(EXIT_SUCCESS);
    } else {
      throw std::invalid_argument(
          "no input file. please give one of -b/--bfile, -p/--pgen, -c/--csv, -g/--bgen or -G/--beagle");
    }
    genetic = (file_t == FileType::PLINK || file_t == FileType::BGEN || file_t == FileType::PGEN);
    // guard the option values; popl checks only their type
    require(threads >= 1, "-n/--threads must be at least 1");
    require(verbose <= 3, "-v/--verbose supports only 0, 1, 2 or 3");
    require(memory >= 0, "-m/--memory must be >= 0 (0 for in-core mode)");
    require(k >= 1, "-k/--pc must be at least 1");
    require(scale == SCALE_STANDARDIZE_GENETIC || (scale >= 0 && scale <= 4),
            "-C/--scale supports only -9, 0, 1, 2, 3 or 4");
    require(maxp >= 1, "--maxp must be at least 1");
    // the window-based RSVD doubles the band size every epoch up to -w (Halko.cpp)
    require(bands >= 4 && (bands & (bands - 1)) == 0, "-w/--batches must be a power of 2 and at least 4");
    require(scaleFactor > 0, "--scale-factor must be > 0");
    require(buffer >= 1, "--buffer must be at least 1 (GB)");
    require(imaxiter >= 1, "--imaxiter must be at least 1");
    require(itol > 0, "--itol must be > 0");
    require(rand <= 1, "--rand supports only 0 (uniform) or 1 (gaussian)");
    require(tol > 0, "--tol-rsvd must be > 0");
    require(tolem > 0, "--tol-em must be > 0");
    require(tolmaf > 0, "--tol-maf must be > 0");
    require(maf >= 0 && maf < 0.5, "--maf has to be in [0, 0.5); 0 disables the filter");
    require(project >= 0 && project <= 3, "--project supports only 0, 1, 2 or 3 in this release");
    require(project_bootstrap != 1, "--project-bootstrap needs at least 2 replicates");
    require(!project_bootstrap_save || project_bootstrap > 0, "--project-bootstrap-save requires --project-bootstrap");
    require(inbreed == 0 || inbreed == 1, "--inbreed supports only 0 or 1");
    require(selection >= 0 && selection <= 2, "--selection supports only 0, 1 or 2");
    require(evaladmix_k >= 0, "--evaladmix-k must be >= 0 (0 uses all computed PCs)");
    require(evaladmix_k == 0 || evaladmix, "--evaladmix-k requires --evaladmix");
    require(evaladmix_k <= (int)k, "--evaladmix-k cannot be larger than -k/--pc");
    require(ld_r2 >= 0 && ld_r2 <= 1, "--ld-r2 has to be in [0, 1]; 0 disables the pruning");
    require(ld_bp >= 1, "--ld-bp must be at least 1");
    require(ld_stats == 0 || ld_stats == 1, "--ld-stats supports only 0 or 1");
    require(clump_p1 > 0 && clump_p1 <= 1, "--clump-p1 has to be in (0, 1]");
    require(clump_p2 > 0 && clump_p2 <= 1, "--clump-p2 has to be in (0, 1]");
    // only the SNPs with p <= --clump-p2 are read (map_index_snps), so a larger
    // --clump-p1 would silently act as --clump-p2
    require(clump_p1 <= clump_p2, "--clump-p1 cannot be larger than --clump-p2");
    require(clump_r2 > 0 && clump_r2 <= 1, "--clump-r2 has to be in (0, 1]");
    require(clump_bp >= 1, "--clump-bp must be at least 1");
    require(three_names(assoc_colnames),
            "--clump-names needs 3 comma-separated column names for chr, pos and pvalue, e.g. CHR,BP,P");

    // handle PI, i.e U,S,V
    if (usvprefix->is_set()) {
      if (fileU.empty()) fileU = usvprefix->value() + ".eigvecs";
      if (fileE.empty()) fileE = usvprefix->value() + ".eigvals";
      if (fileS.empty()) fileS = usvprefix->value() + ".sigvals";
      if (fileV.empty()) fileV = usvprefix->value() + ".loadings";
      if (filebim.empty()) filebim = usvprefix->value() + ".mbim";
    }

    // handle LD. no PCA is run: the genotypes are read again and the PCs of
    // -P/--USV are removed from them as they are read (see run_ld_stuff).
    // dopca stays on, so the allele frequencies are estimated from these
    // genotypes and --maf filters them, as in a PCA run.
    if (print_r2 || ld_r2 > 0 || !clump.empty()) {
      ld = true;
      if (file_t != FileType::PLINK && file_t != FileType::PGEN)
        throw std::invalid_argument("--print-r2, --ld-r2 and --clump support only --bfile/--pgen input");
      if (ld_stats == 0 && fileU.empty())
        throw std::invalid_argument(
            "the ancestry adjusted LD (--ld-stats 0, the default) removes the PCs of a previous run of the same "
            "samples. please give its prefix with -P/--USV, or use --ld-stats 1 for the standard LD");
      memory /= 2.0;  // two blocks of genotypes are held at a time (LDColumns)
    }

    // input types each analysis can read. checked here, before anything is read,
    // because the readers that lack a mode crash or silently do something else:
    // a reader without P (genotype likelihoods) segfaults in the PCAngsd update,
    // one without C (missingness) segfaults in the EMU update, and a reader whose
    // read_block_update() is empty returns the last block for every block.
    const bool plink_or_pgen = (file_t == FileType::PLINK || file_t == FileType::PGEN);
    if (project > 0) {
      require(plink_or_pgen || file_t == FileType::BEAGLE,
              "--project supports only --bfile, --pgen and --beagle input");
      require(project != 3 || file_t == FileType::BEAGLE,
              "--project 3 requires BEAGLE genotype likelihood input (-G/--beagle)");
    }
    require(inbreed == 0 || plink_or_pgen || file_t == FileType::BEAGLE,
            "--inbreed supports only --bfile, --pgen and --beagle input");
    require(!evaladmix || plink_or_pgen, "--evaladmix supports only --bfile and --pgen input");
    require(!pcangsd || file_t == FileType::PLINK || file_t == FileType::BEAGLE,
            "--pcangsd supports only --beagle (genotype likelihoods) and --bfile input");
    require(!emu || file_t != FileType::CSV, "--emu supports only --bfile, --pgen and --bgen input");
    require(!(emu && file_t == FileType::BGEN && memory > 0),
            "--emu with --bgen is not supported with -m (out-of-core) yet. please run it in-core");

    // handle projection
    if (project > 0) {
      if (fileV.empty() || fileS.empty()) throw std::invalid_argument("please use --USV together with --project");
      if (project_bootstrap > 0 && project != 2)
        throw std::invalid_argument("--project-bootstrap currently supports only --project 2");
      dopca = false, missme = true, out_of_core = false;
      memory = 0;
    } else if (project_bootstrap > 0 || project_bootstrap_save) {
      throw std::invalid_argument("--project-bootstrap requires --project 2");
    }

    // handle selection
    if (selection > 0) {
      if (file_t != FileType::PLINK && file_t != FileType::PGEN)
        throw std::invalid_argument("only supports --bfile/--pgen for now");
      if (fileU.empty() || fileE.empty()) throw std::invalid_argument("please use --USV together with --selection");
      dopca = true;  // we need this to init F
    }

    // handle inbreeding
    if (inbreed > 0) {
      dopca = false, center = false;
      if (fileU.empty() || fileV.empty() || fileS.empty())
        throw std::invalid_argument("please use --USV together with --inbreed");
    }

    // handle memory and misc options
    // --ncv was always overwritten before. Spectra needs k < ncv
    if (!ncv_opt->is_set())
      ncv = 20 > (2 * k + 1) ? 20 : (2 * k + 1);
    else
      require(ncv > k, "--ncv must be greater than -k/--pc");
    oversamples = oversamples > k ? oversamples : k;
    if (haploid && genetic) ploidy = 1;
    if (memory > 0 && svd_t != SvdType::FULL) out_of_core = true;

    // PCAngsd filters MAF < 0.05 by default. Its covariance divides by 2f(1-f)
    // per site, so rare sites from genotype likelihoods dominate it, and a site
    // whose EM frequency is 0 made the whole .cov diagonal inf.
    if (file_t == FileType::BEAGLE && dopca && !out_of_core && !maf_opt->is_set()) maf = 0.05;  // -m is refused below
    filterSNP = maf > 0 ? true : false;  // filter SNP if MAf applied
    if (filterSNP) {
      if (out_of_core) throw std::invalid_argument("does not support --maf filters for out-of-core mode yet! ");
    }

    // handle EM-PCA
    if (dopca && file_t == FileType::BEAGLE) pcangsd = true;
    // both would run their own EM update on the same matrix (Data::fit_with_pi)
    require(!(emu && pcangsd), "--emu cannot be used with --pcangsd or BEAGLE input (which implies --pcangsd)");
    if (emu || pcangsd) {
      missme = true;
    } else if (dopca) {
      maxiter = 0;
    }
    // The full SVD (--svd 3) is a single exact eigendecomposition with no EM
    // loop around it, so it cannot fit individual allele frequencies. Refuse
    // the combination instead of returning mean-imputed PCA in files that are
    // indistinguishable from real EM-PCA output.
    if (missme && svd_t == SvdType::FULL)
      throw std::invalid_argument(
          "EM-PCA is not supported with --svd 3 (full SVD). please use --svd 0, 1 or 2 instead. note that "
          "--emu, --pcangsd and BEAGLE input (which implies --pcangsd) all request EM-PCA");
    if (out_of_core && pcangsd && (file_t == FileType::BEAGLE))
      throw std::invalid_argument("not supporting -m option (out-of-core) for PCAngsd and BEAGLE input yet!");

    // Shuffling is for the winSVD of a PCA run only. LD walks the sites in
    // .bim/.pvar order, and --inbreed, --project and --selection decompose
    // nothing: they read a reference PCA and pass over the genotypes in file
    // order. For PGEN the permutation is logical and is only built for a PCA
    // run (Main.cpp), so leaving the flag on made those reads index an empty
    // permutation -- out-of-core --inbreed on PGEN segfaulted.
    if (svd_t == SvdType::PCAoneAlg2 && !noshuffle && !ld && inbreed == 0 && project == 0 && selection == 0)
      perm = true;

  } catch (const popl::invalid_option& e) {
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
  } catch (const std::exception& e) {
    std::cerr << "Exception: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
}

Param::~Param() {}
