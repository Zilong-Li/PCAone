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
  bool noloadings = false;
  std::string copyr{"PCA All In One (v" + (std::string)VERSION + ")        https://github.com/Zilong-Li/PCAone\n" +
                    "(C) 2021-2026 Zilong Li        GNU General Public License v3\n" +
                    "\n" +
                    "Usage: PCAone <input> [options]\n" +
                    "\n" +
                    "Common usage:\n" +
                    "  PCAone -b plink -k 10 -o pcs                  PCA of PLINK files with winSVD (default)\n" +
                    "  PCAone -b plink -k 10 -m 4 -o pcs             the same out-of-core, with 4 GB for the data blocks\n" +
                    "  PCAone -p plink2 -k 10 -m 4 -o pcs            PLINK2 PGEN, with dosages when present\n" +
                    "  PCAone -g data.bgen -k 10 -m 4 -o pcs         BGEN, e.g. imputed genotypes\n" +
                    "  PCAone -b plink -d 3 -k 10 -o pcs             exact PCA, fastest for a few thousand samples\n" +
                    "  PCAone -b plink --emu -k 10 -o pcs            EMU, for genotypes with many missing calls\n" +
                    "  PCAone -G data.beagle.gz -k 10 -o pcs         PCAngsd, for genotype likelihoods\n" +
                    "\n" +
                    "Then reuse the PCs in pcs.* with -P:\n" +
                    "  PCAone -b new -P pcs --project 2 -o new       project new samples onto the PCs\n" +
                    "  PCAone -b plink -P pcs -k 3 --ld-r2 0.2       prune on the LD left after the top 3 PCs\n" +
                    "  PCAone -b plink -P pcs -k 3 --evaladmix       kinship of every pair of samples\n" +
                    "  PCAone -b plink -P pcs -k 3 --inbreed 1       HWE test of each site under structure\n" +
                    "  PCAone -b plink -P pcs --selection 2          selection scan along the PCs\n" +
                    "\n" +
                    "Other data:\n" +
                    "  PCAone -c data.csv.zst -k 10 -C 2 -S          CSV counts, e.g. single-cell RNA-seq\n" +
                    "\n" +
                    "Common options are below; --help adds the advanced ones. Documentation: https://zilongli.org/PCAone\n"};
  OptionParser opts(copyr);
  opts.add<Value<std::string>, Attribute::headline>("","PCAone","General options:");
  auto help_opt = opts.add<Switch>("h", "help", "print all options, including the advanced ones");
  opts.add<Value<double>>("m", "memory", "memory in GB for out-of-core mode; 0 runs in-core. also sets the blocks\n"
                                        "of --svd 3 and the stripes of --evaladmix-kin", memory, &memory);
  opts.add<Unsigned>("n", "threads", "number of threads", threads, &threads);
  opts.add<Unsigned>("v", "verbose", "verbosity: 0 silent, 1 concise, 2 verbose, 3 debug", verbose, &verbose);

  opts.add<Value<std::string>, Attribute::headline>("","PCA","PCA methods:");
  auto svd_opt = opts.add<Unsigned>("d", "svd", "PCA method:\n"
                                                "0: IRAM, the implicitly restarted Arnoldi method;\n"
                                                "1: sSVD, single-pass randomized SVD with power iterations;\n"
                                                "2: winSVD, window-based randomized SVD, for large data;\n"
                                                "3: exact PCA from the N x N sample GRM, for small N (no EM-PCA)", 2);
  auto k_opt = opts.add<Unsigned>("k", "pc", "number of PCs. with -P/--USV, the number of leading PCs of the\n"
                                             "reference to use (default: all)", k, &k);
  opts.add<Value<int>>("C", "scale", "scaling after centering:\n"
                                     "-9: standardize genotypes by sqrt(ploidy*f*(1-f));\n"
                                     " 0: none;\n"
                                     " 1: center and scale each feature, as R's scale();\n"
                                     " 2: count per median and log (CPMED), then standardize;\n"
                                     " 3: log1p, then standardize;\n"
                                     " 4: relative counts, then standardize", scale,  &scale);
  opts.add<Unsigned>("", "maxp", "maximum number of power iterations of --svd 1 and 2", maxp, &maxp);
  opts.add<Switch>("S", "no-shuffle", "do not shuffle the features for --svd 2, e.g. for gene counts", &noshuffle);
  opts.add<Unsigned, Attribute::advanced>("w", "batches", "number of mini-batches of --svd 2", bands, &bands);
  opts.add<Value<int>>("", "seed", "random seed", seed, &seed);
  opts.add<Switch>("", "emu", "EM-PCA of EMU, for genotypes with missing calls", &emu);
  opts.add<Switch>("", "pcangsd", "EM-PCA of PCAngsd, for genotype likelihoods; implied by -G", &pcangsd);
  auto em_k_opt = opts.add<Unsigned>("", "em-k", "number of PCs that model the allele frequencies in --emu and --pcangsd\n");
  opts.add<Unsigned, Attribute::advanced>("", "M", "number of features (e.g. SNPs), if known", 0, &nsnps);
  opts.add<Unsigned, Attribute::advanced>("", "N", "number of samples, if known", 0, &nsamples);
  opts.add<Value<double>, Attribute::advanced>("", "scale-factor", "multiply the normalized counts of each sample by this value", 1.0, &scaleFactor);
  opts.add<Switch, Attribute::advanced>("", "no-prefetch", "out-of-core: do not read the next block during the computation", &noprefetch);
  opts.add<Switch, Attribute::advanced>("", "bed-copy", "out-of-core --svd 2 on BED: write and read a shuffled copy, <out>.perm.*;\n"
                                                        "faster for a BED larger than the memory on a spinning disk", &bedcopy);
  opts.add<Unsigned, Attribute::advanced>("", "buffer", "memory in GiB for the genotypes when shuffling to a copy", buffer, &buffer);
  opts.add<Unsigned, Attribute::advanced>("", "imaxiter", "maximum number of IRAM iterations", imaxiter, &imaxiter);
  opts.add<Value<double>, Attribute::advanced>("", "itol", "tolerance of IRAM", itol, &itol);
  auto ncv_opt = opts.add<Unsigned, Attribute::advanced>("", "ncv", "number of Lanczos vectors of IRAM", ncv, &ncv);
  opts.add<Unsigned, Attribute::advanced>("", "oversamples", "number of oversampling columns of the RSVD", oversamples, &oversamples);
  opts.add<Unsigned, Attribute::advanced>("", "rand", "random matrix of the RSVD: 0 uniform, 1 Gaussian", rand, &rand);
  opts.add<Unsigned, Attribute::advanced>("", "maxiter", "maximum number of EM iterations", maxiter, &maxiter);
  opts.add<Value<double>, Attribute::advanced>("", "tol-rsvd", "tolerance of the RSVD", tol, &tol);
  opts.add<Value<double>, Attribute::advanced>("", "tol-em", "tolerance of the EM iterations", tolem, &tolem);
  opts.add<Value<double>, Attribute::advanced>("", "tol-maf", "tolerance of the EM for allele frequencies", tolmaf, &tolmaf);

  opts.add<Value<std::string>, Attribute::headline>("","INPUT","Input options:");
  auto plinkfile = opts.add<Value<std::string>>("b", "bfile", "prefix of PLINK .bed/.bim/.fam", "", &filein);
  opts.add<Switch, Attribute::advanced>("", "haploid", "the PLINK files hold haploid data", &haploid);
  auto pgenfile = opts.add<Value<std::string>>("p", "pgen", "prefix of PLINK2 .pgen/.pvar/.psam", "", &filein);
  opts.add<Switch, Attribute::advanced>("", "hardcall", "use the hard calls of a PGEN instead of its dosages", &hardcall);
  // removed in v0.8.0; kept hidden only to say what replaces them
  auto binfile = opts.add<Value<std::string>, Attribute::hidden>("B", "binary", "removed. LD now reads the genotypes directly.");
  auto bgenfile = opts.add<Value<std::string>>("g", "bgen", "BGEN file (layout 1 or 2)", "", &filein);
  auto beaglefile = opts.add<Value<std::string>>("G", "beagle", "BEAGLE genotype likelihoods, gzip-compressed", "", &filein);
  auto csvfile = opts.add<Value<std::string>>("c", "csv", "comma-separated values, zstd-compressed", "", &filein);
  auto usvprefix = opts.add<Value<std::string>>("P", "USV", "prefix of a previous PCAone run (.eigvecs, .sigvals, .loadings, .mbim)");
  opts.add<Value<std::string>>("F", "match-bim", "the .mbim to match the variants with (allele frequencies in column 7)", "", &filebim);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-U", "path of file with left singular vectors (.eigvecs).", "", &fileU);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-V", "path of file with right singular vectors (.loadings).", "", &fileV);
  opts.add<Value<std::string>, Attribute::hidden>("", "read-S", "path of file with sigular values (.sigvals).", "", &fileS);
  auto maf_opt = opts.add<Value<double>>("", "maf", "exclude variants with MAF below this (default: 0.05 for -G, else 0)", maf, &maf);

  opts.add<Value<std::string>, Attribute::headline>("","OUTPUT","Output options:");
  opts.add<Value<std::string>>("o", "out", "prefix of the output files", fileout, &fileout);
  // removed in v0.8.0, when the .loadings became the default; kept hidden only to say so
  auto printv_opt = opts.add<Switch, Attribute::hidden>("V", "printv", "removed. the .loadings are written by default.");
  opts.add<Switch>("", "no-loadings", "do not write the .loadings and .mbim (skips the second pass of --svd 3)", &noloadings);
  auto ld_opt = opts.add<Switch, Attribute::hidden>("D", "ld", "removed. LD no longer needs a residual matrix.");

  opts.add<Value<std::string>, Attribute::headline>("","MISC","Analyses of a previous PCA, given with -P/--USV:");
  opts.add<Value<int>>("", "project", "project new samples onto the PCs:\n"
                                      "0: off;\n"
                                      "1: multiply by the loadings, missing calls at the mean;\n"
                                      "2: least squares on the called sites;\n"
                                      "3: EM over the genotype likelihoods of -G", project, &project);
  opts.add<Unsigned>("", "project-bootstrap", "number of SNP bootstrap replicates of --project 2", project_bootstrap, &project_bootstrap);
  opts.add<Switch>("", "project-bootstrap-save", "also write the replicates to .proj.bootstrap.eigvecs", &project_bootstrap_save);
  opts.add<Value<int>>("", "inbreed", "inbreeding under population structure:\n"
                                      "0: off;\n"
                                      "1: per-site F and HWE test (.hwe);\n"
                                      "2: per-sample F (.inbred)", inbreed, &inbreed);
  opts.add<Switch>("", "evaladmix", "kinship from the correlation of residuals (evalAdmix): .kinship and .corres.\n"
                                    "-P must hold the same samples in the same order", &evaladmix);
  auto kin_opt = opts.add<Value<double>>("", "evaladmix-kin", "with --evaladmix, write only the pairs with kinship >= this (.kin0) and an\n"
                                         "unrelated set (.unrelated), within -m. for biobanks, e.g. 0.0442 (3rd degree)");
  auto unrel_opt = opts.add<Value<double>>("", "evaladmix-unrelated", "kinship cutoff of the .unrelated set, at least that of --evaladmix-kin\n"
                                           "(default: the same)");
  opts.add<Switch>("", "evaladmix-ibd", "with --evaladmix, also estimate the IBD sharing k0, k1, k2 (.k0, .k2, or .kin0 columns)", &evaladmix_ibd);
  opts.add<Value<int>>("", "selection", "selection scan along the PCs:\n"
                                        "0: off;\n"
                                        "1: Galinsky et al. (FastPCA);\n"
                                        "2: pcadapt", selection, &selection);

  opts.add<Value<std::string>, Attribute::headline>("","LD","LD, ancestry-adjusted with -P/--USV:");
  opts.add<Switch>("R", "print-r2", "write the R2 of the SNP pairs within --ld-bp to .ld.gz", &print_r2);
  opts.add<Value<double>>("", "ld-r2", "prune to R2 below this cutoff, e.g. 0.2 (.ld.prune.in, .ld.prune.out)", ld_r2, &ld_r2);
  opts.add<Unsigned>("", "ld-bp", "LD window in bases", ld_bp, &ld_bp);
  opts.add<Value<int>>("", "ld-stats", "LD statistic:\n"
                                       "0: ancestry-adjusted, the R2 of the residuals of the PCs of -P/--USV;\n"
                                       "1: standard", ld_stats, &ld_stats);
  auto clumpfile = opts.add<Value<std::string>>("", "clump", "association files to clump, comma-separated", "", &clump);
  auto assocnames = opts.add<Value<std::string>>("", "clump-names", "columns of the chromosome, position and p-value", "CHR,BP,P", &assoc_colnames);
  opts.add<Value<double>>("", "clump-p1", "p-value cutoff of the index variants", clump_p1, &clump_p1);
  opts.add<Value<double>>("", "clump-p2", "p-value cutoff of the clumped variants", clump_p2, &clump_p2);
  opts.add<Value<double>>("", "clump-r2", "R2 cutoff of clumping", clump_r2, &clump_r2);
  opts.add<Unsigned>("", "clump-bp", "clumping window in bases", clump_bp, &clump_bp);
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

    if (printv_opt->is_set())
      throw std::invalid_argument(
          "-V/--printv was removed in v0.8.0. the .loadings and .mbim are now written by default; drop -V/--printv, "
          "or use --no-loadings to turn them off");
    printv = !noloadings;

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
    // the two-stage analyses (-P/--USV: LD, --evaladmix, --project, --selection,
    // --inbreed) use the leading -k PCs of the reference, or all of them without
    // -k: its default of 10 is not a choice made for that reference (ref_pcs())
    if (k_opt->is_set()) ref_k = k;
    // EM-PCA (--emu, --pcangsd) models the individual allele frequencies with
    // --em-k PCs, and writes the leading -k PCs of the final matrix, as EMU's
    // --eig and --eig-out. Every other path fits and writes -k PCs: em_k = k.
    // BEAGLE input implies --pcangsd (see "handle EM-PCA" below).
    em_k = k;
    if (em_k_opt->is_set()) {
      require(emu || pcangsd || file_t == FileType::BEAGLE, "--em-k requires EM-PCA: --emu, --pcangsd or BEAGLE input");
      require(em_k_opt->value() >= 1, "--em-k must be at least 1");
      em_k = em_k_opt->value();
    }
    require(scale == SCALE_STANDARDIZE_GENETIC || (scale >= 0 && scale <= 4),
            "-C/--scale supports only -9, 0, 1, 2, 3 or 4");
    require(maxp >= 1, "--maxp must be at least 1");
    // the window-based RSVD doubles the band size every epoch up to -w (Halko.cpp)
    require(bands >= 4 && (bands & (bands - 1)) == 0, "-w/--batches must be a power of 2 and at least 4");
    require(scaleFactor > 0, "--scale-factor must be > 0");
    require(buffer >= 1, "--buffer must be at least 1 (GiB)");
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
    require(inbreed >= 0 && inbreed <= 2, "--inbreed supports only 0, 1 or 2");
    require(selection >= 0 && selection <= 2, "--selection supports only 0, 1 or 2");
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
    // --inbreed pairs the target with the sites of the reference run, and no
    // reader filters them here, so --maf was silently ignored
    require(inbreed == 0 || !maf_opt->is_set(),
            "--inbreed uses the sites of the reference run. apply --maf in that run instead");
    require(inbreed == 0 || !haploid,
            "--inbreed models the heterozygosity of diploid genotypes, so it cannot be used with --haploid");
    require(!evaladmix || plink_or_pgen, "--evaladmix supports only --bfile and --pgen input");
    // Each analysis has its own early-return path in Main.cpp.
    require(!evaladmix || (project == 0 && selection == 0 && inbreed == 0 && !ld),
            "--evaladmix cannot be combined with --project, --selection, "
            "--inbreed, --print-r2, --ld-r2 or --clump");
    require(!(evaladmix && haploid),
            "--evaladmix takes each sample's variance from its heterozygosity, so it needs diploid genotypes");
    if (kin_opt->is_set()) {
      require(evaladmix, "--evaladmix-kin requires --evaladmix");
      evaladmix_pairs = true;
      evaladmix_kin = kin_opt->value();
      require(evaladmix_kin >= -0.5 && evaladmix_kin <= 0.5, "--evaladmix-kin is a kinship cutoff, in [-0.5, 0.5]");
      evaladmix_unrel = evaladmix_kin;
    }
    require(!evaladmix_ibd || evaladmix, "--evaladmix-ibd requires --evaladmix");
    if (unrel_opt->is_set()) {
      require(evaladmix_pairs, "--evaladmix-unrelated requires --evaladmix-kin");
      evaladmix_unrel = unrel_opt->value();
      // the .unrelated set is found among the pairs that are written
      require(evaladmix_unrel >= evaladmix_kin && evaladmix_unrel <= 0.5,
              "--evaladmix-unrelated has to be in [--evaladmix-kin, 0.5]");
    }
    require(!pcangsd || file_t == FileType::PLINK || file_t == FileType::BEAGLE,
            "--pcangsd supports only --beagle (genotype likelihoods) and --bfile input");
    require(!emu || file_t != FileType::CSV, "--emu supports only --bfile, --pgen and --bgen input");
    require(!(emu && file_t == FileType::BGEN && memory > 0),
            "--emu with --bgen is not supported with -m (out-of-core) yet. please run it in-core");

    // evalAdmix uses existing scores, and estimates frequencies from this input.
    // Keep centring on for mean imputation of missing calls; neither reader
    // standardizes here. No EM masks/probabilities or PCA permutation are needed.
    if (evaladmix) {
      require(!fileU.empty(), "please use -P/--USV or --read-U together with --evaladmix");
      dopca = true;
      center = true;
      emu = pcangsd = missme = false;
      em_k = k;
    }

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
    // --ncv was always overwritten before. Spectra needs k < ncv, for the EM
    // decomposition (--em-k PCs) and the final one (-k PCs) alike
    if (!ncv_opt->is_set())
      ncv = 20 > (2 * max_k() + 1) ? 20 : (2 * max_k() + 1);
    else
      require(ncv > max_k(), em_k > k ? "--ncv must be greater than --em-k" : "--ncv must be greater than -k/--pc");
    // each RSVD oversamples by at least its own rank: the EM one by --em-k
    em_oversamples = oversamples > em_k ? oversamples : em_k;
    oversamples = oversamples > k ? oversamples : k;
    if (haploid && genetic) ploidy = 1;
    if (memory > 0 && (svd_t != SvdType::FULL || evaladmix)) out_of_core = true;

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
    // .bim/.pvar order; evalAdmix, inbreeding, projection and selection decompose
    // nothing: they read a reference PCA and pass over the genotypes in file
    // order. For PGEN the permutation is logical and is only built for a PCA
    // run (Main.cpp), so leaving the flag on made those reads index an empty
    // permutation -- out-of-core --inbreed on PGEN segfaulted.
    if (svd_t == SvdType::PCAoneAlg2 && !noshuffle && !ld && !evaladmix && inbreed == 0 && project == 0 &&
        selection == 0)
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
