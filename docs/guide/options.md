# Command-line options

Run `PCAone` without arguments to print the common options below, `PCAone --help` to add the
advanced ones, or `PCAone --groff > pcaone.1 && man ./pcaone.1` for a man page.

```text
General options:
  -h, --help                     print all options including hidden advanced options
  -m, --memory arg (=0)          RAM usage in GB unit for out-of-core mode. default is in-core mode.
                                 with --svd 3, it sets the blocks the GRM is streamed in;
                                 with --evaladmix-kin, the stripes of the N x N matrix
  -n, --threads arg (=12)        the number of threads to be used
  -v, --verbose arg (=1)         verbosity level for logs. Options are
                                 0: silent, no messages on screen;
                                 1: concise messages to screen;
                                 2: more verbose information;
                                 3: enable debug information.

PCA algorithms:
  -d, --svd arg (=2)             SVD method to be applied. default 2 is recommended for big data. Options are
                                 0: the Implicitly Restarted Arnoldi Method (IRAM);
                                 1: the Yu's single-pass Randomized SVD with power iterations;
                                 2: the accurate window-based Randomized SVD method (PCAone);
                                 3: exact PCA by eigendecomposition of the sample GRM, streamed block by block when N <= M,
                                    in N x N memory (no EM-PCA support).
  -k, --pc arg (=10)             top k principal components (PCs) to be calculated. with -P/--USV, the number
                                 of leading PCs of the reference to use (default all)
  -C, --scale arg (=-9)          do normalization or scaling for input file. Options are
                                 -9: standardize genetic data by sqrt(ploidy*f*(1-f));
                                  0: do nothing and proceed to SVD;
                                  1: do direct standardization, as the scale(x, center=TRUE, scale=TRUE) function in R;
                                  2: do first count per median log transformation (CPMED), then standardization;
                                  3: do first log1p transformation, then standardization;
                                  4: do first relative counts, then standardization.
  --maxp arg (=20)               maximum number of power iterations for RSVD algorithm.
  -S, --no-shuffle               do not shuffle columns of data for --svd 2 (if not locally correlated).
  --seed arg (=112)              seeds for reproducing results.
  --emu                          use EMU algorithm for genotype input with missingness. not with --svd 3.
  --pcangsd                      use PCAngsd algorithm for genotype likelihood input. not with --svd 3.
  --em-k arg                     the number of PCs that model the individual allele frequencies in the EM iterations
                                 of --emu and --pcangsd. -k PCs of the final matrix are written. default is -k

Input options:
  -b, --bfile arg                prefix of PLINK .bed/.bim/.fam files.
  -p, --pgen arg                 prefix of PLINK2 .pgen/.pvar/.psam files.
  -c, --csv arg                  path of comma seperated CSV file compressed by zstd.
  -g, --bgen arg                 path of BGEN file (layout 1 or 2; zlib, zstd or no compression).
  -G, --beagle arg               path of BEAGLE file compressed by gzip.
  -F, --match-bim arg            the .mbim file to be matched, where the 7th column is allele frequency.
  -P, --USV arg                  prefix of PCAone .eigvecs/.sigvals/.loadings/.mbim.

Output options:
  -o, --out arg (=pcaone)        prefix of output files. default [pcaone].
  --no-loadings                  do not output the right eigenvectors (.loadings) and .mbim, which are written by default.
                                 with --svd 3, this also skips the second pass over the data
  -R, --print-r2                 print LD R2 to *.ld.gz file for pairwise SNPs within a window controlled by --ld-bp.

Misc options:
  --maf arg (=0)                 exclude variants with MAF lower than this value. default is 0.05 for
                                 BEAGLE input, as in PCAngsd, and 0 (no filter) otherwise
  --project arg (=0)             project the new samples onto the existing PCs. Options are
                                 0: disabled;
                                 1: by multiplying the loadings with mean imputation for missing genotypes;
                                 2: by solving the least squares system Vx=g. skip sites with missingness;
                                 3: by EM to account for genotype uncertainty (BEAGLE input).
  --project-bootstrap arg (=0)   run SNP bootstrap diagnostics for --project 2 using this many replicates.
  --project-bootstrap-save       save raw bootstrap projection coordinates to *.proj.bootstrap.eigvecs.
  --inbreed arg (=0)             compute the inbreeding coefficient accounting for population structure. Options are
                                 0: disabled;
                                 1: compute per-site inbreeding coefficient and HWE test;
                                 2: compute per-sample inbreeding coefficient.
  --evaladmix                    compute the correlation of residuals (evalAdmix) given existing PCs from -P/--USV
                                 (same samples in the same order).
  --evaladmix-kin arg            with --evaladmix, write only the pairs with kinship >= this cutoff to .kin0, and a
                                 maximal unrelated set to .unrelated, instead of the two N x N matrices. computes the
                                 matrix in stripes that fit in -m, for biobank-scale samples (e.g. 0.0442, 3rd degree)
  --evaladmix-unrelated arg      kinship cutoff for the .unrelated set of --evaladmix-kin, at least that cutoff.
                                 default is the --evaladmix-kin cutoff
  --evaladmix-ibd                with --evaladmix, also estimate the probabilities of sharing 0, 1 and 2 alleles IBD
                                 (k0, k1, k2): .k0 and .k2 matrices, or K0 K1 K2 columns in the .kin0 of --evaladmix-kin.
                                 one more Gram product of the size of the kinship one
  --selection arg (=0)           compute selection statistics. Options are
                                 0: disabled;
                                 1: perform selection scan using Galinsky et al method;
                                 2: perform selection scan using PCAdapt method.
  --ld-r2 arg (=0)               R2 cutoff for LD-based pruning (usually 0.2).
  --ld-bp arg (=1000000)         physical distance threshold in bases for LD window.
  --ld-stats arg (=0)            statistics to compute LD R2 for pairwise SNPs. Options are
                                 0: the ancestry adjusted, i.e. correlation between the residuals after
                                    removing the PCs given by -P/--USV (reads .eigvecs of the same samples);
                                 1: the standard, i.e. correlation between two alleles.
  --clump arg                    assoc-like file with target variants and pvalues for clumping.
  --clump-names arg (=CHR,BP,P)  column names in assoc-like file for locating chr, pos and pvalue.
  --clump-p1 arg (=0.0001)       significance threshold for index SNPs.
  --clump-p2 arg (=0.01)         secondary significance threshold for clumped SNPs.
  --clump-r2 arg (=0.5)          r2 cutoff for LD-based clumping.
  --clump-bp arg (=250000)       physical distance threshold in bases for clumping.
```
