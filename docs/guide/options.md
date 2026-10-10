# Command-line options

Run `PCAone` without arguments to print the common options below, `PCAone --help` to add the
advanced ones, or `PCAone --groff > pcaone.1 && man ./pcaone.1` for a man page.

```text
General options:
  -h, --help                     print all options, including the advanced ones
  -m, --memory arg (=0)          memory in GB for out-of-core mode; 0 runs in-core. also sets the blocks
                                 of --svd 3 and the stripes of --evaladmix-kin
  -n, --threads arg (=12)        number of threads
  -v, --verbose arg (=1)         verbosity: 0 silent, 1 concise, 2 verbose, 3 debug

PCA methods:
  -d, --svd arg (=2)             PCA method:
                                 0: IRAM, the implicitly restarted Arnoldi method;
                                 1: sSVD, single-pass randomized SVD with power iterations;
                                 2: winSVD, window-based randomized SVD, for large data;
                                 3: exact PCA from the N x N sample GRM, for small N (no EM-PCA)
  -k, --pc arg (=10)             number of PCs. with -P/--USV, the number of leading PCs of the
                                 reference to use (default: all)
  -C, --scale arg (=-9)          normalization:
                                 -9: standardize genotypes by sqrt(ploidy*f*(1-f));
                                  0: none;
                                  1: center and scale each feature, as R's scale();
                                  2: count per median and log (CPMED), then standardize;
                                  3: log1p, then standardize;
                                  4: relative counts, then standardize
  --maxp arg (=20)               maximum number of power iterations of --svd 1 and 2
  -S, --no-shuffle               do not shuffle the features for --svd 2, e.g. for gene counts
  --seed arg (=112)              random seed
  --emu                          EM-PCA of EMU, for genotypes with missing calls
  --pcangsd                      EM-PCA of PCAngsd, for genotype likelihoods; implied by -G
  --em-k arg                     number of PCs that model the allele frequencies in --emu and
                                 --pcangsd; -k PCs are written (default: -k)

Input options:
  -b, --bfile arg                prefix of PLINK .bed/.bim/.fam
  -p, --pgen arg                 prefix of PLINK2 .pgen/.pvar/.psam
  -g, --bgen arg                 BGEN file (layout 1 or 2)
  -G, --beagle arg               BEAGLE genotype likelihoods, gzip-compressed
  -c, --csv arg                  comma-separated values, zstd-compressed
  -P, --USV arg                  prefix of a previous PCAone run (.eigvecs, .sigvals, .loadings, .mbim)
  -F, --match-bim arg            the .mbim to match the variants with (allele frequencies in column 7)
  --maf arg (=0)                 exclude variants with MAF below this (default: 0.05 for -G, else 0)

Output options:
  -o, --out arg (=pcaone)        prefix of the output files
  --no-loadings                  do not write the .loadings and .mbim (skips the second pass of --svd 3)

Analyses of a previous PCA, given with -P/--USV:
  --project arg (=0)             project new samples onto the PCs:
                                 0: off;
                                 1: multiply by the loadings, missing calls at the mean;
                                 2: least squares on the called sites;
                                 3: EM over the genotype likelihoods of -G
  --project-bootstrap arg (=0)   number of SNP bootstrap replicates of --project 2
  --project-bootstrap-save       also write the replicates to .proj.bootstrap.eigvecs
  --inbreed arg (=0)             inbreeding under population structure:
                                 0: off;
                                 1: per-site F and HWE test (.hwe);
                                 2: per-sample F (.inbred)
  --evaladmix                    kinship from the correlation of residuals (evalAdmix): .kinship and .corres.
                                 -P must hold the same samples in the same order
  --evaladmix-kin arg            with --evaladmix, write only the pairs with kinship >= this (.kin0) and an
                                 unrelated set (.unrelated), within -m. for biobanks, e.g. 0.0442 (3rd degree)
  --evaladmix-unrelated arg      kinship cutoff of the .unrelated set, at least that of --evaladmix-kin
                                 (default: the same)
  --evaladmix-ibd                with --evaladmix, also estimate the IBD sharing k0, k1, k2 (.k0, .k2, or .kin0 columns)
  --selection arg (=0)           selection scan along the PCs:
                                 0: off;
                                 1: Galinsky et al. (FastPCA);
                                 2: pcadapt

LD, ancestry-adjusted with -P/--USV:
  -R, --print-r2                 write the R2 of the SNP pairs within --ld-bp to .ld.gz
  --ld-r2 arg (=0)               prune to R2 below this cutoff, e.g. 0.2 (.ld.prune.in, .ld.prune.out)
  --ld-bp arg (=1000000)         LD window in bases
  --ld-stats arg (=0)            LD statistic:
                                 0: ancestry-adjusted, the R2 of the residuals of the PCs of -P/--USV;
                                 1: standard
  --clump arg                    association files to clump, comma-separated
  --clump-names arg (=CHR,BP,P)  columns of the chromosome, position and p-value
  --clump-p1 arg (=0.0001)       p-value cutoff of the index variants
  --clump-p2 arg (=0.01)         p-value cutoff of the clumped variants
  --clump-r2 arg (=0.5)          R2 cutoff of clumping
  --clump-bp arg (=250000)       clumping window in bases
```
