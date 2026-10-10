# Change log

## Legacy command migration

### LD workflows before v0.8.0

The `-D/--ld` and `-B/--binary` options and the intermediate `.residuals`
file have been removed. Compute a reference PCA, then pass its prefix with
`-P/--USV` while reading the original PLINK or PGEN genotypes:

```shell
PCAone -b data -k 3 -o ref
PCAone -b data -P ref --ld-r2 0.2 -o adj
```

The old `-D/--ld` workflow did not standardize sites. Add `--scale 0` to the
reference PCA command to reproduce that scaling. See v0.8.0 below for
result changes and compatibility details.

### Loadings before v0.8.0

`-V/--printv` has been removed, and PCAone stops if it is given. The `.loadings`
and `.mbim` are now written by default, so drop `-V/--printv` from old commands;
add `--no-loadings` to skip them.

### Projection syntax before v0.4.8

Older examples supplied `--read-V`, `--read-S`, and `--match-bim` separately.
Since v0.4.8, `-P/--USV` accepts the reference prefix for `.eigvecs`,
`.sigvals`, `.loadings`, and `.mbim`.

### Variance explained with `--svd 3` before v0.7.0

Older examples ran `--svd 3` to get the proportion of variance explained,
dividing each eigenvalue by the sum of `.eigvals`. Up to v0.6.0 `--svd 3` ran a
full SVD and wrote every eigenvalue. Since v0.7.0 it writes only the top `-k`,
so that sum no longer covers all the variance.

## Release history

### v0.8.0

**Breaking changes**, with the new commands in [Legacy command migration](#legacy-command-migration)

- `-D/--ld`, `-B/--binary` and the `.residuals` file are removed. LD (`-R`, `--ld-r2`, `--clump`)
  reads the PLINK or PGEN genotypes and removes the PCs of a previous run given with `-P/--USV`.
  BGEN, BEAGLE and CSV input no longer support LD.
- `-V/--printv` is removed: the `.loadings` and `.mbim` are written by default. Use
  `--no-loadings` to skip them.
- `--evaladmix` no longer runs a PCA; pass the PCs of a previous run with `-P/--USV`.
- Invalid option values and unsupported combinations stop with an error naming the option,
  instead of running on or crashing: e.g. `-k 0`, `-d 7`, more than one input file, `--clump-p1`
  above `--clump-p2`, EM-PCA with `--svd 3`, `--inbreed` with `--maf`, or `-k` not below the
  number of samples and sites. Errors print `Error: <message>` and exit with status 1.

**New**

- [Relatedness](small-n/relatedness.md): `--evaladmix` writes the kinship of every pair of
  samples (`.kinship`, `.corres`) from the residuals of the PCs, with no admixture run.
  `--evaladmix-ibd` adds the IBD sharing probabilities `k0`, `k2`, which tell parent–offspring
  from full sibs.
- [Relatedness at biobank scale](biobank/relatedness.md): `--evaladmix-kin <cutoff>` writes only
  the pairs above a kinship cutoff (`.kin0`) and a maximal unrelated set (`.unrelated`), within
  the memory set by `-m`; `--evaladmix-unrelated` gives the set its own cutoff.
- `--inbreed 2` estimates the inbreeding coefficient of each sample under population structure
  (`.inbred`). See [HWE and inbreeding](small-n/hwe.md).
- `--em-k` sets the number of PCs that model the allele frequencies in `--emu` and `--pcangsd`,
  separately from the `-k` PCs written.
- In every analysis that reads a reference with `-P/--USV` (LD, `--evaladmix`, `--project`,
  `--selection`, `--inbreed`), `-k` selects the leading PCs, and all of them are used by default,
  so one reference run serves any smaller number.
- BGEN layouts 1 and 2, every compression and bit depth, and phased files are tested, and BGEN
  input no longer warns that its support is limited.
- `.sigvals` records how the PCA scaled the data, so analyses with `-P/--USV` put new genotypes
  on the same scale. A `.sigvals` from an older version gives a warning.
- Out-of-core PCA reads the shuffled variants straight from BED, PGEN and BGEN, so it needs no
  extra disk space; `--bed-copy` writes a shuffled copy first, which is faster for a large BED on
  a spinning disk. The next block is read while the current one is used; `--no-prefetch` turns
  this off.
- Warnings when EM-PCA reaches `--maxiter` or the RSVD reaches `--maxp` without converging.

**Faster**

- BGEN is read by all threads: with 20 threads, 15x faster in-core and 7-8x with `-m`.
- `--svd 3` streams the data into the N x N GRM, in N x N memory instead of N x M (0.11 GB
  instead of 2.1 GB on 400 x 656,281), and is faster.
- The default `--svd 2` is about 2x faster, `--svd 1` 4x, `--emu` 3x, `--pcangsd` 3x and BEAGLE
  input up to 30x.
- Ancestry-adjusted LD is up to 14x faster for pruning and 38x for `-R` with `-m`, and in-core
  and `-m` now run equally fast.
- Out-of-core PGEN is up to 3.7x faster when the file is not in the page cache.

**Results change**

- The sign of each PC from `--svd 1` and `--svd 2` follows the rule of `--svd 0` and `--svd 3`:
  the largest entry of each PC is positive.
- `--seed` now fixes every random draw, on every platform and thread count, so the draws differ
  from earlier versions.
- BEAGLE input filters `--maf 0.05` by default, as PCAngsd does.
- `--inbreed 1` removed far too many sites as out of HWE; it now agrees with PCAngsd.
- Ancestry-adjusted LD with `-P` and more than one PC used the wrong subspace, and in-core
  `--clump` used the standard LD; both are fixed. The PCs of a default PCA are standardized, so
  add `--scale 0` to the PCA to reproduce the old `-D/--ld`.
- `--ld-r2` follows PLINK `--indep-pairwise` and prunes fewer variants; `--clump` no longer
  depends on the row order of the association file.
- `--project 3` put targets on the wrong scale and matched BEAGLE alleles the wrong way round.
  `--project-bootstrap` underestimated the uncertainty. `--project 2` writes NA for a sample with
  fewer called sites than PCs.
- `--selection` with `--maf` writes one row per input site, NA for the removed ones.
- Phased BGEN gave wrong dosages.
- In-core `--svd 0 --emu` computed the final PCs from the imputation of the previous iteration.
- CSV with `-m` is now transformed as in-core: not at all with the default `--scale` or with
  `--scale 0`.

**Fixes**

- Crashes: in-core `--svd 2` with fewer than 4,096 sites; `--clump` across chromosomes out-of-core;
  `--inbreed` with PGEN and `-m`; a missing `.mbim` with `-P`; several malformed CSV files;
  EM-PCA on very large matrices.
- `--selection 1` and `--project` with `-k` below the reference's number of PCs.
- LD `R2` of monomorphic variants is 0 instead of `NaN`.
- The BEAGLE `.eigvecs2` has one tab-separated column per PC.
- `--ncv` is used when given.

The full notes, with measurements, are in
[`dev/changes-v0.8.0.md`](https://github.com/Zilong-Li/PCAone/blob/main/dev/changes-v0.8.0.md).

### v0.7.2

- bug fix for out_of_range issue #23
- add `--project-bootstrap`
- optimization
### v0.7.1

- bug fix for sign flipping with `--project 3`
- bug fix for output mbim with PGEN input and out-of-core
### v0.7.0

- add support for [PLINK2 PGEN](https://www.cog-genomics.org/plink/2.0/input#pgen) input via `--pgen`
- bug fixes for projection with missing data
- use dosages from `.pgen` when available, with `--hardcall` to force hard-call genotypes
- add genotype-likelihood aware projection for [BEAGLE](http://www.popgen.dk/angsd/index.php/Input#Beagle_format) input via `--project 3`
- projection now matches overlapping markers against the reference `.mbim` file and corrects flipped alleles when needed
- add genome-wide selection scans via `--selection 1` (Galinsky/FastPCA) and `--selection 2` (pcadapt)
- `--svd 3` eigendecomposes the N x N sample GRM when N <= M instead of running a full SVD, which is much faster when N << M. `.eigvals` now holds only the top `-k` values
- always write `.mbim` together with `.loadings` to support downstream projection, HWE, and selection workflows
### v0.6.0

- **fix HWE! Rework pcangsd and emu plugins**
### v0.5.4

- **fix HWE! There is a bug in the previous release for HWE**
- add option `--seed`
### v0.5.3

- add support for linux aarch64 system. PR [#10](https://github.com/Zilong-Li/PCAone/pull/10) contributed by [@SauersML](https://github.com/SauersML)
### v0.5.2

- fix Makefile for bioconda
### v0.5.1

- bug fix: there is no standardization step for CSV
- feat: more options for scaling data via `--scale`
- feat: new `--scale-factor` option for adjusting scale
### v0.5.0

- optimization: reducing the binary size
- fix Makefile: use clang and Accelerate framework on MacOS (faster than any other solutions)
- fix a small bug in FileCSV
### v0.4.9

- use C+17 standard
- **breaking**: `--verbose` takes integer value as the level of verbose
- new short option: `-f`
### v0.4.8

- HWE test accounting for population structure via `--inbreed 1`
  - BEAGLE file can work with in-core mode
  - PLINK file can work with both in-core and out-of-core mode
- new option `--USV` as the prefix of pcaone `.eigvecs/.eigvals/.loadigns/.mbim` files
- simplied usage because both `--project` and `--inbreed` can work together with `--USV`
- the output `.eigvecs2` stores the eigenvectors of the covariance matrix, which also works for BEALGE file with genotype likelihood input.
- the `.eigvecs,.eigvals,.loadings` stores the U,S,V matrix for reconstructing the PI matrix
### v0.4.7

- **breaking**: rename option `--ld-bim` as `--match-bim`
- projection support via `--project`
  - 1: projection by multiplying the loadings with mean imputation for missing genotypes
  - 2: projection by solving the least squares system $Vx=g$. sites with missingness are skipped like smartPCA.
- new options: `--read-U, --read-V, --read-S`
- can do LD for subset plink file
### v0.4.6

- give warnings when sample standard deviation is zero.
- checks unknow options in CLI
- fix parsing bim file ending with empty line
- adjust memory estimator for LD
### v0.4.5

- **can do LD prunning and clumpint out-of-core**
- add `--ld-bims, --print-r2` options
### v0.4.4

- add `--clump-names` option
### v0.4.3

- new output `eigvecs2` for plink input
- fix ld clumping and sort snps by pvalues
### v0.4.2

- fix slowness due to the ld logics
- add ld-clump functions
### v0.4.1

- fix ld windows
- keep sites with higher maf in *prune.in*
- support standard pruning using `--svd 0` when PCAone is compiled with *DEBUG=1*
- add `--ld-snps` to output pairwise r2 for given SNPs
### v0.4.0

- add ld pruning
- fix `--maf` being af not maf
### v0.3.9

- pre-lease of ld feature
- verbose cli check
### v0.3.4

- fix the order of SNP loadings for plink input
### v0.3.3

- breaking change option `--svd` for choosing different SVD methods. default is PCAone algorithm2.
- add permuting `BGEN` file with multithreading.
- add Full SVD support by `--svd 3`
- add log transformation support for bulk RNA-seq data by `--scale 1`
- add binary file support
- recover the original order for SNP `loadings`.
- recover the eigenvalues for diploid genotype (0,1,2) data.
### v0.3.2

- use `tab` as separator for output files
### v0.3.1

- fix makefile for bioconda can't find zlib
### v0.3.0

- version after first manuscript revision corresponding to the third version manuscript on biorxiv
- fix algorithm2 to use one more omega updates between epoch.
- parameter `--band` changed to `--windows`
- support `CXXSTD` compiling option to use `c++17` standard. default `CXXSTD=c++11`
- add `PCAoneR`, which implement the idea of PCAone in R but without out-of-core
- colorful warning and error message. more checks and warnings.
### v0.2.1

- bug fix for PCAngsd
### v0.2.0

- add `--maf` option for SNPs filtering
- faster parser for beagle file
- optimization for PCAngsd
- fix printing loadings for PCAoneA blockwise mode
### v0.1.9

- command line options are re-designed
- upgrade to zstd v1.5.2
- recover the original order of SNPs loading for fancy batch mode
- add `--no-shuffle` option, remove `--shuffle` option
- default PCAoneF (fancy RSVD) algorithm is chose.
### v0.1.8

- add -a, –tmp options
- optimize makefile for conda
### v0.1.7

- two releases `x64` and `avx2` for linux
- add `conda install -c bioconda pcaone`
- upgrade Spectra to v1.0.1
- change oversamples default as max(10, k)
- add logger
### v0.1.6

- add –cpmed option to support raw counts for scRNAs
- output elapsed I/O time and total wall time
- remove pgenlib
### v0.1.5

- publish PCAone on mac with `libiomp5` support
- update documentation
- add –shuffle option
- optimization of FileCsv
### v0.1.4

- add CSV format support for scRNAs data
- add –printv to print out another eigen vectors
- add -N and -M options
### v0.1.3

- refactor Halko implementation
- use Arnoldi as default method
- add external/zlib
- upgrade to Eigen 3.4.0
- add structured permutation
### v0.1.2

- use double instead of float to improve numerical accuracy
- add –no-shuffle, –oversamples, –buffer options
- disable checking padding region of plink bed
- fix denominator too small
- set band range as 4,8,16,32,64
### v0.1.1

- supports bfile, bgen as input for both batch and blockwise mode
- implement super power iteration for Halko
- automatically permute plink bed file for fast halko
- port [jeremymcrae/bgen](https://github.com/jeremymcrae/bgen) as bgen parser.
- remove external BGEN v1.3 lib dependence
- use libiomp5 instead of libgomp as multithreading RTL
- fix bugs for PCAngsd algorithm
### v0.1.0

- supports bfile as input for both batch and blockwise mode
- supports bgen as input for batch mode
- supports both EMU and PCAngsd algorithm
- supports beagle as input for PCAngsd
- external dependecy: Eigen v3.3.8, Spectra v0.9.0, BGEN v1.3
- release two pre-compiled binary for Linux and Mac OSX (libiomp5 required).
