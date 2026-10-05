# Input and output

## Input formats

PCAone is designed to be extensible to accept many different formats.
Currently, PCAone can work with SNP major genetic formats to study
population structure, such as [PLINK](https://www.cog-genomics.org/plink/1.9/formats#bed), [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format) and [Beagle](http://www.popgen.dk/angsd/index.php/Input#Beagle_format). Also, PCAone supports
a comma delimited CSV format compressed by zstd, which is useful for other
datasets requiring specific normalization such as single cell RNAs data. PCAone also supports [PLINK2 PGEN](https://www.cog-genomics.org/plink/2.0/input#pgen)
input via `--pgen`. If dosages are stored in the `.pgen` file, PCAone uses
them by default; add `--hardcall` to force hard-call genotypes instead. The
current `BGEN` support is limited, so for large production workflows we
recommend converting `BGEN` to `PGEN` when possible.

## Output files

### Eigen vectors

Eigen vectors are saved in file with suffix `.eigvecs`. Each row represents
a sample and each col represents a PC.

### Eigen/Singular values

Eigenvalues and singular values are saved in file with suffix `.eigvals` and
`.sigvals` respectively. Each row represents the eigenvalue/singularvalue of
corresponding PC.

### Features loadings

Features Loadings are saved in file with suffix `.loadings`. Each row
represents a feature and each column represents a corresponding PC. Use
`--printv` option to output it.
To plot them along the genome, see
[the plotting guide](../plotting.md#snp-loadings).

### Variant information

A plink-like bim file named with `.mbim` is used to store the variants list
with extra information. Currently, the `mbim` file has 7 columns with the 7th
being the allele frequency. PCAone writes this file automatically whenever
it outputs `.loadings` via `--printv`.

### LD R2

The LD R2 for pairwise SNPs within a window can be outputted to a file
with suffix `ld.gz` via `--print-r2` option. This file uses the same long
format as the one [plink](https://www.cog-genomics.org/plink/1.9/ld#r) used.
