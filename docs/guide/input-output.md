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

### Example data

The examples on this page and in the rest of the guide use the example data.
Download it, or run `make data` in the source tree:

```shell
wget http://popgen.dk/zilong/datahub/pca/example.tar.gz
tar -xzf example.tar.gz && rm -f example.tar.gz
```

Commands using `./PCAone` assume a local binary; use `PCAone` if it is on your PATH.

### PLINK

Run exact PCA with `--svd 3` and plot the sample coordinates in R.

```shell
./PCAone --bfile example/plink -d 3
```

Then, we can make a PCA plot in R.

```r
pcs <- read.table("pcaone.eigvecs2",h=F)
plot(pcs[,3:4], col=factor(pcs[,1]), xlab="PC1", ylab="PC2", cex.lab=1.5)
legend("topleft", legend=levels(factor(pcs[,1])), col=1:4, pch = 21, cex=1.2)
```

### PLINK2 PGEN

PCAone can also read PLINK2 `.pgen/.pvar/.psam` input directly. By default it
uses dosages when they are present in the `.pgen` file.

```shell
./PCAone --pgen example/plink2 -k 10 -m 2
```

If you want to ignore dosages and use hard-call genotypes instead, add
`--hardcall`.

```shell
./PCAone --pgen example/plink2 -k 10 -m 2 --hardcall
```

### BGEN (limited support)

**NB:** `BGEN support is very limited. Please convert BGEN to PGEN and use PGEN input!`

Imputation tools usually generate the genotype probabilities or dosages in
BGEN format. To do PCA with the imputed genotype probabilities, we can
work on BGEN file with `--bgen` option instead.

```shell
./PCAone --bgen example/test.bgen -k 10 -m 2
```

### CSV (single-cell RNA-seq)

In this example, we run PCA for the scRNAs-seq data using CSV format with a
normalization method called count per median log transformation (CPMED).
Since the features (genes) tend to be not correlated locally, we use `-S`
option to disable permutation for winSVD.

```shell
./PCAone --csv example/BrainSpinalCord.csv.zst -k 10 -m 2 --scale 2 -S
```

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
[the plotting guide](plotting.md#snp-loadings).

### Variant information

A plink-like bim file named with `.mbim` is used to store the variants list
with extra information. Currently, the `mbim` file has 7 columns with the 7th
being the allele frequency. PCAone writes this file automatically whenever
it outputs `.loadings` via `--printv`.

### LD R2

The LD R2 for pairwise SNPs within a window can be outputted to a file
with suffix `ld.gz` via `--print-r2` option. This file uses the same long
format as the one [plink](https://www.cog-genomics.org/plink/1.9/ld#r) used.
