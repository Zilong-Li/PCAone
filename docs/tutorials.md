# Tutorials

Commands using `./PCAone` assume a local binary; use `PCAone` if it is on your PATH.

Let's download the example data first if you haven't done so.

```shell
wget http://popgen.dk/zilong/datahub/pca/example.tar.gz
tar -xzf example.tar.gz && rm -f example.tar.gz
```

## Genotype data (PLINK)

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

## PLINK2 PGEN input

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

## BGEN input (limited support)

**NB:** `BGEN support is very limited. Please convert BGEN to PGEN and use PGEN input!`

Imputation tools usually generate the genotype probabilities or dosages in
BGEN format. To do PCA with the imputed genotype probabilities, we can
work on BGEN file with `--bgen` option instead.

```shell
./PCAone --bgen example/test.bgen -k 10 -m 2
```

## Single cell RNA-seq data (CSV)

In this example, we run PCA for the scRNAs-seq data using CSV format with a
normalization method called count per median log transformation (CPMED).
Since the features (genes) tend to be not correlated locally, we use `-S`
option to disable permutation for winSVD.

```shell
./PCAone --csv example/BrainSpinalCord.csv.zst -k 10 -m 2 --scale 2 -S
```
