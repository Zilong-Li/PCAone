# Genome-wide selection scan

PCAone can test selection for differentiated variants after a good reference
PCA has been computed with good sites after site-level QC, e.g. LD pruning.

First run PCA on the subset genotype matrix with good sites

```shell
PCAone -b example/subset -k 10 -o goodsites
```

Then reuse the saved `USV` outputs with `--selection`. Method `1` computes the
Galinsky/FastPCA statistic for each variant and PC, while method `2` computes
pcadapt-style z-scores, chi-square statistics, p-values, and a genomic
inflation factor.

```shell
PCAone -b example/fullset \
       --USV goodsites \
       --selection 1 \
       -o sel
```

```shell
PCAone -b example/fullset \
       --USV goodsites \
       --selection 2 \
       -o sel
```

The selection outputs are written per variant in the same marker order as the
plink `.bim` file.
