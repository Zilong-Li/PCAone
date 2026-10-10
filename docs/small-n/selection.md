# Selection scan

A selection scan looks for variants that are more differentiated along the PCs
than the rest of the genome, a signature of local adaptation. PCAone
computes two such statistics from a reference PCA.

First compute the PCs from a set of good sites, after site-level QC such as
LD pruning, so that the PCs reflect genome-wide structure rather than a few
regions:

```shell
PCAone -b example/subset -k 10 -o goodsites
```

Then scan every variant with `--selection`, passing the reference prefix with
`-P/--USV`. Method `1` computes the Galinsky/FastPCA statistic for each variant
and PC; method `2` computes pcadapt-style z-scores, a chi-square statistic per
variant, p-values and a genomic inflation factor.

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

The scan reads the genotypes of the same samples, in the same order, as the
reference PCA; the variants can differ. It uses every PC of the reference, or
the leading `-k`. It runs out-of-core with `-m` too.

## Output

Each output has one row per variant, in the order of the `.bim`/`.pvar`, and no
variant IDs. Variants removed by `--maf`, and variants with no variance
(monomorphic, or fully explained by the PCs) are written as NA, so the rows
always line up with the input.

| `--selection` | file | contents |
|---|---|---|
| `1` | `.galinsky` | statistic per variant and PC, chi-square with 1 df |
| `1` | `.galinsky.pval` | its p-value per variant and PC |
| `2` | `.zscore` | z-score per variant and PC |
| `2` | `.pcadapt` | squared robust Mahalanobis distance per variant |
| `2` | `.pcadapt.chi2` | the same divided by the genomic inflation factor |
| `2` | `.pcadapt.pval` | its p-value, chi-square with `k` df |
| `2` | `.pcadapt.gif` | the genomic inflation factor |

Use method `1` to ask which PC a signal is on, and method `2` for one test per
variant across all PCs.

## Method details

Both statistics start from the regression of each variant on the PCs. With `g`
the standardized genotypes of a variant (scaled as the reference PCA was, as
recorded in its `.sigvals`) and `U` the `N x k` orthonormal PCs, the
coefficients are `b = U' g`.

**Galinsky (FastPCA)** (Galinsky et al. 2016, *AJHG* 98:456). Dividing `b` by
the singular values gives the variant's loading `v_jk` on each PC, i.e. its
entry in the right singular vector. Under neutrality `M v_jk^2` follows a
chi-square distribution with 1 df, where `M` is the number of variants.

**pcadapt** (Luu et al. 2017, *Mol Ecol Resour* 17:67). The z-score of
variant `j` on PC `k` is `b_jk / sigma_j`, where `sigma_j^2 = (|g_j|^2 - |b_j|^2) / (N - k)`
is the residual variance of the regression. The statistic is the squared
Mahalanobis distance of the vector of z-scores, with a robust estimate of their
mean and covariance: the orthogonalized Gnanadesikan-Kettenring estimator
(Maronna & Zamar 2002), with the defaults of `bigutilsr::dist_ogk()` that
pcadapt uses. The genomic inflation factor is its median divided by the median
of a chi-square with `k` df, and the p-values come from that chi-square after
dividing by it.
