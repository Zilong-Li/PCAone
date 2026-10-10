# Exact PCA and EM-PCA

## Exact PCA

For up to a few thousand samples, the exact PCA is the fastest method and has
no approximation error. It computes the top `-k` eigenvectors of the `N x N`
sample GRM:

```shell
./PCAone --bfile example/plink --svd 3 -k 10
```

Plot the first two PCs in R, coloured by the population in the FID column:

```r
pcs <- read.table("pcaone.eigvecs2", h = F)
plot(pcs[, 3:4], col = factor(pcs[, 1]), xlab = "PC1", ylab = "PC2", cex.lab = 1.5)
legend("topleft", legend = levels(factor(pcs[, 1])), col = 1:4, pch = 21, cex = 1.2)
```

What to know:

- Memory is the GRM, `8 N^2` bytes, plus one block of genotypes. The genotype
  matrix is never held, so 400 samples x 656,281 SNPs take 0.11 GB. `-m` sets
  the size of the blocks.
- The `.loadings` take a second pass over the data. Add `--no-loadings` to
  skip it when you only need the sample PCs.
- With `--maf`, or with more samples than sites, the exact PCA runs in-core
  and ignores `-m`, with a warning.
- It does not support EM-PCA (`--emu`, `--pcangsd`); use `--svd 0`, `1` or `2`
  for those.
- Only the top `-k` eigenvalues are written, so the sum of `.eigvals` is not
  the total variance. See the [change log](../changelog.md#variance-explained-with---svd-3-before-v070).

## EM-PCA: EMU and PCAngsd

`--emu` runs the EM-PCA of [EMU](https://github.com/Rosemeis/emu) for
genotypes with many missing calls, and `--pcangsd` (implied by `--beagle`)
that of [PCAngsd](https://github.com/Rosemeis/pcangsd) for genotype
likelihoods from low-depth sequencing.

```shell
# PLINK genotypes with missing calls
./PCAone --bfile example/plink --emu -k 3
# genotype likelihoods from ANGSD
./PCAone --beagle example/beagle.gz -k 3
```

BEAGLE input filters sites with MAF below 0.05 by default, as PCAngsd does;
`--maf 0` turns the filter off. Its output includes the `.cov`, PCAngsd's
covariance matrix of the samples.

By default the same `-k` PCs model the allele frequencies and are written.
`--em-k` sets the number of PCs that model them, and `-k` the number written
from the final matrix, as EMU's `--eig` and `--eig-out` do. This lets a few PCs
that capture the population structure fit the model, while more PCs are
reported:

```shell
# model the allele frequencies with 2 PCs, write 10 PCs
./PCAone --bfile example/plink --emu --em-k 2 -k 10
./PCAone --beagle example/beagle.gz --em-k 2 -k 10
```

The EM iterations stop when the PCs converge (`--tol-em`) or after `--maxiter`
iterations, with a warning if they did not converge. `--emu` reads PLINK, PGEN
and BGEN input, and runs out-of-core with `-m` except for BGEN. BEAGLE input
runs in-core only.

## Method details

### Exact PCA

PCAone standardizes each block of sites as it is read, by
`sqrt(ploidy * f * (1 - f))` by default, and adds the block to the lower
triangle of the GRM with a multithreaded rank update. With N <= M this streams
the data, in-core and with `-m`, so the memory is the GRM plus one block (about
64 MB, or as much as `-m` allows). The top `-k` eigenpairs of the GRM come
from Lanczos iterations (Spectra, tolerance 1e-12) for N > 1000, and from the
full eigendecomposition below that. The sign of each PC is set so that the
entry of largest magnitude is positive. The loadings, `V = G' U S^-1`, take a
second pass over the data.

Times with 2 threads:

| data | `--svd 3` | other |
|---|---|---|
| `example/plink` (400 x 656,281), `-k 10` | 2.1 s | PLINK 2 `--pca` 2.9 s, default `--svd 2` 8.5 s |
| 2,000 x 100,000, with loadings | 6.5 s | |
| 4,000 x 20,000, with loadings | 6.7 s | |

### EM-PCA

Each iteration estimates the individual allele frequencies from the top PCs,
rebuilds the matrix from them (imputing the missing calls for EMU, taking the
expected genotypes given the likelihoods for PCAngsd), then recomputes the PCs,
until they converge. Each RSVD in the iterations starts from the previous PCs.

With `--em-k`, the iterations cost what a run with `-k` equal to `--em-k`
costs. Only one final decomposition computes `-k` PCs, of the matrix rebuilt
from the `--em-k` PCs: standardized for EMU, the expected genotypes for
PCAngsd. `.eigvecs`, `.eigvals`, `.sigvals` and `.loadings` hold the `-k` PCs.
For BEAGLE input, the `.cov` is PCAngsd's covariance under the `--em-k` model,
the matrix PCAngsd writes with `--eig` set to it, so it does not depend on
`-k`, and `.eigvecs2` holds its top `-k` eigenvectors. `--em-k` works with
`--svd 0`, `1` and `2`, in-core and out-of-core (`-m`) where the input allows
it. Without `--em-k`, or with it equal to `-k`, the output is that of the plain
EM-PCA.

Checked against an exact EMU and PCAngsd in numpy, with a full SVD in each
iteration, at a matched `--tol-em`, for `--em-k` below and above `-k`: the PCs
of IRAM agree to 0.05 degrees and those of the RSVDs to 0.45 degrees, and the
`.cov` to 1e-5.

The PCAngsd covariance divides each site by `2f(1-f)`, so rare sites dominate
it, which is why BEAGLE input filters `--maf 0.05` by default. With `--maf 0`,
a site whose EM frequency is 0 is left out of the diagonal of the `.cov`.
