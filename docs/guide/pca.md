# Choosing a PCA method

PCAone has four ways to compute the PCs, selected with `-d/--svd`, and runs
each of them either with the data in memory (in-core, the default) or block by
block from the file (out-of-core, `-m`). Which one to use depends mostly on the
number of samples (`N`), the number of variants or features (`M`) and the
number of PCs (`k`).

| Method        | Scenario             | Accuracy  | Speed                            |
|---------------|----------------------|-----------|----------------------------------|
| Exact (-d 3)  | small `N`, `N << M`  | Exact     | fastest for `N` below ~2000      |
| winSVD (-d 2) | `M or N >> 500000`   | Very high | fast (only 7 iterations used)    |
| IRAM (-d 0)   | speed insensitive    | Very high | depends on `N` and # iterations |
| sSVD (-d 1)   | accuracy insensitive | High      | depends on # iterations         |

In short:

- **Up to a few thousand samples**, use the exact PCA, `--svd 3`. It needs the
  `N x N` sample GRM in memory and nothing else. See
  [Exact PCA and EM-PCA](../small-n/pca.md).
- **Tens of thousands of samples or more**, use the default winSVD, `--svd 2`,
  and add `-m` when the genotypes do not fit in memory. See
  [Out-of-core PCA](../biobank/pca.md).
- **Many missing calls, or genotype likelihoods** (low-depth sequencing), use
  the EM-PCA of `--emu` or `--pcangsd` (implied by BEAGLE input). See
  [EM-PCA](../small-n/pca.md#em-pca-emu-and-pcangsd).

## In-core and out-of-core

IRAM (`--svd 0`), sSVD (`--svd 1`) and winSVD (`--svd 2`) each run in-core or
out-of-core, and the exact PCA (`--svd 3`) streams the data into the GRM, so
there are 7 ways to compute the PCs.

In-core is the default and the fastest, as long as the data fit in memory. You
can gain some speed in-core by limiting `-n` to half of the available threads.
If the data do not fit, give `-m/--memory` a value in GB greater than 0 to read
them block by block. `-m 1` or `-m 2` is usually enough for large datasets;
PCAone allocates more when it needs to, so the total RAM can exceed `-m`.

```shell
# winSVD (default), in-core
./PCAone --bfile example/plink

# winSVD, out-of-core
./PCAone --bfile example/plink -m 2

# sSVD, out-of-core
./PCAone --bfile example/plink --svd 1 -m 2

# IRAM, out-of-core
./PCAone --bfile example/plink --svd 0 -m 1

# exact PCA; streams the GRM when N <= M
./PCAone --bfile example/plink --svd 3
```

winSVD shuffles the variants, in-core and out-of-core, so that each of its
mini-batches is a random sample of the genome rather than a run of variants in
LD. Use `-S/--no-shuffle` for data without local correlation, such as gene
counts. `--seed` fixes the shuffle and every other random draw, on every
platform and thread count.

## Data normalization

Use `-C/--scale` to choose the normalization. Genetic data are centred and
standardized by allele frequency by default, `sqrt(ploidy * f * (1 - f))`. CSV
input is left as it is (neither centred nor scaled) with the default and
`--scale 0`; modes 1–4 standardize, after the optional count or log
transformation described under [`--scale`](options.md).

The scaling is recorded in the `.sigvals`, so analyses that read the PCs back
with `-P/--USV` (projection, selection scans, HWE) put new genotypes on the
same scale without being told.
