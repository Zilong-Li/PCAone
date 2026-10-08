# PCA methods and memory

## Which SVD method to use

This depends on your datasets, particularly the relationship between number
of samples (`N`) and the number of variants / features (`M`) and the top PCs
(`k`). Here is an overview and the recommendation.

| Method        | Scenario             | Accuracy  | Speed                            |
|---------------|----------------------|-----------|----------------------------------|
| Exact (-d 3)  | small `N`, `N << M`  | Exact     | fastest for `N` below ~2000      |
| winSVD (-d 2) | `M or N >> 500000`   | Very high | fast (only 7 iterations used)    |
| IRAM (-d 0)   | speed insensitive    | Very high | depends on `N` and # iterations |
| sSVD (-d 1)   | accuracy insensitive | High      | depends on # iterations         |

## Performance and memory

PCAone has both **in-core** and **out-of-core** mode for 3 different partial SVD
algorithms, which are IRAM (`--svd 0`), sSVD (`--svd 1`) and winSVD (`--svd 2`).
Also, PCAone supports an exact PCA (`--svd 3`). When N <= M, it streams the
data block by block and holds only the `N x N` GRM (`-m` sets the block size);
with `--maf`, or when `N > M`, it runs in-core.
Therefore, there are **7** ways for doing PCA in PCAone. In default PCAone uses
**in-core** mode, which is the fastest way (**NOTE**: you can gain some speedup for
in-code computation by limiting the `-n threads` to half of the available
threads of your machine). However, in case the server runs out of memory,
you can trigger `out-of-core mode` by specifying the amount of memory using
`-m/--memory` option with a value greater than 0. Normally, use `-m 1` is enough
for large dataset and PCAone will allocate more RAM when needed.

### Run winSVD method (default) with in-core mode

```shell
./PCAone --bfile example/plink
```

### Run winSVD method with out-of-core mode

```shell
./PCAone --bfile example/plink -m 2
```

Out-of-core BED winSVD writes a shuffled copy of the BED before the PCA. It
assigns the SNPs at random (`--seed`) to the `-w` bands of read blocks and
keeps their source order within each band. winSVD only updates its test
matrix between bands, so the PCs equal those of a full random permutation in
exact arithmetic. Floating-point rounding can differ slightly. The rewrite
reads the BED once and writes it once, in large writes whatever its size.
It needs disk space for one more BED. Samples keep their order; `-S` turns
the shuffle off.

`--buffer` (default 2 GiB) bounds the two genotype buffers of the rewrite;
it does not change the order. SNP indices and the BIM lines take additional
memory. The permuted BED/BIM/FAM (`<out>.perm.*`) are removed after the PCA
unless `-v 3` is used, and `-o` may not point them at the input.

### Run sSVD method with out-of-core mode

```shell
./PCAone --bfile example/plink --svd 1 -m 2
```

### Run IRAM method with out-of-core mode

```shell
./PCAone --bfile example/plink --svd 0 -m 2
```

### Run the exact PCA (streams the GRM when N <= M)

```shell
./PCAone --bfile example/plink --svd 3
```

## Data normalization

Use `--scale` to choose normalization. Genetic data are centred and
standardized by allele frequency by default. CSV input is left as it is
(neither centred nor scaled) with the default and `--scale 0`; modes 1–4
standardize, after the optional count or log transformation described under
[`--scale`](options.md).

## EM-PCA: EMU and PCAngsd

`--emu` runs the EM-PCA of [EMU](https://github.com/Rosemeis/emu) for
genotypes with many missing calls, and `--pcangsd` (implied by `--beagle`)
that of [PCAngsd](https://github.com/Rosemeis/pcangsd) for genotype
likelihoods. Each iteration estimates the individual allele frequencies from
the top PCs, rebuilds the matrix from them (imputing the missing calls, or
taking the expected genotypes), then recomputes the PCs, until they converge
(`--tol-em`, `--maxiter`).

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

The EM iterations cost what a run with `-k` equal to `--em-k` costs; only one
final decomposition computes `-k` PCs, of the matrix rebuilt from the
`--em-k` PCs: standardized for EMU, the expected genotypes for PCAngsd.
`.eigvecs`, `.eigvals`, `.sigvals` and `.loadings` hold the `-k` PCs.
For BEAGLE input, the `.cov` is PCAngsd's covariance under the `--em-k` model,
the matrix PCAngsd writes with `--eig` set to it, and `.eigvecs2` holds its top
`-k` eigenvectors. `--em-k` works with `--svd 0`, `1` and `2`, in-core and
out-of-core (`-m`) where the input allows it.
