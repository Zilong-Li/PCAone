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

Out-of-core winSVD shuffles the SNPs. It assigns them at random (`--seed`)
to the `-w` bands of read blocks and keeps their source order within each
band. winSVD only updates its test matrix between bands, so the PCs equal
those of a full random permutation in exact arithmetic. Floating-point
rounding can differ slightly. Samples keep their order; `-S` turns the
shuffle off.

The shuffled SNPs are read from the input itself, so no disk space is
needed: BED, PGEN and BGEN read the variants of each block in file order.
On a spinning disk, a file that is not in the page cache costs a seek per
SNP, since each block takes about one SNP in `-w` from the whole file.
PCAone warns when the BED is on a spinning disk (Linux).

For a BED on a spinning disk that is larger than the free memory, use
`--bed-copy`. It writes the shuffled BED to `<out>.perm.*` before the PCA and
reads that copy in sequence: one read and one write of the BED, in large
writes whatever its size, and disk space for one more BED. On 20,000 samples
x 1.39M SNPs (7 GB) read from a spinning disk without the page cache, the PCA
took 224 s with `--bed-copy` and 1,840 s without; from an SSD, 209 s and
215 s; from the page cache, 168 s and 160 s. With 134,400 samples, a record
is 34 KB and a seek costs less: 136 s and 170 s on the spinning disk. Both
give the same bytes of every output.

With `--bed-copy`, `--buffer` (default 2 GiB) bounds the two genotype
buffers of the rewrite; it does not change the order. SNP indices and the
BIM lines take additional memory. The permuted BED/BIM/FAM (`<out>.perm.*`)
are removed after the PCA unless `-v 3` is used, and `-o` may not point them
at the input.

While the PCA works on one block, the next block is read in the background,
so reading overlaps the computation instead of alternating with it. BED is
read into a second buffer of one block of packed records (1/32 of the decoded
block), as long as reading a block takes more than 5% of the time the PCA
spends on one; a file in the page cache is read on the spot, without a second
thread competing for the cores, and the buffer is given back. For BGEN, PGEN
and the binary copy of CSV, the kernel is asked to read the next block into
the page cache (`POSIX_FADV_WILLNEED`, `F_RDADVISE` on macOS). That holds no
memory in PCAone, so the blocks take what `-m` gives them, as before. PGEN
and BGEN read the variants of a block in file order and ask for the next
block's records in file order (PGEN unless they are in the page cache
already), which matters most for the scattered reads of their logical
permutation on a cold disk. A file system that ignores the request is read as without it. The
results are the same; `--no-prefetch` reads in the foreground.

### Run sSVD method with out-of-core mode

```shell
./PCAone --bfile example/plink --svd 1 -m 2
```

### Run IRAM method with out-of-core mode

```shell
./PCAone --bfile example/plink --svd 0 -m 1
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
