# Out-of-core PCA

For a cohort whose genotypes do not fit in memory, run the default winSVD
(`--svd 2`) out-of-core by giving `-m/--memory` the memory in GB for the
blocks of data:

```shell
# PLINK BED
PCAone -b cohort -k 20 -m 16 -o pcs
# PLINK2 PGEN, using dosages when the file has them
PCAone -p cohort -k 20 -m 16 -o pcs
# BGEN, e.g. imputed genotypes
PCAone -g cohort.bgen -k 20 -m 16 -o pcs
```

`-m` sizes the working blocks, so the total RAM can exceed it; `-m 1` or
`-m 2` is usually enough even for large data. winSVD converges in about 7
passes (epochs) over the data. IRAM (`--svd 0`) and sSVD (`--svd 1`) also run
out-of-core; see [Choosing a PCA method](../guide/pca.md).

What to know:

- **Threads.** `-n` sets the threads. Build PCAone with MKL or OpenBLAS (see
  [Install](../installation.md)) for the fastest products on many cores.
- **Variant filters.** `--maf` is not available out-of-core. Filter rare
  variants, and LD-prune if you like, beforehand, e.g. with
  `plink2 --maf 0.01 --make-bed`, or with
  [ancestry-adjusted LD pruning](ld.md).
- **Loadings.** The `.loadings` and `.mbim` are written by default and are
  needed by [projection](projection.md). Add `--no-loadings` if you only want
  the sample PCs.
- **PGEN dosages.** PCAone uses the dosages in a `.pgen` when it has them;
  add `--hardcall` for the hard calls.
- **BGEN.** Layouts 1 and 2 (v1.1 to v1.3), uncompressed or compressed by zlib
  or zstd, any bit depth, phased or unphased. PCAone uses the dosage of the
  minor allele of each biallelic variant, and a sample with the missing flag
  counts as missing. `--emu` with BGEN runs in-core only.
- **Reproducibility.** The shuffle and every other random draw follow
  `--seed`, on every platform and thread count.

## The shuffle and the disk

winSVD works through the variants in mini-batches, so it shuffles them so that
each mini-batch is a random sample of the genome rather than a run of
variants in LD. Out-of-core, it assigns the variants at random (`--seed`) to
the `-w` bands of read blocks and keeps their source order within each band.
`-S/--no-shuffle` turns the shuffle off, for data without local correlation;
the samples always keep their order.

The shuffled variants are read from the input itself, so no extra disk space is
needed: BED, PGEN and BGEN read the variants of each block in file order. On a
spinning disk, a file that is not in the page cache costs a seek per variant,
since each block takes about one variant in `-w` from the whole file. PCAone
warns when the BED is on a spinning disk (Linux).

For a BED on a spinning disk that is larger than the free memory, use
`--bed-copy`. It writes the shuffled BED to `<out>.perm.*` before the PCA and
reads that copy in sequence: one read and one write of the BED, in large
writes whatever its size, and disk space for one more BED.

| 20,000 x 1.39M SNPs (7 GB) | without `--bed-copy` | with `--bed-copy` |
|---|---|---|
| spinning disk, not in the page cache | 1,840 s | 224 s |
| SSD | 215 s | 209 s |
| page cache | 160 s | 168 s |

With 134,400 samples a record is 34 KB and a seek costs less: 170 s without
and 136 s with `--bed-copy` on the spinning disk. Both give the same bytes of
every output.

With `--bed-copy`, `--buffer` (default 2 GiB) bounds the two genotype buffers
of the rewrite; it does not change the order. SNP indices and the BIM lines
take additional memory. The permuted BED/BIM/FAM (`<out>.perm.*`) are removed
after the PCA unless `-v 3` is used, and `-o` may not point them at the input.

## Method details

### winSVD

winSVD is the window-based randomized SVD of
[Li et al. (2023)](https://genome.cshlp.org/content/early/2023/10/05/gr.277525.122).
It splits the variants into `-w` mini-batches (64 by default) and updates the
randomized range estimate after a band of them, not only after a full pass. The
band doubles every epoch (2, 4, 8, ... mini-batches) until it covers all the
variants, so the estimate improves many times within the first passes over the
data, and it converges in about 7 epochs (`log2(64) + 1`). The sign of each PC
is set so that the entry of largest magnitude is positive.

Because winSVD only updates the range estimate between bands, assigning the
variants at random to the bands and keeping their order within a band gives the
same PCs as a full random permutation, in exact arithmetic. Floating-point
rounding can differ slightly. Shuffling runs of neighbouring variants instead
would keep the reads large, but winSVD then converged more slowly (9 to 21
epochs instead of 7, with PCs further from the exact ones).

### Reading ahead

While the PCA works on one block, the next block is read in the background,
so reading overlaps the computation instead of alternating with it.

- **BED** is read into a second buffer of one block of packed records (1/32 of
  the decoded block), as long as reading a block takes more than 5% of the time
  the PCA spends on one. A file in the page cache is read on the spot, without
  a second thread competing for the cores, and the buffer is given back.
- **BGEN, PGEN** and the binary copy of CSV: the kernel is asked to read the
  next block into the page cache (`POSIX_FADV_WILLNEED`, `F_RDADVISE` on
  macOS). That holds no memory in PCAone, so the blocks take what `-m` gives
  them. PGEN and BGEN ask for the next block's records in file order (PGEN
  unless they are in the page cache already), which matters most for the
  scattered reads of the shuffle on a cold disk.

A file system that ignores the request is read as without it. The results are
the same; `--no-prefetch` reads in the foreground. BED 10,000 x 200,000, `-S
-n 2`, on an emulated 200 / 100 MB/s disk: 53.8 / 67.6 s, against 89.6 / 122.6 s
with `--no-prefetch`.

### Decoding

- **PGEN.** Each thread reads a run of nearby records of the block, in file
  order, through pgenlib. Centred hard calls go from the 2-bit calls straight
  into the block, dosages over them. 2,000 x 500,000 with dosages (1.7 GB),
  with the page cache capped below the file: 76 s; in the page cache: 45 s
  (hard calls 39 s).
- **BGEN.** The variants are indexed by one pass over their headers when the
  file is opened (0.4 s for 2.7 GB). Each thread then reads, inflates and
  decodes variants of its own, in-core and with `-m`. With 20 threads, on
  simulated data, 10,000 x 200,000 with `-m 2` takes 59 s, and
  20,000 x 500,000 with `-m 4` 344 s. The dosages divide by the largest
  probability value, so certain genotypes are exactly 0, 1 and 2 at any bit
  depth. For a phased file the dosage is the sum over the haplotypes, so it
  gives the PCs of the same genotypes unphased.

### Threading the products

The product of each block with the test matrix is threaded over panels of rows
of the block, on all of `-n`. With 192 samples or more, the PCs are those of
one thread whatever `-n` is. BED 10,000 x 200,000, `-n 2`, in the page cache:
52-54 s. Builds with MKL, OpenBLAS or Accelerate use the BLAS product instead.
