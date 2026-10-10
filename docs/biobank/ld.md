# LD pruning and clumping at biobank scale

The [ancestry-adjusted LD](../small-n/ld.md) analyses run out-of-core with
`-m`, from the PCs of an out-of-core PCA. They read the genotypes again and
remove the PCs as each block is read, so the residuals are never written to
disk and the memory stays within `-m`:

```shell
# the PCs
PCAone -b cohort -k 20 -m 16 -o pcs
# prune on the LD left after the leading 10 PCs
PCAone -b cohort -P pcs -k 10 --ld-r2 0.1 --ld-bp 1000000 -m 16 -o pruned
# clump association results
PCAone -b cohort -P pcs -k 10 --clump gwas.assoc --clump-p1 5e-8 -m 16 -o clumped
```

Pruning writes `pruned.ld.prune.in` and `pruned.ld.prune.out`, for
`plink2 --extract`; clumping writes `clumped.p0.clump`. The options, the
outputs and the choice of PCs are described in
[Ancestry-adjusted LD](../small-n/ld.md). PLINK (`-b`) and PGEN (`-p`) input
are supported, and in-core and `-m` give identical output.

What to know:

- **The PCs** must come from the same samples, in the same order, but not from
  the same variants. Compute them once, out-of-core, and select the leading ones
  with `-k`.
- **Windows and `-m`.** A window of `--ld-bp` must fit in two blocks of `-m`;
  otherwise PCAone asks for a larger `-m` or a smaller `--ld-bp`.
- **`--maf`** applies in-core only. Out-of-core, filter rare variants
  beforehand, e.g. with `plink2 --maf 0.01 --make-bed`.
- **`-R/--print-r2`** writes every pair within `--ld-bp`, and the number of
  pairs grows with the density of the variants times the window: 104k SNPs in
  1 Mb windows give 25M pairs and a 140 MB file. Thin dense data first, e.g.
  with `plink --thin`, and bin large files with
  [`summarise_ld_r2bin`](../guide/plotting.md#large-files-summarise_ld_r2bin)
  to plot the LD decay.

## Method details

The residuals of a site are `(I - Q Q') g`, with `Q` from a QR decomposition of
the `.eigvecs`; see [the derivation](../small-n/ld.md#method-details).

`LDColumns` (`src/LD.cpp`) holds the residuals of the sites in memory: the whole
matrix in-core, or with `-m` the two consecutive blocks the current windows
need. Each column is scaled to unit norm once, so the correlation of two sites
is the dot product of their columns. A site with no variance (monomorphic, or
fully explained by the PCs) becomes a zero column, and its R2 is reported as 0.

The correlations of a batch of windows come from one matrix product,
`C = X_rows' X_cols`, rather than one dot product per pair. That reuses every
column while it is in cache.

- **Pruning** takes the next windows whose lead site is still kept. It
  multiplies them against the kept sites of their span, then applies the greedy
  rule (remove the lower-MAF site of each pair above `--ld-r2`) in window order.
  This makes the same decisions as one window at a time. The batch grows while
  leads survive and shrinks while they are pruned away.
- **`-R`** multiplies consecutive windows over their span. It then formats and
  gzip-compresses the lines on all threads, one gzip member per thread. The
  `.ld.gz` is a series of members, which `zcat`, `gzip -d`, R and Python read as
  one file.
- **Clumping** reads its target sites and takes dot products.

In-core and `-m` go through the same code. On 4 threads and simulated data
with N x M = 1e8 (N = 2,000 to 20,000, 1 Mb windows, up to 50M pairs), the
genotypes take 16x less to read than the `.residuals` file of PCAone before
v0.8.0 (2 bits against 4 bytes per genotype), and the batched products make it
faster still: pruning is 5-14x faster with `-m` than before, and `-R` up to
38x, so `-m` now runs as fast as in-core. The output is unchanged.
