# Ancestry-adjusted LD

LD patterns vary across diverse ancestry and structured groups, and
conventional LD statistics, e.g. the implementation in `plink --ld`, failed to
model the LD in admixed populations. Thus, we can use the so-called
ancestry-adjusted LD statistics to account for population structure in
LD. See our [paper](https://doi.org/10.1101/2024.05.02.592187) for more details.

The ancestry-adjusted LD is the correlation between the residuals of the
genotypes after removing the top principal components, which capture the
population structure. We first figure out the number of PCs (`-k/--pc`) that
capture population structure and run the PCA. In this example, assuming that 3
PCs can account for population structure:

```shell
./PCAone -b example/plink -k 3 -o adj
```

Then pass its prefix to `-P/--USV` in any of the LD analyses below. PCAone
reads the genotypes again, removes the PCs in `adj.eigvecs` from them as they
are read, and computes the LD from the residuals, so no residual matrix is
written to disk. The PCs must come from the same samples in the same order, but
the variants need not be the ones the PCA used (e.g. a different `--maf`). This
works in-core and out-of-core (`-m`), for PLINK (`-b`) and PGEN (`-p`) input.
By default every PC in `.eigvecs` is removed; `-k/--pc` removes only the
leading ones, so one reference PCA run with e.g. `-k 10` serves any smaller
number of PCs (`-P adj -k 3`). Use `--ld-stats 1` without `-P` for the standard
LD instead. See
[the LD methods guide](../methods/ld-ancestry-adjusted.md) for the derivation and implementation.

## Report LD statistics

Currently, the LD R2 for pairwise SNPs within a window can be outputted via `--print-r2` option.

```shell
./PCAone -b example/plink \
         -P adj \
         --ld-bp 1000000 \
         --print-r2 \
         -o adj
```

To plot LD decay curves from the `.ld.gz` files, e.g. standard (`--ld-stats 1`)
against ancestry-adjusted LD, use [plot-ld-decay.R](https://github.com/Zilong-Li/PCAone/blob/main/scripts/plot-ld-decay.R):

```shell
Rscript scripts/plot-ld-decay.R adj.ld.gz std.ld.gz --labels Adjusted,Standard -n example/plink.fam
```

The nextflow workflow [ld.nf](https://github.com/Zilong-Li/PCAone/blob/main/workflows/ld.nf) runs the whole comparison. See
[the plotting guide](plotting.md#ld-decay) for both.

## Pruning

Given the PCs, we can do pruning based on user-defined thresholds and windows.
Of each pair of variants above the R2 cutoff, the one with the lower MAF is
removed. The kept and removed variants are written to `.ld.prune.in` and
`.ld.prune.out`.

```shell
./PCAone -b example/plink \
         -P adj \
         --ld-r2 0.8 \
         --ld-bp 1000000 \
         -o adj
```

## Clumping

Likewise, we can do clumping based on the ancestry-adjusted LD and
user-defined association results

```shell
./PCAone -b example/plink \
         -P adj \
         --clump example/plink.pheno0.assoc,example/plink.pheno1.assoc  \
         --clump-p1 0.01 \
         --clump-p2 0.05 \
         --clump-r2 0.1 \
         --clump-bp 10000000 \
         -o adj
```
