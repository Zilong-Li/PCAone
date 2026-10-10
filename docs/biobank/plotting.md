# Plotting large cohorts

A scatter plot of hundreds of thousands of samples is a solid blob, slow to
draw and to open. The scripts below summarise instead of drawing every point,
so they stay fast and readable at biobank scale. They are in `scripts/` and
need R >= 4.0 and the `data.table` package.

## Sample PCA density

`scripts/plot-pca.R` plots large cohorts, including 500,000 samples, as a
2D density raster. Every sample contributes to a bin; no subsampling is used.
Colours show log-scaled sample counts, with empty bins white. This shows
population density rather than individual points or population labels.

```bash
Rscript -e 'install.packages("data.table")'
Rscript scripts/plot-pca.R pcaone.eigvecs2 -o pca.png
Rscript scripts/plot-pca.R pcaone.eigvecs --pcs 3,4 -o pca.pdf
```

The script reads only the two selected PCs. It recognises headers containing
`PC1`, `PC2`, etc., and PCAone's `.eigvecs2` format with FID/IID columns.
For other headerless files with two ID columns, specify `--format ids`;
for numeric matrices, use `--format matrix`. Non-finite coordinates are
excluded with a warning.

`--bins 700` controls the raster resolution per axis (default 700; maximum
2000). PNG and PDF are supported; PDF embeds the density raster while keeping
axes and text as vectors. `--title "My cohort"` sets the title.
The raster preserves the full coordinate range, so extreme outliers can
compress the central cloud. Bin counts depend on the chosen resolution.

From R, source the script and call it on an open graphics device:

```r
source("scripts/plot-pca.R")
png("pca.png", width = 1800, height = 1500, res = 200)
result <- plot_pca("pcaone.eigvecs2", pcs = c(1, 2), bins = 700)
dev.off()
sum(result$counts)  # number of plotted samples
```

## SNP loadings of millions of variants

[`scripts/plot-loadings.R`](../guide/plotting.md#snp-loadings) plots the
loadings of each PC along the genome, as the largest |loading| in bins of
consecutive variants, so millions of variants plot in seconds: about 16 s for
6.6M SNPs and 12 PCs. It reads only the PCs asked for with `--pcs`, 8 bytes per
SNP per PC, e.g. about 0.5 GB for 6.6M SNPs and 10 PCs. Use it to find PCs
driven by a single region, such as the HLA, a centromere or an inversion,
before using them as covariates.

## Per-sample inbreeding

[`scripts/plot-inbred.R`](../small-n/hwe.md#plot-per-sample-results) draws the
`.inbred` output of `--inbreed 2` with log-count density rasters for its
scatter panels, so it handles hundreds of thousands of samples without
subsampling.

## Method details

`plot-pca.R` reads only the two selected PC columns. It pads their range by 2%
on each side, splits each axis into `--bins` equal intervals, and counts the
samples in each of the `bins x bins` cells in one pass. A cell's colour is
`log1p(count) / log1p(max count)` on a 256-colour scale, so a cell with one
sample is visible next to one with thousands; empty cells stay white. The
memory is the two columns plus the `bins x bins` counts, independent of how
the samples are distributed, and `plot_pca()` returns the counts so you can
check that every sample was plotted.
