# Plotting

Scripts for plotting PCAone output. They are in `scripts/` and need R >= 4.0
and the `data.table` package.

- [SNP loadings](#snp-loadings): `scripts/plot-loadings.R`
- [Sample PCA density](#sample-pca-density): `scripts/plot-pca.R`
- [Compare PCA runs](#compare-pca-runs): `scripts/compare-pca.R`
- [Per-site HWE](hwe.md#plot-per-site-results): `scripts/plot-hwe.R`
- [Per-sample inbreeding](hwe.md#plot-per-sample-results): `scripts/plot-inbred.R`
- [LD decay](#ld-decay): `scripts/plot-ld-decay.R`, `scripts/summarise_ld_r2bin.cpp`
  and the Nextflow workflow `workflows/ld.nf`

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

## Compare PCA runs

Compare winSVD and IRAM sample coordinates with:

```bash
Rscript scripts/compare-pca.R winsvd.eigvecs2 iram.eigvecs2 \
    --k 10 --labels winSVD,IRAM -o winsvd-vs-iram
# Numeric matrices require confirmed identical sample order:
Rscript scripts/compare-pca.R winsvd.eigvecs iram.eigvecs \
    --k 10 --row-order -o winsvd-vs-iram
```

The first K PCs are compared (default 10). Inputs with IDs are matched by
FID/IID, retaining their intersection and reporting excluded sample counts.
Duplicate IDs and non-finite coordinates cause an error. Use `--format-a ids`
or `--format-b ids` for headerless files with two ID columns under other names;
`matrix` selects numeric-only input. Both runs should use the same variants
and preprocessing for a solver comparison.

The script writes a three-panel PDF and three tables:

- `.pcs.tsv`: signed and absolute same-PC correlation, sign correction,
  and relative coordinate errors after sign alignment or an orthogonal
  Procrustes rotation of B toward A.
- `.angles.tsv`: principal angles between the two K-dimensional subspaces,
  ordered from best to worst agreement; these directions are not individual PCs.
- `.summary.tsv`: shared/excluded sample counts, mean squared principal-angle
  cosine (subspace overlap), largest angle, and overall Procrustes error.

All coordinates are centred and each PC is normalised to unit length, so
errors are scale independent. Zero error, zero angles and overlap one indicate
agreement. The correlation heatmap also reveals PC swaps. Similar eigenvalues
can allow PCs to rotate while spanning the same subspace, so inspect subspace
agreement alongside individual-PC correlation. Rotation errors describe the
normalised coordinates; principal angles measure subspaces independently of
column scaling. This does not compare eigenvalues or establish which solver
is more accurate relative to the original data matrix.

Work scales as O(N K^2), memory as O(N K); no sample-by-sample N x N matrix is
constructed. Only selected PCs and, when present, IDs are read. Large K can
still be expensive. From R, source the script and call
`compare_pca(file_a, file_b, k = 10)`; use `plot_comparison(result)` on an open
graphics device.

## SNP loadings

`scripts/plot-loadings.R` plots the SNP loadings of each PC along the genome.
Use it to see which PCs reflect genome-wide structure and which are driven by
a single region, such as a centromere, the HLA or an inversion. Such PCs are
often worth dropping, or fixing with LD pruning, before downstream analyses.

### Input

Run PCAone as usual (the loadings are written by default):

```bash
PCAone -b example/plink -k 10 -o pcaone
```

This writes `pcaone.loadings`, with one row per variant and one column per PC,
and `pcaone.mbim`, with the position of each variant. The script also reads
`pcaone.sigvals` and, by default, multiplies each PC’s loadings by its singular
value before plotting and reporting peak loadings. Use `--sigvals F` for a
custom singular-values file, or `--no-weights` to plot raw loadings. In R, use
`read_loadings(..., sigvals = "custom.sigvals")` or `weighted = FALSE`.

### Quick start

```bash
# all PCs in one plot
Rscript scripts/plot-loadings.R -p pcaone

# one row per PC
Rscript scripts/plot-loadings.R -p pcaone --mode panel --pcs 1-8

# colour PC1 and PC2, and colour each PC's peak by category
Rscript scripts/plot-loadings.R -p pcaone --highlight 1,2 \
    --groups "Population structure:1,2,4;HLA:6,12"

# zoom in on the HLA
Rscript scripts/plot-loadings.R -p pcaone --region 6:25M-35M --mode panel -o hla.pdf
```

The plot is written to `pcaone.loadings.png` unless `-o` is given. The
extension of `-o` sets the format: `.png`, `.pdf`, `.svg`, `.jpg` or `.tiff`.

### What is plotted

Each PC is summarised as the maximum |loading| in bins of `--window`
consecutive variants. By default the window gives about 4000 bins over the
plotted variants, e.g. 1,650 SNPs per bin for 6.6M SNPs. Bins never cross
chromosomes, and the maximum keeps every peak, so millions of variants plot in
seconds. Use `--window 1` to plot every variant, e.g. in a small region.

Raw loadings are about `1/sqrt(M)` for M SNPs; singular-value weighting changes
their scale. Small plotted values are shown in
multiples of a power of ten, e.g. `(x10^-3)`.

The x axis follows the variant order (`--xaxis index`), with chromosomes in
alternating shades. With `--xaxis bp` it uses physical positions, so gaps such
as centromeres take up space. A single chromosome or region uses `bp` by
default, labelled in Mb.

### Two layouts

**overlay** (default). All PCs are in one plot as grey lines, and the
`--highlight` PCs are drawn in colour. Each PC gets a dot with its number at
its highest bin. With `--groups`, the dots are coloured by category. This is a
compact overview, but dots overlap when many PCs peak in the same region.

**panel** (`--mode panel`). One row per PC with its own y axis (`--fixed-y` for
a shared one), the peak position labelled (e.g. `chr6:32.59 Mb`), and the group
name in the right margin. Use it when overlay dots pile up, or to inspect a
region.

### Groups

`--groups` puts PCs into categories. Give them inline, separated by `;`:

```bash
--groups "Population structure:1,2,4;Centromere:3,5,8-11;HLA:6,12"
```

or as a file with one group per line, either `label: PCs` or `PC label`:

```
# label: PCs
Population structure: 1,2,4
Centromere: 3,5,8-11
HLA: 6,12
```

The categories of the original 40-PC loadings figure are:

```bash
--groups "Population structure:1,2,4;Centromere:3,5,8-11,13,14,19,20,22,29,32,34,35,37;HLA:6,12,15,16,27,30,38;Other structure:7,17,18,23-26,28,31,33,36,39,40;Inversion:21"
```

Colours come from a fixed 8-colour, colour-blind-checked palette. Groups get
colours first, in the order given, then the highlighted PCs, so at most 8
groups and highlighted PCs together. In the overlay, PCs that are in no group
get grey dots, listed as "Ungrouped".

### Peak table

The script prints each PC's top variant (also marked in the plot).
`--peaks-out FILE` saves it as a tsv.

| column      | meaning                                                          |
|-------------|------------------------------------------------------------------|
| `PC`        | the PC                                                           |
| `chr`, `snp`, `bp` | the variant with the largest \|loading\|                  |
| `row`       | its row in `.loadings` / `.mbim` (1-based)                       |
| `loading`   | its loading, with sign                                           |
| `chr_share` | share of the PC's sum of squared loadings on that chromosome     |
| `group`     | its `--groups` category                                          |

`chr_share` helps sort PCs into groups. For a PC driven by genome-wide
structure it is close to the chromosome's share of the SNPs, e.g. about 0.06
for chr6. A PC driven by one region has a much larger share on that
region's chromosome. It is `NA` when only one chromosome is plotted.

### All options

| option | meaning | default |
|---|---|---|
| `-p, --prefix P` | read `P.loadings` and `P.mbim` (or `P.bim`) | |
| `-l, --loadings F` | loadings file, instead of `-p` | |
| `-b, --bim F` | `.mbim`/`.bim` with the same variants; without one the x axis is the variant index | `P.mbim` |
| `--pcs S` | PCs to plot, e.g. `1-10` or `1,3,5` | all |
| `--chr S` | chromosomes to keep, e.g. `1-22` or `6` | all |
| `--region R` | one region, e.g. `6:25M-35M` or `chr6:25000000-35000000` | |
| `-o, --out F` | output file | `P.loadings.png` |
| `--mode M` | `overlay` or `panel` | `overlay` |
| `--highlight S` | overlay: PCs drawn in colour | none |
| `--groups G` | PC categories, inline or a file | none |
| `--window N` | variants per bin | about 4000 bins |
| `--xaxis X` | `index` or `bp` | `bp` for one chromosome, else `index` |
| `--fixed-y` | panel: same y axis for all PCs | |
| `--no-peaks` | do not mark each PC's peak | |
| `--title T` | plot title | |
| `--width W`, `--height H` | size in inches | 11 x 4.5 (overlay), 11 x (0.9 + 1 per PC) (panel) |
| `--res N` | resolution of png/jpg/tiff | 150 |
| `--cex X` | scale text and markers | 1 |
| `--peaks-out F` | also save the peak table | |

`Rscript scripts/plot-loadings.R -h` prints the same list.

### From R

Sourcing the script defines the functions without running it:

```r
source("scripts/plot-loadings.R")
d <- read_loadings("pcaone.loadings", "pcaone.mbim", pcs = 1:10, chr = 1:22)
pdf("loadings.pdf", 11, 4.5)
pk <- plot_loadings(d, highlight = 1:2, groups = list(HLA = c(6, 12)))
dev.off()
pk                  # the peak table
loading_peaks(d)    # the same, without plotting
```

`plot_loadings()` takes the same settings as the command line: `mode`,
`highlight`, `groups`, `window`, `xaxis`, `peaks`, `fixed_y`, `title` and
`cex`. It draws on the current device, so in overlay mode it can be one
panel of a larger figure, e.g. after `par(mfrow = c(2, 1))`.

### Notes

- Variants must be grouped by chromosome, as PCAone writes them. Otherwise the
  script warns and labels each block of a split chromosome separately.
- The selected PCs are read into memory: 8 bytes per SNP per PC, e.g. about
  0.5 GB for 6.6M SNPs and 10 PCs. `--pcs` reads only the PCs asked for.
- Speed: about 16 s for 6.6M SNPs and 12 PCs.
- The `.loadings` and `.mbim` must come from the same run. The script stops if
  their row counts differ.

## LD decay

`scripts/plot-ld-decay.R` plots the mean r2 of SNP pairs against their
distance. Comparing the standard LD with the ancestry-adjusted LD shows how
much of the long-range LD comes from population structure: in a structured
sample the standard curve levels off well above the r2 of unlinked SNPs, while
the adjusted curve decays to it.

### Input

Run PCAone twice with `-R/--print-r2`, once with `-P` for the ancestry-adjusted
LD and once with `--ld-stats 1` for the standard LD:

```bash
PCAone -b plink -k 5 -o pcs                                    # the PCs
PCAone -b plink -P pcs -R --ld-bp 1000000 -o adj               # adjusted LD
PCAone -b plink --ld-stats 1 -R --ld-bp 1000000 -o std         # standard LD
```

Each run writes all pairs within `--ld-bp` to `<out>.ld.gz`. The number of
pairs grows with the SNP density times the window, so thin dense data first,
e.g. with `plink --thin`. For reference, 104k SNPs in 1 Mb windows give 25M pairs
and a 140 MB file.

The script also reads `plink --ld` output (use `--ld-window-r2 0`, or plink
drops the pairs below 0.2).

### Quick start

```bash
# both curves in one plot, with the r2 expected without LD, 1/(N-1)
Rscript scripts/plot-ld-decay.R adj.ld.gz std.ld.gz --labels Adjusted,Standard \
    -n plink.fam -o ld-decay.png

# save the bins, then re-plot from them in a second
Rscript scripts/plot-ld-decay.R adj.ld.gz std.ld.gz --save-bins -o ld-decay.png
Rscript scripts/plot-ld-decay.R adj.decay.tsv std.decay.tsv --labels Adjusted,Standard \
    --ylog --by-chr -o ld-decay-log.pdf
```

Each input is one curve. It prints a summary per curve:

```
    curve    pairs first_bin_r2 plateau_r2 half_decay_bp unlinked_r2
 Adjusted 35703027       0.4093   0.009298          9304    0.008986
 Standard 35703027       0.4116   0.019970          9334    0.018930
```

This example uses the 1000 Genomes panel of the tutorial: 120 samples from
6 populations, every 5th SNP of chr20–22, 5 PCs and 1 Mb windows. With
1/(N-1) = 0.0084, the adjusted plateau is close to the no-LD value, while the
standard one is more than twice as high.

| column | meaning |
|---|---|
| `pairs` | pairs in the plotted bins |
| `first_bin_r2` | mean r2 of the closest pairs, below `--min` |
| `plateau_r2` | mean r2 of the bins that reach into the last 20% of the distance range |
| `half_decay_bp` | distance where r2 falls halfway from the first bin to the plateau (log-interpolated) |
| `unlinked_r2` | r2 of unlinked pairs, from `--baseline` or from inter-chromosomal pairs in the input |

With `--correct`, the summary uses corrected r2.

### What is plotted

The pairs are binned by distance: one bin for pairs closer than `--min`
(default 1 kb), then `--bins` log-spaced bins (default 40) up to `--max`. By
default `--max` is the largest distance in the first `.ld.gz`, rounded up to
one significant digit, i.e. the `--ld-bp` of the run. All inputs share the same
bins. `--linear` uses equal-width bins from 0 and a linear x axis.

Each point is one bin: the mean distance of its pairs against their mean r2.
The chromosomes are pooled, weighted by their number of pairs. `--chr` pools
only the chosen chromosomes, and `--by-chr` also draws each chromosome as a
thin line.

The r2 of a sample has an upward bias, about 1/N for unlinked SNPs.

- `-n N` draws the expected r2 without LD, 1/(N-1), as a dotted line. N can be
  a number or a `.fam`/`.psam` file.
- `--correct` removes the bias instead, as in LD score regression:
  r2 - (1 - r2) / (N - 2).

`--ylog` gives a log y axis, which separates the plateaus better.

### Unlinked pairs

`--baseline 0.0090,0.0189` gives the r2 of unlinked pairs, one value per input.
The script draws these as short lines to the right of the plot, labelled
"unlinked". If an input contains pairs on different chromosomes (e.g.
`plink --ld --inter-chr`), their mean r2 is used without `--baseline`.

PCAone only pairs SNPs on the same chromosome. The `cross` step of the
workflow below gets unlinked pairs anyway: it gives the SNPs of two chromosomes
random chromosome labels, lets plink re-sort them, and averages the r2 of the
pairs whose SNPs came from different chromosomes.

### Large files: summarise_ld_r2bin

`plot-ld-decay.R` reads the whole `.ld.gz` into memory: 25M pairs take about
14 s and 2.4 GB. For larger files, bin them with `summarise_ld_r2bin`. It
streams the file (25M pairs in 9 s, 10 MB of memory) and writes the same bins,
byte for byte, as `--save-bins`:

```bash
make -C scripts summarise_ld_r2bin                   # needs zlib
scripts/summarise_ld_r2bin -i adj.ld.gz --max 1000000 -o adj.decay.tsv
scripts/summarise_ld_r2bin -i std.ld.gz --max 1000000 -o std.decay.tsv
Rscript scripts/plot-ld-decay.R adj.decay.tsv std.decay.tsv --labels Adjusted,Standard -n plink.fam
```

It takes `--bins`, `--min`, `--max` (default 1000000; set it to the run's
`--ld-bp`) and `--linear`, as in the R script. `-i -` reads from stdin.

The bin table (tsv) has one row per chromosome and bin: `chr`, `lo`, `hi` (bin
edges), `n` (pairs), `mean_dist`, `mean_r2` and `sd_r2`. Pairs on different
chromosomes form one row, with chr `inter` and `NA` edges.

### All options

| option | meaning | default |
|---|---|---|
| `-o, --out F` | `.png`, `.pdf`, `.svg`, `.jpg` or `.tiff` | `ld-decay.png` |
| `--labels L` | comma-separated curve names | input file names |
| `--bins N` | log-spaced bins between `--min` and `--max` | 40 |
| `--min D` | upper edge of the closest-pairs bin, bp (`1k`, `2.5M` also work) | 1000 |
| `--max D` | upper edge of the last bin, bp | largest distance, rounded up |
| `--linear` | equal-width bins and a linear x axis | |
| `--save-bins` | write the bins of each `.ld.gz`, e.g. `adj.ld.gz` to `adj.decay.tsv` | |
| `--chr S` | chromosomes to pool, e.g. `1-22` | all |
| `-n, --nsamples N` | sample size or `.fam`/`.psam`, one or one per input | |
| `--correct` | correct r2 for sample size (needs `-n`) | |
| `--baseline V` | r2 of unlinked pairs, one per input | inter-chromosomal pairs |
| `--by-chr` | also draw each chromosome | |
| `--ylog` | log y axis | |
| `--ymax Y` | upper y limit | |
| `--legend P` | legend position, e.g. `bottomleft` | `topright` |
| `--title T` | plot title | |
| `--width W`, `--height H` | size in inches | 6 x 4.5 |
| `--res N` | resolution of png/jpg/tiff | 150 |
| `--cex X` | scale text | 1 |
| `--table F` | also write the plotted curves (tsv) | |

The binning options do not apply to bin tables, which keep their bins.

### From R

```r
source("scripts/plot-ld-decay.R")
adj <- ld_bins("adj.ld.gz")              # or a .decay.tsv
std <- ld_bins("std.ld.gz", max = 1e6)
pdf("ld-decay.pdf", 6, 4.5)
curves <- plot_ld_decay(list(Adjusted = adj, Standard = std), nsamples = 120)
dev.off()
decay_summary(curves)
write_bins(adj, "adj.decay.tsv")
```

`ld_bins()` takes `bins`, `min`, `max` and `linear`. `plot_ld_decay()` takes
the plot options above: `chr`, `nsamples`, `correct`, `baseline`, `by_chr`,
`ylog`, `ymax`, `title`, `cex` and `legend_pos`. `pool_bins()` gives one pooled
curve and `correct_r2()` the correction.

### Nextflow workflow

`workflows/ld.nf` runs the whole comparison for one or more PLINK data sets:
QC and thinning, the PCA, standard and adjusted LD, the bins and the plots.

```bash
make -C scripts summarise_ld_r2bin
nextflow run workflows/ld.nf --data data --pops giraffe --K 5,10 --run_step plot
```

It needs PCAone and plink 1.9 in `PATH`, and Rscript with data.table. It was
tested with Nextflow 25.04.

| parameter | meaning | default |
|---|---|---|
| `--data` | folder with `<pop>.{bed,bim,fam}` | `data` |
| `--pops` | data sets to run, comma-separated | `giraffe` |
| `--K` | numbers of PCs for the adjusted LD, comma-separated | `10` |
| `--run_step` | `curve`, `cross` or `plot` (both) | `curve` |
| `--maf`, `--thin` | plink `--maf` and `--thin` before the analyses | 0.05, 0.05 |
| `--ld_bp` | LD window, bp | 1000000 |
| `--ld_bins` | log-spaced distance bins | 40 |
| `--perm_chr` | chromosomes shuffled for the unlinked pairs | `1,2` |
| `--correct` | correct r2 for sample size (N from the `.fam`) | `true` |
| `--results` | output folder | `results` |

Outputs go to `<results>/<pop>/thin_<thin>/<run_step>/`:

- `curve`: `adj.<K>.ld.gz`, `std.ld.gz`, their `.decay.tsv` bins, and `ld_curve_k<K>.png`.
- `cross`: `cross_mean_adj_<K>.txt` and `cross_mean_std.txt`, the r2 of unlinked pairs.
- `plot`: all of these, and `ld_combine_k<K>.png` with the unlinked r2 on the right.

The adjusted LD of the unlinked pairs uses the PCs of the whole data set, the
same ones as the curve.
