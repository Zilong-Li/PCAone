# HWE and inbreeding

In a structured sample, a test of Hardy-Weinberg equilibrium (HWE) against one
allele frequency per site flags every differentiated site: the Wahlund effect
shows up as a deficit of heterozygotes. As in
[PCAngsd](https://github.com/Rosemeis/pcangsd), PCAone tests each site against
the individual allele frequencies, the allele frequency of each sample predicted
from the PCs. It reconstructs them from the `.eigvecs`, `.sigvals`,
`.loadings` and `.mbim` of a reference PCA, so first run

```shell
PCAone -b example/plink -k 3 -o pcaone --pcangsd
## alternative for PLINK input with missingness
## PCAone -b example/plink -k 3 -o pcaone --emu 
```

Then `--inbreed 1` tests HWE and estimates the inbreeding coefficient F of
each site. The results are written to `.hwe`, one row per site with the
columns `#ID`, `HWE_P`, `LRT` and `Inbreeding_coefficient`.

```shell
PCAone -b example/plink \
       --USV pcaone \
       --inbreed 1 \
       -o inbreed
```

The target must have the samples and sites of the reference run, in the same
order and with the same alleles; PCAone checks this against the `.mbim`. So
apply `--maf` in the reference run, not here. Choose `-k` to capture the
population structure; too few PCs leave a Wahlund effect at the differentiated
sites.

## Plot per-site results

`scripts/plot-hwe.R` reads the `.hwe` output from `--inbreed 1`, with columns
`#ID`, `HWE_P`, `LRT` and `Inbreeding_coefficient`. It needs R and data.table.

```shell
Rscript scripts/plot-hwe.R inbreed.hwe -o hwe.png
# Add genomic positions from the matching variant file:
Rscript scripts/plot-hwe.R inbreed.hwe --bim pcaone.mbim -o hwe.pdf
```

The plots show the HWE Q-Q curve, P-value histogram, and per-site F histogram.
With `--bim`, a fourth panel shows the largest -log10(P) in each genomic
position bin, with a Bonferroni line at `--alpha / valid variants` (alpha
defaults to 0.05). The BIM/MBIM must have exactly the same variant IDs and row
order as the HWE file. Chromosomes appear in their input order.

For large files, the Q-Q plot uses exact sorted P values at logarithmically
spaced ranks, including the most extreme value. The genome plot retains the
strongest signal per bin. `--max-points 5000` controls Q-Q ranks and genomic
bins. These summaries bound plotting work; all valid variants contribute to
histograms and summary counts. Reading uses O(M) memory and sorting O(M log M)
work for M variants.

Invalid P values are excluded with a warning. Zero and extremely small P values
are capped at `--p-floor 1e-300` for plotting, with their count reported;
significance counts use the original P values. Non-finite F values are omitted
from the F histogram. The script also writes `OUTPUT.summary.tsv`, reporting
valid/excluded/capped P counts, the Bonferroni cutoff, significant variant
count, and median finite F. The LRT column is required as part of the input
format but is not plotted. These are per-site results, distinct from the
per-sample `.inbred` output.

From R:

```r
source("scripts/plot-hwe.R")
d <- read_hwe("inbreed.hwe", bim = "pcaone.mbim")
pdf("hwe.pdf", width = 11, height = 8)
summary <- plot_hwe(d)
dev.off()
```

Sourcing defines the public functions without running the command-line code:
`read_hwe(file, bim = NULL, p_floor = 1e-300)` reads and validates results;
`plot_hwe(d, alpha = 0.05, title = "HWE diagnostics", max_points = 5000)`
draws on the current device and invisibly returns the summary table;
`save_hwe()` manages a PNG/PDF device and optionally saves the summary table.
The plotting function restores graphics settings after drawing.

```r
summary <- save_hwe(d, out = "hwe.pdf", alpha = 0.05,
                    width = 11, height = 8,
                    summary_out = "hwe-summary.tsv")
# Skip writing a summary file:
save_hwe(d, out = "hwe.png", res = 200, summary_out = NULL)
```

## Per-sample inbreeding

With `--inbreed 2` the same reference gives the inbreeding coefficient of each
sample instead, written to a file with suffix `.inbred`:

```shell
PCAone -b example/plink \
       --USV pcaone \
       --inbreed 2 \
       -o inbreed
```

| column  | meaning                                                                  |
|---------|--------------------------------------------------------------------------|
| FID IID | from the `.fam`, the `.psam` (FID 0 if it has none) or the BEAGLE header |
| N_SITES | sites with a call, or with a likelihood that is not flat                 |
| O_HET   | heterozygotes observed (posterior expectation for `-G`)                  |
| E_HET   | heterozygotes expected at F = 0                                          |
| F       | `1 - O_HET/E_HET`; NA for a sample with no data                          |

F is `1 - O_HET/E_HET`: positive for fewer heterozygotes than expected from
the sample's ancestry, negative for more. The output is the same in-core and
with `-m`, and for any `-n`. For PGEN input only the hard calls are used, as
for `--inbreed 1`: when the file holds dosages, the uncertain ones become
missing calls. Those are more often heterozygotes, which biases F upwards. The
target must have the samples and sites of the reference run, so apply `--maf`
in that run, not here.

### Plot per-sample results

`scripts/plot-inbred.R` reads the `.inbred` file from `--inbreed 2`:

```shell
Rscript scripts/plot-inbred.R inbreed.inbred -o samples.png
Rscript scripts/plot-inbred.R inbreed.inbred -o samples.pdf --title "My cohort"
```

Four panels show the F histogram, informative-site count histogram,
observed versus expected heterozygosity (`O_HET/N_SITES` versus
`E_HET/N_SITES`), and F versus informative-site count. The two scatter panels
use log-count density rasters to support hundreds of thousands of samples
without subsampling. Dashed reference lines mark equal heterozygosity and F=0.
The full ranges are shown, so extreme values can compress the main cloud.

Samples with non-finite F, zero informative sites or nonpositive/non-finite
expected heterozygosity are omitted from F panels and counted in the summary.
Site coverage includes all samples. The script writes `OUTPUT.summary.tsv`
with sample counts, median site count, median F, and 1st/99th F percentiles.
It does not designate samples for removal. `--breaks N` controls histogram
binning (default 100). R and data.table are required.

The sourceable functions do not run command-line code:

```r
source("scripts/plot-inbred.R")
d <- read_inbred("inbreed.inbred")
summary <- save_inbred(d, "samples.pdf")
# Draw on an existing device:
plot_inbred(d)
# Custom file dimensions, without writing a summary TSV:
save_inbred(d, "samples.png", width = 11, height = 8, res = 200,
            summary_out = NULL)
```

## Method details

### Individual allele frequencies

The reference PCA approximates the scaled genotypes by `U S V'`. Undoing the
scaling recorded in `.sigvals` gives the individual allele frequency of sample
`i` at site `j`, `pi_ij = f_j + (U S V')_ij * c_j`, bounded to
`[1e-4, 1 - 1e-4]`. The factor `c_j` is `sd_j / sqrt(ploidy)` for the default
standardization, 1 for an unscaled PCA and 1/2 for a PCAngsd run, which
decomposes centred dosages. A reference with `--scale 1` to `4` cannot be
mapped back to allele frequencies and is refused.

### Per-site HWE test (`--inbreed 1`)

The test is that of Meisner & Albrechtsen (2019, *Mol Ecol Resour* 19:1144).
Under HWE with a per-site inbreeding coefficient `F_j`, the genotype
probabilities of sample `i` are

```
P(g = 0) = (1 - pi)^2   + pi (1 - pi) F
P(g = 1) = 2 pi (1 - pi) (1 - F)
P(g = 2) = pi^2         + pi (1 - pi) F
```

with `pi = pi_ij`. `F = 0` is HWE, so the null model is nested in the
alternative. PCAone estimates `F_j` by EM, accelerated with SQUAREM (`--maxiter`,
`--tol-em`), bounded to `[-1, 1]`, and computes the likelihood ratio statistic

```
LRT = 2 (log L(F_j) - log L(0))
```

with the same floor of 1e-12 on the genotype probabilities in both models.
`HWE_P` is its p-value from a chi-square with 1 df.

For called genotypes the EM reaches `F = 1 - O/E` in one step, the moment
estimator that matches the heterozygote count `O` to its expectation `E`. It
is not the maximum likelihood estimate, so at some sites it fits worse than
`F = 0` and the statistic comes out negative: at 5.45% of the sites of
`example/plink.chr1` (400 samples, `-k 3`). These are reported as 0, never
significant, and PCAone logs their count. With the constrained maximum
likelihood estimate instead, the sites at p < 1e-6 on that data went from 280
to 278, at several times the run time, so the moment estimator is kept.

On a simulated cohort of 10,000 samples from 20 populations (chromosome 22,
6,455 sites at MAF >= 5%, `-k 3`), a test stratified by population, with no
Wahlund effect, finds 0.30% of the sites out of HWE. PCAone removes 0.05% at
p < 1e-6 and PCAngsd 0.12%, and their F correlate at 0.75. Neither `--emu` nor
`--pcangsd` is needed for this.

### Per-sample inbreeding (`--inbreed 2`)

The estimator is the one of PCAngsd's `--inbreed-samples`, i.e. the moment
estimator of `plink --het` with each sample's own allele frequency $\pi_{ij}$
in place of one frequency per site: $F_i = 1 - O_i / E_i$, where $O_i$ is the
number of heterozygotes observed and $E_i = \sum_j 2\pi_{ij}(1-\pi_{ij})$ the
number expected without inbreeding, over the same sites. For called genotypes
(`-b`, `-p`) this takes a single pass over the data. For genotype likelihoods
(`-G`) $O_i$ is the posterior expected number of heterozygotes, which depends
on $F_i$, so it is solved by EM (`--maxiter`, `--tol-em`). Missing calls and
flat likelihoods are left out, which gives the same F as keeping them, in fewer
iterations. F is bounded to [-1, 1] and can be negative (excess heterozygosity).
