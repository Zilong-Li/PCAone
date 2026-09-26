# Reproducing the `--evaladmix` benchmark

Everything in [evaladmix.md](evaladmix.md) is reproducible from two scripts in
[`scripts/benchmark/`](../scripts/benchmark/). Expect about ten minutes, most of
it compiling the comparison tools.

```bash
bash scripts/benchmark/run_all.sh  ~/evaladmix-benchmark  ./PCAone
Rscript scripts/benchmark/benchmark.R ~/evaladmix-benchmark
```

The first script fetches the data, builds the other methods and runs everything;
the second reads the outputs and prints every table. Both are idempotent — rerun
either without redoing the first. Both paths may be relative. `THREADS` (default
8) sets the threads for every tool; keep the default to match the numbers below
exactly, because RelateAdmix's estimates for unrelated pairs depend on it (RMSE
0.00054 with 8 threads, 0.00046 with 2; related pairs are unaffected).

## Nothing is simulated

The dataset is **not** generated here — the missing-data runs below only mask
some of its calls. It ships with
[relateAdmix](https://github.com/aalbrechtsen/relateAdmix) in its `data/`
directory and is used as-is:

| file | what it is |
|---|---|
| `smallPlink.{bed,bim,fam}` | 126 individuals, 104,290 autosomal SNPs, no missing genotypes |
| `smallPlink.2.P` | ancestral allele frequencies, ADMIXTURE output for K=2 |
| `smallPlink.2.Q` | admixture proportions, same run |

The individuals are simulated from two source populations with a known pedigree.
**ADMIXTURE is never run**: the `.P` and `.Q` files are provided, and every method
that needs `q̂, f̂` is handed the same ones, so no method is advantaged by a better
admixture fit. The PCA-based methods get one PC (`K-1 = 1`) from the same
genotypes.

## Where "truth" comes from

Not from any method's output. The pedigree is visible in the data — individuals
are ordered so relatives are consecutive even–odd pairs, ten pairs per class —
and `benchmark.R` reconstructs it by index:

| individuals (0-based) | relationship | theoretical φ |
|---|---|---|
| 0–5 | unrelated | 0 |
| 6–25 | duplicate / MZ | 0.5 |
| 26–45 | parent–offspring | 0.25 |
| 46–65 | full sib | 0.25 |
| 66–85 | half sib | 0.125 |
| 86–105 | first cousin | 0.0625 |
| 106–125 | second cousin | 0.015625 |

with the remaining 7,815 pairs unrelated. "Truth" throughout is the
**theoretical** kinship for the relationship, never an estimate. That is worth
stressing: an earlier version of this comparison used RelateAdmix's own output as
the reference, which flatters RelateAdmix and obscures that it attenuates badly
at half sibs and cousins.

## What gets run

`run_all.sh` does the following, all with the same PLINK files and the same
`.P`/`.Q`.

**RelateAdmix** — ML estimation of `(k0,k1,k2)` (built from source):

```bash
./relateAdmix -plink smallPlink -f smallPlink.2.P -q smallPlink.2.Q -P 8
```

**evalAdmix**, both estimators — requires v1.0 or later for `-method corrected`:

```bash
# leave-one-out frequency correction (EM), the default
./evalAdmix -plink smallPlink -fname smallPlink.2.P -qname smallPlink.2.Q \
            -P 8 -o evaladmix_em.corres
# analytic projection estimator of van Waaij et al. 2023
./evalAdmix -plink smallPlink -fname smallPlink.2.P -qname smallPlink.2.Q \
            -method corrected -P 8 -o evaladmix_proj.corres
```

**PCAone** — the implementation under test, plus three consistency runs:

```bash
PCAone -b smallPlink -k 1 -d 0 --evaladmix --maf 0.05 -o pcaone
PCAone -b smallPlink -k 4 -d 0 --evaladmix --evaladmix-k 1 --maf 0.05 -o pcaone_k4
PCAone -b smallPlink -k 1 -d 0 --evaladmix -m 0.002 -o pcaone_ooc   # out-of-core
PCAone -b smallPlink -k 1 -d 0 --evaladmix          -o pcaone_ic    # in-core
```

`-d 0` selects IRAM for a deterministic run. `--maf` is omitted from the
out-of-core run because PCAone rejects that combination.

**Missing genotypes** — `make_missing.R` writes five copies of the data with calls
set to missing, and PCAone (`-k 1 -d 0 --evaladmix --maf 0.05`) and evalAdmix EM
(same `.P`/`.Q`) run on each:

| design | what is missing |
|---|---|
| `mcar0.05`, `mcar0.1`, `mcar0.2` | each call with probability 5, 10, 20% |
| `varying` | each call with a per-sample probability drawn from U(0, 0.4) |
| `batch` | the samples split at random into two batches, each missing its own random 25% of the sites (`batch.batch` records the split) |

The seeds are fixed, so the masks are the same on every run.

**PCA + projection, R reference** — `evalPCA()` from
[popgenDK/evalPopStructure](https://github.com/popgenDK/evalPopStructure), called
by `benchmark.R`:

```r
pca <- makePCA(g, method = "standard", center = TRUE, scale = FALSE)
evalPCA(pca, k = K - 1)$corres
```

**PC-Relate** — optional, and the only step needing a package install:

```r
BiocManager::install("GENESIS")
```

`benchmark.R` skips the PC-Relate rows with a message if it is absent, so the
rest still runs. Per the design of the comparison, **PC-AiR is deliberately not
used**: the PCs come from plain PCA on all 126 individuals
(`training.set = NULL`), so the unrelated-training-set question is held constant
across methods. PC-Relate is run twice, with `small.samp.correct` on (the
default) and off, and at one to four PCs for the sensitivity table.

## Scales

Two conversions matter, and getting either wrong changes the conclusions:

- **evalAdmix returns the correlation of residuals, which estimates `2φ`.** All
  evalAdmix-family numbers are divided by two before comparison. PC-Relate and
  RelateAdmix are already on the `φ` scale.
- **PCAone codes genotypes as `{0, 0.5, 1}`**, i.e. `g/2`. This cancels inside
  `--evaladmix` because both `b̂` and `ĉ` become correlations, but it matters if
  you compute anything from the genotypes yourself.

## Checks built into the run

After the accuracy tables, `benchmark.R` prints three cross-checks that should
all pass:

| check | expected |
|---|---|
| `--evaladmix-k 1` (4-column eigvecs) vs a `-k 1` run | max diff 0 |
| in-core vs out-of-core | max diff 0 |
| PCAone vs the `evalPCA()` R reference | r = 0.999955 |

The first is a regression test for the
[`read_usv` bug](read-usv-fix.md): it reads a multi-column `.eigvecs` and
projects on its first column, which must equal using one PC directly. Before that
fix it differed by 8.8e-02.

## Expected output

The headline table, RMSE against theoretical kinship:

```
                     RMSE all RMSE unrelated RMSE related bias related
PCAone --evaladmix    0.00385        0.00386      0.00251     -0.00055
evalAdmix (EM)        0.00378        0.00378      0.00292     -0.00120
PC-Relate             0.00331        0.00330      0.00358      0.00055
evalAdmix (proj)      0.00382        0.00381      0.00429      0.00069
PCA + projection (R)  0.00387        0.00386      0.00443      0.00084
RelateAdmix           0.00089        0.00054      0.00812     -0.00571
PC-Relate (corr off)  0.00824        0.00822      0.00962     -0.00789
```

Numbers should match to the digits shown. The PCAone rows are deterministic;
GENESIS depends on `snpgdsPCA`, so PC-Relate may move in the last digit across
versions. (Rerun 2026-09-26 with GENESIS 2.32.0 and SNPRelate 1.36.1: every row
matched.)

Then the detection tables — KING degree bins, false positives among the unrelated
pairs, and AUC against them — and the missing-data table, which should read

```
                                     RMSE related est/truth degree right
5% at random | PCAone --evaladmix         0.00263   0.99148      1.00000
5% at random | evalAdmix (EM)             0.00301   0.98785      1.00000
10% at random | PCAone --evaladmix        0.00261   0.99328      1.00000
10% at random | evalAdmix (EM)            0.00302   0.98964      1.00000
20% at random | PCAone --evaladmix        0.00303   0.99213      0.96667
20% at random | evalAdmix (EM)            0.00337   0.98948      0.96667
0-40% by sample | PCAone --evaladmix      0.00325   0.99093      1.00000
0-40% by sample | evalAdmix (EM)          0.00340   0.98521      1.00000
25% by batch | PCAone --evaladmix         0.00328   0.99133      0.96667
25% by batch | evalAdmix (EM)             0.00331   0.98668      1.00000
```

`est/truth` is the mean ratio over the related pairs other than second cousins;
`degree right` counts all 60. See
[evaladmix.md](evaladmix.md#missing-genotypes) for what PCAone gave before it
imputed from the PCs and rescaled by sites in common.

## Caveats on scope

One dataset: 126 samples, K=2, discrete and well-separated source populations,
PLINK input, complete or with calls masked at random or by batch. That is the
best case for the admixture-model methods and says nothing about continuous
ancestry, misspecified `K`, or the regime where `O(N²)` memory matters. With ten
pairs per class, detection is saturated on complete data: every method calls
every degree right down to first cousins. The differences between methods here
are in the third or fourth decimal — see
[evaladmix.md](evaladmix.md#getting-the-number-of-pcs-wrong-is-the-dominant-risk)
for the one choice that is worth two orders of magnitude more than the choice of
method.
