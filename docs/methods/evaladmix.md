# `--evaladmix`: correlation of residuals from the top PCs

`--evaladmix` computes the evalAdmix statistic (Garcia-Erill & Albrechtsen 2020,
*Mol Ecol Resour* 20:936) using the analytic projection estimator of van Waaij
et al. (2023, *Genetics* 225:iyad157), with **PCA rather than an admixture model**
as the front-end.

```bash
PCAone -b plink -k <K-1> --maf 0.05 -o pcs
PCAone -b plink -P pcs --evaladmix --maf 0.05 -o out
```

writes

- `out.corres` — the N x N correlation of residuals, and
- `out.kinship` — the same divided by 2, which is the kinship scale,

both with sample IDs from the `.fam` or `.psam` on the header line.
`-k/--pc` selects how many of the computed PCs enter the projection, so one
PCA run can produce many PCs for other purposes while the statistic uses `K-1`.
The analysis reads `pcs.eigvecs` without rerunning PCA; `--read-U path` can
supply that file directly. By default it uses all reference columns; `-k`
selects a leading subset, as in the other two-stage analyses (`-P/--USV`).
The reference must contain the same samples in the same order as the genotype
input. The IIDs in `pcs.eigvecs2`, which a PCA of PLINK or PGEN input writes
beside `pcs.eigvecs`, are compared with the `.fam`/`.psam`, and a mismatch is an
error. Without that file only the number of rows is checked, with a warning.
Apply the desired `--maf` filter
again in the analysis stage.

For biobank-scale cohorts, where an `N x N` matrix does not fit,
`--evaladmix-kin <cutoff>` writes only the pairs above a kinship cutoff, plus an
unrelated set, in memory bounded by `-m`. See
[Biobank scale](#biobank-scale---evaladmix-kin).

`--evaladmix-ibd` adds the probabilities of sharing 0, 1 and 2 alleles IBD
(`k0`, `k1`, `k2`), which tell parent–offspring from full sibs. See
[IBD sharing probabilities](#ibd-sharing-probabilities---evaladmix-ibd).

## What it estimates

Predict each genotype from ancestry alone; whatever is left over is the
residual. Two individuals who share recent ancestors deviate from their
predictions in the same direction, so the correlation of their residuals
measures relatedness. Since a genotype is the sum of two alleles,

```
Cov(g_i, g_j) = 4 phi_ij pi (1-pi),    Var(g_i) = 2 pi (1-pi)
```

so the correlation of residuals estimates `2 * phi`. Near-zero entries also mean
the model fits, which is what evalAdmix was originally written to test — the two
readings are the same number.

With `P` the projection onto `[PC_1..PC_k, 1]` and `R = G(I-P)` the residuals,

```
corres = bhat - chat
bhat   = cov2cor( Rtilde' Rtilde ),          Rtilde = column-centred R
chat   = cov2cor( (I-P) Dhat (I-P) ),        Dhat = diag(mean heterozygosity)
```

`chat` is the correlation the fitted model implies on its own; subtracting it
removes the structure that the projection induces even when nobody is related.

**Use `K-1` PCs.** The intercept is included alongside them, so `K-1` PCs span
the same dimension as `Q` for a `K`-population admixture model. This is the single
most important setting — see
[Getting the number of PCs wrong](#getting-the-number-of-pcs-wrong-is-the-dominant-risk).

## Why no residual matrix is needed

With missing calls imputed (see [Missing genotypes](#missing-genotypes)),
`Rtilde' Rtilde` can be written without ever forming `R`:

```
Rtilde' Rtilde = (I-P) [ G'G - M gbar gbar' ] (I-P)
```

so the statistic depends on the data only through three quantities:

- `G'G`, the N x N genotype Gram matrix,
- `gbar`, the per-sample mean genotype,
- `Dhat`, the per-sample mean heterozygosity.

All three are accumulations over sites, so **one streaming pass** suffices.
Memory is one `N x N` matrix, independent of the number of sites, which is what
makes the approach usable out-of-core. Cost is one `O(MN^2)` accumulation — the
same order as the covariance matrix a PCA already forms. Everything after it is
`O(N^2 k)`: `P` has rank `k+1`, so each product with `I-P` is a rank-`k+1`
correction and `I-P` itself is never formed.

Measured with 2 threads on simulated genotypes, whole runs (PCA included, so
peak memory includes the PCA's), against the first implementation, which
multiplied by a dense `N x N` `I-P` and held four `N x N` matrices:

| samples x sites | before | now |
|---|---|---|
| 4,000 x 20,000 | 24.2 s, 1.15 GB | 10.1 s, 0.90 GB |
| 4,000 x 20,000, 5% missing | 24.7 s, 1.15 GB | 17.4 s, 1.16 GB |
| 8,000 x 10,000 | 109.8 s, 2.67 GB | 17.8 s, 1.65 GB |

At 8,000 samples, everything after the pass took 95 s and now takes 3 s, most
of it writing the two 600 MB files. The output is byte-identical on complete
data. Missing genotypes double the pass, because of the pair counts.

The identity was verified numerically against forming the residuals directly
(agreement to 1e-15).

## Accuracy

On the 126-sample example data distributed with
[relateAdmix](https://github.com/aalbrechtsen/relateAdmix) (K=2, 102,185 sites
after `--maf 0.05`, a known pedigree with ten pairs per relationship class):

| relationship | truth | `--evaladmix` |
|---|---|---|
| duplicate | 0.5 | 0.5000 |
| parent–offspring | 0.25 | 0.2496 |
| full sib | 0.25 | 0.2475 |
| half sib | 0.125 | 0.1227 |
| first cousin | 0.0625 | 0.0622 |
| second cousin | 0.0156 | 0.0177 |
| unrelated (7815 pairs) | 0 | −0.0016 |

## Detecting related pairs

The tables above measure how close the estimates are. A relatedness screen asks
something coarser: is the pair related, and to what degree? With the KING degree
bins (Manichaikul et al. 2010), extended to a 5th-degree bin, 0.0110–0.0221, for
the second cousins:

| method | degree right, duplicate … first cousin | second cousins | unrelated called ≥ 5th degree | max unrelated |
|---|---|---|---|---|
| **`PCAone --evaladmix`** | 50/50 | **10/10** | 0 | 0.0093 |
| evalAdmix (EM) | 50/50 | 10/10 | 0 | 0.0092 |
| evalAdmix (proj) / PCA + projection (R) | 50/50 | 10/10 | 0 | 0.0092 / 0.0093 |
| RelateAdmix | 50/50 | 7/10 | 0 | 0.0038 |
| PC-Relate | 50/50 | 7/10 | 1 of 7815 | 0.0120 |
| PC-Relate, `small.samp.correct=FALSE` | 50/50 | 8/10 | 4 of 7815 | 0.0132 |

No unrelated pair reaches the 4th-degree bin under any method, and every method
separates every class from the unrelated pairs (AUC = 1, except 0.997 for second
cousins under uncorrected PC-Relate). **On complete data this benchmark is
saturated:** it cannot rank the methods on detection, only on the second-cousin
boundary, where RelateAdmix's attenuation and PC-Relate's scatter cost three pairs
each. What moves the calls is missing data — see
[Missing genotypes](#missing-genotypes).

## How it compares with the other evalAdmix routes

The same statistic can be reached four ways (RMSE over the 60 related pairs):

| route | RMSE, related | needs ADMIXTURE | genotypes in RAM |
|---|---|---|---|
| ADMIXTURE + evalAdmix EM | 0.00292 | yes | yes |
| ADMIXTURE + projection | 0.00429 | yes | yes |
| PCA (K−1) + projection, R | 0.00443 | no | yes |
| **`PCAone --evaladmix`** | **0.00251** | **no** | **no** |

Accuracy is essentially identical across all four — correlations with the truth
differ in the fourth decimal. Most of the RMSE spread is the `[-1,1]` clip, which
this implementation applies and the others leave to the user: apply it uniformly
and the three projection routes land on 0.00251–0.00255 against the EM route's
0.00292. So the margin above is a sensible default, not a better estimator.

The real differences are cost and dependencies. The EM route refits `F` with each
individual excluded, so its work grows with sample size on top of the per-pair
cost. The two PCA routes need no admixture run at all — usually the dominant cost
of a real analysis. Only this one avoids holding the genotype matrix in RAM.

All four overshoot equally at second cousins (~1.13x), which is a property of the
statistic near the noise floor rather than of any front-end.

## How it compares with PC-Relate and RelateAdmix

These are genuinely different methods, not other routes to the same statistic.
**Read the split, not the total.** 99.24% of the 7875 pairs are unrelated, so an
RMSE over all pairs is very nearly a measurement of performance on unrelated
pairs alone:

| method | RMSE, all | RMSE, unrelated | **RMSE, related** | bias, related |
|---|---|---|---|---|
| RelateAdmix | **0.00089** | **0.00054** | 0.00812 | −0.00571 |
| PC-Relate (GENESIS) | 0.00331 | 0.00330 | 0.00358 | +0.00055 |
| evalAdmix (EM) | 0.00378 | 0.00378 | 0.00292 | −0.00120 |
| PCA + projection (R) | 0.00387 | 0.00386 | 0.00443 | +0.00084 |
| **`PCAone --evaladmix`** | 0.00385 | 0.00386 | **0.00251** | −0.00055 |
| PC-Relate, `small.samp.correct=FALSE` | 0.00824 | 0.00822 | 0.00962 | −0.00789 |

**The ranking inverts between the columns.** On overall RMSE, RelateAdmix looks
best by a factor of four. On the 60 pairs that actually have relatedness to
estimate it is the *worst* of the five, at more than three times `--evaladmix`.
Its overall win is an artefact of the constrained ML: `(k0,k1,k2)` must lie on
the simplex, so `phi >= 0` by construction. Among unrelated pairs it never
returns a negative value and 93.6% have `k1 < 1e-3` — the estimates are pinned at
the boundary, which happens to be the truth for those pairs. PC-Relate, an
unconstrained estimator, scatters symmetrically about zero.

So: **PC-Relate is better at "are these two unrelated?"; `--evaladmix` is better
at "how related are they?"**. Neither margin is large — both correlate with the
truth at 0.988.

What separates them is not accuracy. PC-Relate is the only one that returns
`(k0,k1,k2)` and inbreeding coefficients, which matters because parent–offspring
and full sibs both have `phi = 1/4` and no kinship estimator can tell them apart.
`--evaladmix` is the only one that needs no admixture run and never holds the
genotype matrix in RAM.

One aside worth knowing if you use PC-Relate: nearly all of its accuracy comes
from `correctKin`, the post-hoc `small.samp.correct=TRUE` step, which regresses
kinship on the PC outer products over the effectively-unrelated pairs and
subtracts the fit. Disable it and PC-Relate lands in the last row above. The raw
estimator is biased downward because the pair's own genotypes sit inside the
fitted allele frequencies; PC-Relate does not avoid that, it measures and
subtracts it afterwards.

Every number in these tables is produced by
[`scripts/benchmark/`](https://github.com/Zilong-Li/PCAone/tree/main/scripts/benchmark) — see
[reproducing-the-benchmark.md](reproducing-the-benchmark.md).

## Missing genotypes

Missing calls are imputed to the site mean by the readers, which keeps the
projection defined. Left there, they cost the estimator its calibration: an
imputed call carries no relatedness, so the covariance of a pair runs over the
sites both are genotyped at while each variance runs over its own, and the
estimate shrinks by `n_ij / sqrt(n_i n_j)` — by the missing fraction when calls
are missing at random. `--evaladmix` therefore

1. **replaces each missing call by its fit from the PCs**, `(QQ'g)_i`. At the
   site mean a missing call keeps a residual `f − π_i`, its ancestry deviation,
   which two samples share when they miss the same sites — a genotyping batch;
2. **rescales each pair by `sqrt(n_i n_j) / n_ij`**, the pair counts accumulated
   in the same pass. Pairwise, not per sample: with a per-sample factor
   `1/sqrt(o_i o_j)` instead, pairs in the same batch came out 23% too high in
   a batch design like the one below.

The benchmark data with calls masked (`run_all.sh` step 6), `K−1 = 1` PC,
related pairs:

| missing calls | RMSE, before | **RMSE, now** | evalAdmix (EM) | degree right, before → now (EM) |
|---|---|---|---|---|
| none | 0.00251 | 0.00251 | 0.00292 | 60/60 → 60/60 (60/60) |
| 5% at random | 0.01215 | **0.00263** | 0.00301 | 60 → 60 (60) |
| 10% at random | 0.02546 | **0.00261** | 0.00302 | 60 → 60 (60) |
| 20% at random | 0.05213 | **0.00303** | 0.00337 | 59 → 58 (58) |
| 0–40% by sample | 0.05886 | **0.00325** | 0.00340 | **49 → 60** (60) |
| 25% by batch | 0.04697 | **0.00328** | 0.00331 | **54 → 58** (60) |

Before, the estimate was 0.79 of the truth at 20% missing, and in the batch
design 0.73 of it for pairs in different batches but 0.98 within one — the
missing fraction a pair does not share. The degree calls failed accordingly:
with 0–40% missing per sample, one duplicate pair was called 1st degree, five
parent–offspring and full-sib pairs 2nd degree, and two pairs each of half sibs
and first cousins a degree too distant. Now every design is within 0.99
of the truth and as accurate as evalAdmix's pairwise-complete EM, and every
remaining miss is a second cousin at 0.0223–0.0233, just over the 0.0221 bin
edge, as EM's are at 20%. RMSE on unrelated pairs rises as it should with fewer
sites, and matches EM's within 6% (0.00389–0.00425 against 0.00382–0.00404).

The cost is a second `N x N` matrix for the pair counts and about twice the time
for the pass; complete data allocates nothing and is unchanged. `--evaladmix-kin`
counts the same pairs from the lists of missing calls instead, which costs
next to nothing at biobank call rates. Sites whose
frequency is exactly 0.5 (0.1–0.2% here) are left out of both steps: there a
heterozygous call and a missing one look the same once centred.

## Getting the number of PCs wrong is the dominant risk

Far more important than the choice of method. With `K=2`, one PC is correct:

| PCs used | 1 (correct) | 2 | 3 | 4 |
|---|---|---|---|---|
| `--evaladmix` / PCA + projection | **0.00443** | 0.03035 | 0.05001 | 0.06917 |
| PC-Relate, for comparison | **0.00358** | 0.03474 | 0.06011 | 0.09459 |
| PC-Relate, mean at duplicates | **0.4972** | 0.5395 | 0.5880 | 0.6561 |

One PC too many costs a factor of seven to ten — two orders of magnitude more
than any difference between methods. With `K` real populations the later PCs fit
*relatedness* rather than ancestry, and projecting them out distorts the very
residual structure the estimator depends on.

Note that PC-Relate is **not** more robust to this, despite needing no explicit
`K`: choosing the number of PCs is the same decision in different clothes, and it
degrades slightly faster. The genuine advantage of a PC front-end is that it does
not assume *discrete* ancestral populations, not that it saves you a decision.

Note also how it fails: the damage concentrates in the closest pairs — duplicates
inflate to 0.66 while parent–offspring stays at 0.251 and first cousins at 0.066
— so it will not announce itself in a summary statistic.

**Validation.** Against `evalPCA()` in
[popgenDK/evalPopStructure](https://github.com/popgenDK/evalPopStructure):
r = 0.999955 over all 7875 pairs, identical RMSE on related non-duplicate pairs.
The only entries differing by more than 1.6e-04 are the ten duplicate pairs,
where this implementation applies the `[-1,1]` clip and the R reference does not.
In-core and out-of-core runs, and `-k 1` on a 4-PC reference against a 1-PC run, agree
exactly.

## Biobank scale: `--evaladmix-kin`

The dense output does not scale: at N = 245,000 (All of Us srWGS) the Gram
matrix alone is 480 GB of RAM, and each text file is about 540 GB.
`--evaladmix-kin <cutoff>` computes the same statistic within `-m` and writes
only the pairs whose kinship reaches the cutoff:

```bash
PCAone -b cohort -k 16 -m 64 -o pcs
PCAone -b cohort -P pcs -k <K-1> --evaladmix --evaladmix-kin 0.0442 -m 64 -o rel
```

writes

- `rel.kin0`, one line per pair with kinship >= the cutoff, ordered by the
  first sample and then the second, both in file order:
  `#FID1 IID1 FID2 IID2 NSNP KINSHIP` (`#IID1 IID2 NSNP KINSHIP` for a `.psam`
  without FID). `KINSHIP` is the entry `.kinship` would have, to the printed
  digit. `NSNP` is the number of sites where both samples are genotyped. A
  missing call at a site whose frequency is exactly 0, 0.5 or 1 cannot be
  detected (see [Missing genotypes](#missing-genotypes)), so it counts as
  genotyped;
- `rel.unrelated`, a maximal set of samples without a pair at or above
  `--evaladmix-unrelated` (default: the cutoff), as `FID IID` lines for
  `--keep`. It uses the greedy rule of Hail's `maximal_independent_set` (which
  All of Us used) and plink2's `--king-cutoff`: drop the sample with the most
  relatives left until none has any, breaking ties by more missing calls, then
  take back every dropped sample whose relatives were all dropped. A sample
  genotyped at none of the sites has no estimate, so it is left out, with a
  warning.

The log reports each stripe with an estimate of the time left, and counts the
pairs per KING degree bin. It also compares the scatter of unrelated pairs with
chance. Kinship of unrelated pairs has an sd of about `0.5 / sqrt(Meff)` when
the PCs fit, where `Meff = (sum v)^2 / sum v^2` with `v = f(1-f)` per site
(divided by the mean share of the sites a sample is genotyped at).
The pairs below 0 contain no relatives, so their root mean square measures the
scatter. On the simulations below the two agree to within 0.3% (0.00761
against 0.00761 at N = 20,000). That gives two warnings:

- the scatter is more than 1.3 times chance: the PCs leave structure, as with
  one PC too few (1.47 times) — see
  [Getting the number of PCs wrong](#getting-the-number-of-pcs-wrong-is-the-dominant-risk);
- chance alone should put more than a tenth as many unrelated pairs above the
  cutoff as were written: the cutoff is within the noise of the sites. At N =
  8,000 with 5,000 sites, a cutoff of 0.0221 is 2.9 sd, and chance predicts
  60,285 of the 62,505 pairs written.

A cutoff of 0.0442 keeps
3rd degree and closer; All of Us used 0.1 for its unrelated set. Choose the
cutoff with the rerun in mind, since a rerun repeats the whole computation:
for example, write the pairs from 0.0221 and build the unrelated set at 0.1.

**How.** In `corres_ij = sb_i sb_j S_ij - (L R')_ij` only `S_ij` is a pair
quantity. Every other term is per sample, including `T = S Q = G(G'Q) - M gbar (gbar'Q)`,
so one pass over the genotypes collects them in `O(NMr)`. The pairs are then
computed in stripes of samples, each holding `S(j, i)` for the stripe's samples
`i` and every `j > i`. Each stripe is one matrix product of the stripe's rows
of `G` with the rows from the stripe down, and one more pass over the genotypes:
from RAM without `-m`, from the file with `-m`. The stripes are sized to fit in
`-m`, so the total work is still the dense path's single `O(N^2 M)` Gram
product, split into pieces. With missing genotypes, the sites a pair has in
common, `n_ij = n - m_i - m_j + mm_ij`, need `mm_ij`, the sites both samples
miss. That count comes from the lists of missing calls, not from a second
`N x N` product. At biobank call rates nearly every site has a missing call,
but only a few samples miss each one. A site missing in a large share of the
stripe falls back to a product of missing indicators.

**Memory.** `-m` bounds the stripes, the genotype block and the working
buffers. Peak RSS was at most 20 MB above `-m` in every run below. Without `-m` the
genotypes are held in RAM, as in the dense path, and the stripes take up to
2 GB. A small `-m` only costs extra passes over the genotypes.

**Checks.** On every test set, the `.kin0` of a cutoff of -0.5 (all pairs)
matches `.kinship` character for character. That held in-core, with `-m` (1 to
9 stripes), with and without missing genotypes (including sites missing 30% of
calls), and for BED and PGEN input: 32 million pairs at N = 8,000, plus the
smaller sets in `tests/test_evaladmix_pairs.py`. The dense output is
byte-identical to the previous version.

Measured with 2 threads on x86-64. The data are simulated 3-way admixture with
5,000 unlinked sites, 1.4% missing calls (one sample in five at three times the
rate) and planted relatives, and `-k 2`. Time and peak RSS are for the
`--evaladmix` run alone:

| samples | dense `--evaladmix` | `--evaladmix-kin 0.0442` |
|---|---|---|
| 8,000 | 12.4 s, 1.45 GB, plus two 0.6 GB files | 6.4 s at `-m 1` (1.01 GB, 2 stripes); 11.7 s at `-m 0.3` (0.32 GB, 9 stripes) |
| 20,000 | needs 6 GB of RAM for its two N x N matrices, and writes two 3.6 GB files | 39.5 s at `-m 1` (1.02 GB, 5 stripes); 32.4 s at `-m 2` (2.00 GB, 3 stripes) |

At N = 20,000 every planted pair is found: 40/40 duplicates (mean 0.4996),
2,400/2,400 parent–offspring and full-sib pairs (0.2499), and 100/100 half sibs
(0.1244). One of the 2 x 10^8 unrelated pairs reaches 0.0442 (0.0446), which is
expected with 5,000 sites. The `.unrelated` set drops 1,041 of the 2,182 samples
with relatives. On the planted families that is the minimum possible: one per
duplicate, and the children of each family.

**Extrapolation (not measured).** The cost is `N^2 M` flops: 9 x 10^15 for
N = 245,000 and M = 150,000 LD-pruned sites. The stripes above ran at about
28 GFLOP/s per core, counting reading, imputation and output. At that rate, 64
cores would take about 1.4 hours, and `-m 64` makes 5 to 10 passes over the
9 GB `.bed`. Linking MKL or OpenBLAS
(`Makefile`) usually speeds up the products further.

## IBD sharing probabilities: `--evaladmix-ibd`

Kinship cannot tell parent–offspring from full sibs: both have `phi = 1/4`.
What differs is how the sharing is distributed. A parent and child share exactly
one allele IBD at every site (`k0, k1, k2 = 0, 1, 0`). Full sibs share none at a
quarter of the sites, one at half and both at a quarter (`1/4, 1/2, 1/4`).
`--evaladmix-ibd` estimates all three:

```bash
PCAone -b plink -k <K-1> --maf 0.05 -o pcs
PCAone -b plink -P pcs --evaladmix --evaladmix-ibd --maf 0.05 -o out
PCAone -b cohort -P pcs --evaladmix --evaladmix-ibd --evaladmix-kin 0.0442 -m 64 -o rel
```

The dense run also writes `out.k2` and `out.k0`, N x N with the same header as
`out.kinship`. `k1 = 4 phi - 2 k2` follows from `out.kinship` and `out.k2`.
As in `.kinship`, the diagonal is not an estimate (0 in `.k2`, 1 in `.k0`).
With `--evaladmix-kin`, `rel.kin0` gains three columns,
`... NSNP KINSHIP K0 K1 K2`. The log also counts how many first-degree pairs have
`k0 < 0.125` (parent–offspring) and how many have more (full sibs).
`.corres`, `.kinship`, `KINSHIP` and `.unrelated` are byte-identical with and
without the option.

**What it estimates.** `k2` comes from the homozygotes. PC-Relate (Conomos et
al. 2016, *AJHG* 98:127) codes each genotype by its dominance deviation. With
`pi` the fitted allele frequency of the sample at the site:

```
h = pi        if x = 0
    0         if x = 1/2
    1 - pi    if x = 1
```

For an outbred sample, `h` has mean `pi(1-pi)` and zero covariance with the
genotype. Between two samples, its covariance is `k2` times its variance. A pair
that shares one allele IBD agrees in `h` no more than by chance. So the
dominance residual `h - pi(1-pi)` is to `k2` what the genotype residual is to
`2 phi`. PCAone therefore computes `k2` as the evalAdmix statistic of the
dominance residuals, with the same steps as for the genotypes:

- the projection onto `[PC_1..PC_k, 1]`;
- per-sample centring;
- `cov2cor`;
- subtracting `chat`;
- rescaling for missing calls.

`pi` is the fit of the imputed genotypes, `f + Q Q' g`, bounded to
`[0.01, 0.99]`, as PC-Relate bounds it. A missing call has dominance residual 0.
Then

```
k0 = 1 - 4 phi + k2,    k1 = 4 phi - 2 k2
```

(`4 phi = k1 + 2 k2` and `k0 + k1 + k2 = 1`), using the clipped kinship of
`.kinship`. The estimates are not forced onto the simplex. They scatter around
the truth and can fall a little below 0 or above 1.

This differs from PC-Relate's `k2` in two ways:

- **Normalization.** PC-Relate divides by the variance its model predicts and
  then subtracts `f_i f_j` for inbreeding. PCAone divides by the observed
  variance, as evalAdmix does, and centres each sample. So a duplicate pair has
  `k2 = 1`. PC-Relate's formula, recomputed here without its
  small-sample correction, gave 0.91–0.99 on the simulations below.
- **Projection.** The projection takes out any dominance deviation shared along
  the PCs. On the 1000 Genomes panel (6 populations, `-k 5`), skipping the
  projection gives every pair `k2 ≈ +0.004`. Some sites depart from
  Hardy–Weinberg across the whole sample, and every pair reads that as shared
  homozygosity. With the projection, the unrelated pairs centre on 0
  (sd 0.004).

**Accuracy.** On the relateAdmix pedigree from [Accuracy](#accuracy) (K = 2,
`-k 1`, 102,185 sites after `--maf 0.05`), the mean per class against
RelateAdmix's constrained ML. RelateAdmix uses the ADMIXTURE `Q` and `P` shipped
with the data.

| relationship | truth k0 / k1 / k2 | PCAone k0 | k1 | k2 | RelateAdmix k0 | k1 | k2 |
|---|---|---|---|---|---|---|---|
| duplicate | 0 / 0 / 1 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 |
| parent–offspring | 0 / 1 / 0 | −0.002 | 1.006 | −0.003 | 0.000 | 1.000 | 0.000 |
| full sib | 0.25 / 0.5 / 0.25 | 0.249 | 0.511 | 0.239 | 0.264 | 0.496 | 0.240 |
| half sib | 0.5 / 0.5 / 0 | 0.510 | 0.488 | 0.001 | 0.557 | 0.442 | 0.000 |
| first cousin | 0.75 / 0.25 / 0 | 0.752 | 0.248 | 0.000 | 0.790 | 0.210 | 0.000 |
| second cousin | 0.9375 / 0.0625 / 0 | 0.931 | 0.066 | 0.002 | 0.954 | 0.045 | 0.001 |
| unrelated (7,815 pairs) | 1 / 0 / 0 | 1.005 | −0.003 | −0.002 | 0.999 | 0.001 | 0.000 |

On the 60 related pairs, the RMSE of `(k0, k1, k2)` is `(0.009, 0.012, 0.006)`
for PCAone and `(0.032, 0.031, 0.005)` for RelateAdmix. RelateAdmix's ML
attenuates `k1` for the distant classes, as it attenuates kinship in
[the comparison above](#how-it-compares-with-pc-relate-and-relateadmix). On the
unrelated pairs it is again near-exact, because it is pinned at the boundary of
the simplex: an RMSE of `(0.002, 0.002, 0.0003)` against PCAone's
`(0.012, 0.012, 0.007)`.

The classes that matter separate cleanly:

- parent–offspring `k0 <= 0.009` and `k2 <= 0.004`;
- full sibs `k0 >= 0.240` and `k2 >= 0.234`.

**Simulations.** These used the relateAdmix allele frequencies, 40,000 unlinked
sites, `-k 1` and 320 samples:

- **Admixed families.** In families of admixed and unadmixed parents, the class
  means of `k0` and `k2` are within 0.008 of the truth. That holds for
  duplicates, parent–offspring, full and half sibs, and first cousins.
- **Inbred samples.** Unrelated inbred samples (F of 1/16, 1/8 and 1/4) have
  `k2 = 0.002`.
- **Missing genotypes.** With 10% of calls missing at random, or 0–40% per
  sample, every class mean moves by at most 0.004.

**Cost.** One more Gram product, of the dominance residuals, of the same size
as the kinship one. The dense path holds a second `N x N` matrix and writes two
more files. In `--evaladmix-kin`, each stripe holds a second matrix of pairs,
so the same `-m` makes more stripes. The table shows N = 3,000 samples and
20,000 sites, 1% missing, with 2 threads:

| | without | with `--evaladmix-ibd` |
|---|---|---|
| dense | 11.3 s, 0.80 GB | 15.7 s, 0.87 GB |
| `--evaladmix-kin 0.05 -m 0.2` | 11.4 s, 0.22 GB, 3 stripes | 17.8 s, 0.22 GB, 4 stripes |

**Checks.** `tests/test_evaladmix_ibd.py` covers:

- **numpy reference.** The dense `.k2` matches a numpy reference to the printed
  digit.
- **Pair columns.** The `K0 K1 K2` columns of `--evaladmix-kin -0.5` match the
  dense matrices in-core and in stripes.
- **Run modes.** `-m`, PGEN and in-core give the same output.
- **Planted relatives.** Duplicates, parent–offspring and full sibs are told
  apart.

On the data above, `.corres`, `.kinship`, `.kin0` and `.unrelated` without the
option are byte-identical to the previous version, in-core and with `-m`.

**Limitations.** `k2` is a smaller signal than kinship, carried by the
homozygotes only. In order of how much they matter:

- **Noise.** `k2` scatters about as much as `corres`, twice as much as the
  kinship: sd 0.007 for the unrelated pairs of the pedigree. `k0` and `k1` add
  it to `4 phi`, so there they scatter by about 0.01.
- **Children of parents of different ancestry.** These depart from
  Hardy–Weinberg at their own fitted `pi`. A single `pi` per sample misses that
  their two alleles come from different populations. In the extreme case,
  children of parents from two different unadmixed populations:
  - half sibs read `k2 ≈ +0.03`;
  - full sibs read 0.27–0.28;
  - unrelated children of such couples read up to `+0.004`.

  PC-Relate's formula gives the same half-sib bias.
- **Small samples.** As with kinship, the fit uses the pair's own genotypes.
  With 285 samples, full sibs of an unadmixed population read `k2 = 0.232`.
  With 825, they read 0.246. `k0` is less affected (0.247 and 0.250), because
  the kinship is attenuated too.
- **Many close relatives.** In a small sample dense with close relatives, the
  projection spreads their signal. Unrelated pairs between members of different
  duplicate pairs read `k2 ≈ −0.05` (kinship −0.02) on the pedigree above.
- **Hidden missing calls.** A missing call at a site whose frequency is exactly
  0.5 cannot be detected (see [Missing genotypes](#missing-genotypes)). It
  counts as a heterozygote.

## Implementation notes

- **Genotype scale.** PCAone codes genotypes as `{0, 0.5, 1}` (`BED2GENO`), i.e.
  `g/2`. Since `bhat` and `chat` are each converted to correlations, the constant
  factor cancels and no rescaling is needed.

- **Centring.** Both paths read *centred* genotypes with missing calls imputed
  to the site mean (centred value 0): with `center = false` the in-core readers
  would leave missing calls as −9. The out-of-core block reader estimates `F` on
  the fly, so `F` and `centered_geno_lookup` are allocated before the loop. The
  centring is undone to recover heterozygosity. `A` and `gbar` need no
  correction: per-site centring subtracts `c_s * 1` from column `s`, and every
  resulting term carries a factor `(I-P)1 = 0` because the projection contains
  the intercept.

- **Finding the missing calls.** An imputed call is an exact 0 in the centred
  genotypes, and an observed one is 0 only if it equals the site frequency, which
  needs `f` in {0, 0.5, 1}. So the pair counts use every other site, with no
  change to the readers.

- **`read_usv()`.** The PC scores are read with `Utils::read_usv()`, which this
  branch also fixes — see [read-usv-fix.md](../dev/read-usv-fix.md). `-k` below the
  reference's PC count is a direct regression test for that bug.

## Limitations

- PLINK bed / PLINK2 pgen input only, diploid (`--haploid` is refused).
- `--evaladmix-ibd` doubles the `N x N` memory of the dense output and adds a
  second matrix of pairs to each stripe of `--evaladmix-kin`. See
  [its limitations](#ibd-sharing-probabilities---evaladmix-ibd) for the
  estimator's.
- **The dense output needs `N x N` memory and two `N x N` text files**; use
  `--evaladmix-kin` beyond a few tens of thousands of samples.
- **Missing genotypes cost the dense output a second `N x N` matrix and double
  the pass.** The
  rescaling by sites in common corrects the covariance, not the variances: the
  imputed calls still leave a small structural residual in each sample's
  variance. That is second order at the missingness tested here (up to 40% per
  sample), but a pair never genotyped at the same site has no estimate and is
  written as `nan` (left out of `.kin0`, with a warning).
- `--evaladmix-kin` holds the pair counts in float, so it refuses more than
  2^24 (16.8 million) sites with missing calls. It can only report pairs at or
  above the cutoff it was run with.
- Without `.eigvecs2` (e.g. scores given with `--read-U` from another program)
  only the number of rows of the reference is checked.
- `Dhat` assumes Hardy–Weinberg within individuals, so inbred samples need care.
  A sample with no heterozygous call is warned about.
- `--maf` with out-of-core is rejected by PCAone itself, unrelated to this
  feature.
- Tested on one dataset: 126 samples, K=2, PLINK and PGEN input, in-core and
  out-of-core, complete and with five missingness patterns. Speed and memory
  are measured on simulated data up to N = 8,000, and for `--evaladmix-kin` up
  to N = 20,000. The All of Us cost above is an extrapolation.
