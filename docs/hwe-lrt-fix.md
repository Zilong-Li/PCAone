# A bug in `--inbreed`: the HWE likelihood ratio test is not on a common scale

`PCAone --inbreed 1` implements the structured-population Hardy–Weinberg test of
Meisner & Albrechtsen (2019, *Mol Ecol Resour* 19:1144). On our test data it
returns a **negative likelihood ratio statistic for 68% of markers**, the most
extreme being −9898. That is impossible: the alternative model nests the null,
so the maximised log-likelihood can only increase and `2(logAlt − logNull) ≥ 0`.

The consequence is a badly over-conservative filter. Against PCAngsd on the same
data, the same markers and the same K, PCAone removed **23%** of common markers
where PCAngsd removed **0.12%**.

## The cause

In `calc_inbreed_site_lrt()` (`src/InbredSites.cpp`), the **alternative** model
floors each genotype probability at `1e-4` and then renormalises:

```cpp
p0 = fmax(1e-4, (1.0 - PI(i,j)) * (1.0 - PI(i,j)) + Fadj);
p1 = fmax(1e-4, 2.0 * PI(i,j) * (1.0 - PI(i,j)) - (2.0 * Fadj));
p2 = fmax(1e-4, PI(i,j) * PI(i,j) + Fadj);
pSum = 1.0 / (p0 + p1 + p2);
p0 *= pSum;  p1 *= pSum;  p2 *= pSum;
```

while the **null** model, a few lines below, is left raw:

```cpp
logNull += log((1.0 - PI(i,j)) * (1.0 - PI(i,j)));
```

Unclamped, the three probabilities sum to exactly 1:

```
[(1−π)² + π(1−π)F] + [2π(1−π)(1−F)] + [π² + π(1−π)F]
  = (1−π)² + 2π(1−π) + π² + 2π(1−π)F − 2π(1−π)F = 1
```

so `pSum` is 1 and the renormalisation does nothing. It only bites when the
floor fires — and with **individual** allele frequencies it fires constantly,
because at a marker differentiated between groups many individuals have π close
to 0 or 1, making `(1−π)²` or `π²` smaller than `1e-4`.

When it fires, `p0 + p1 + p2 > 1`, so `pSum < 1` and **every** alternative
probability is scaled down — including the one belonging to the genotype that
was actually observed, which is usually the common homozygote with probability
near 1. The alternative is then penalised relative to an unpenalised null, and
`logAlt < logNull`.

The two models are no longer nested, because one has been renormalised and the
other has not. The statistic stops being a likelihood ratio.

The same floor is applied in `inbreed_coef_site()`, so the estimate of `F` is
distorted too: PCAone's F spans **−1.000 to +0.435** with many values pinned at
the −1 boundary, against PCAngsd's **−0.044 to +0.165** on identical input.

## The fix in this branch

`src/InbredSites.cpp`, three changes:

1. **Drop the per-term floor and the renormalisation from the alternative.**
   The probabilities sum to one by construction; a shared `PROB_EPS = 1e-12`
   floor is applied only so `log()` stays finite.
2. **Apply the same floor to the null**, so both sides are treated alike.
3. **Guard the statistic with `fmax(0.0, ...)`.** With 1 and 2 in place this
   should never bind; if it ever does it signals a remaining problem rather than
   emitting an impossible value.

The same `PROB_EPS` replaces `1e-4` in `inbreed_coef_site()`.

## Effect

Chromosome 22 of a simulated cohort: 10,000 individuals, 20 populations in 4
continental groups, 6,455 markers at MAF ≥ 5%, K = 3.

| | median F | removed at p < 1e-6 | negative LRT |
|---|---|---|---|
| PCAone v0.7.2 | −0.0436 | 22.97% | **4,399 / 6,455** |
| PCAone + this branch | −0.0436 | 8.60% | **0** |
| PCAngsd v1.36.4 | +0.0074 | 0.12% | 1,500 (small, near F = 0) |
| independent reimplementation | +0.0074 | 0.05% | 0 |

Ground truth for this dataset is **0.3%**: testing HWE within each of the 20
populations separately and summing the per-population chi-squares (20 df) — same
total sample size as the pooled test, so comparable power, but no Wahlund
effect. Only 0.3% of these markers genuinely deviate.

## A second, separate difference: π is estimated once, not iterated

The fix above removes the impossible statistics and most of the excess, but 8.6%
against 0.05% remains. That residual is **not** in `InbredSites.cpp`.

With called genotypes the posterior probability of heterozygosity is simply the
indicator — you either observed a heterozygote or you did not — so the EM update
reduces to

```
F = 1 − obsHet / expHet,     expHet = Σᵢ 2πᵢ(1−πᵢ)
```

`obsHet` is an observed count and is identical in both programs. Therefore **any
difference in F is a difference in π, and nothing else.**

PCAone forms π once, in `FileUSV::read_all()`:

```cpp
G = U * S.asDiagonal() * V.transpose();
G(j,i) = (G(j,i) + 2.0 * F(i)) * 0.5;
G(j,i) = fmin(fmax(G(j,i), 1e-4), 1.0 - 1e-4);
```

a single rank-K reconstruction of a matrix in which missing genotypes were
filled with the site mean. PCAngsd instead iterates: it estimates π, refills the
missing entries from the current π, refits, and repeats to convergence (its log
reports *"converged in 2 iterations using 3 eigenvectors"*). Mean-imputed entries
pull the fitted surface toward the population average; iterating replaces each
missing genotype with that individual's own expected frequency.

The `[1e-4, 1−1e-4]` clip compounds it: a π pinned at `1e-4` contributes
≈ 0.0002 to `expHet`, and because `F = 1 − obsHet/expHet`, a deflated
denominator inflates F.

Fixing this means iterating π in the SVD step, which is outside the scope of
this branch.

## Reproducing

```bash
# common markers only -- PCAngsd defaults to --maf 0.05 for good reason:
# below ~5% MAF the rank-K approximation of pi is poor and this test breaks down
# (median F was +0.83 with 96.6% spurious significance at MAF < 1%)
plink2 --bfile data --maf 0.05 --make-bed --out common

PCAone -b common -k 3 --svd 0 --printv -o svd
PCAone -b common --inbreed 1 -P svd -k 3 --svd 0 -o hwe
awk 'NR>1 && $3 < 0' hwe.hwe | wc -l        # negative LRTs: many before, 0 after
```

## Notes for anyone else running `--inbreed`

Four invocation requirements that are currently undocumented, three of which
fail by crashing rather than by erroring:

- `--inbreed` requires `--USV`; without it, a clean exit with a message.
- `-k` must match the SVD's k, or **segfault** (after only a warning).
- `-m` must not be used when the data fits in one block, or **abort**.
- The default window-based SVD needs more markers than its window; on small
  inputs, **floating point exception**. Use `--svd 0`.
