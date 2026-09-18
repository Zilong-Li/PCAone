# Two bugs in `--inbreed`: the HWE test removes 23% of markers instead of 0.3%

`PCAone --inbreed 1` implements the structured-population Hardy–Weinberg test of
Meisner & Albrechtsen (2019, *Mol Ecol Resour* 19:1144). On our test data v0.7.2
removed **22.97%** of common markers where PCAngsd removed **0.12%** and the
truth is **0.30%**, and it returned a **negative likelihood ratio statistic for
68% of markers** (most extreme −9898), which is impossible when the alternative
nests the null.

There are two independent defects. With both fixed, PCAone agrees with PCAngsd.

| | median F | removed at p<1e-6 | negative LRT | r(F, PCAngsd) |
|---|---|---|---|---|
| v0.7.2 as shipped | −0.0436 | **22.97%** | 4,399 | 0.168 |
| + fix 1 (LRT) | −0.0436 | 8.60% | 0 | 0.168 |
| + fix 2 (π scale) | **+0.0103** | **0.05%** | **0** | **0.751** |
| PCAngsd v1.36.4 | +0.0074 | 0.12% | — | — |
| stratified truth | | 0.30% | | |

Test data: chromosome 22 of a simulated cohort, 10,000 individuals, 20
populations in 4 continental groups, 6,455 markers at MAF ≥ 5%, K = 3. Ground
truth is a stratified test — HWE within each of the 20 populations, chi-squares
summed over 20 df — which has the same total sample size as the pooled test but
no Wahlund effect. Only 0.3% of these markers genuinely deviate.

Both fixes work with the plain default SVD; neither `--emu` nor `--pcangsd` is
needed.

---

## Fix 1 — `calc_inbreed_site_lrt()`: the two models are not on a common scale

In `src/InbredSites.cpp` the **alternative** model floors each genotype
probability at `1e-4` and then renormalises:

```cpp
p0 = fmax(1e-4, (1.0 - PI(i,j)) * (1.0 - PI(i,j)) + Fadj);
p1 = fmax(1e-4, 2.0 * PI(i,j) * (1.0 - PI(i,j)) - (2.0 * Fadj));
p2 = fmax(1e-4, PI(i,j) * PI(i,j) + Fadj);
pSum = 1.0 / (p0 + p1 + p2);
p0 *= pSum;  p1 *= pSum;  p2 *= pSum;
```

while the **null**, a few lines below, is left raw:

```cpp
logNull += log((1.0 - PI(i,j)) * (1.0 - PI(i,j)));
```

Unclamped the three sum to exactly 1:

```
[(1−π)² + π(1−π)F] + [2π(1−π)(1−F)] + [π² + π(1−π)F] = 1
```

so `pSum` is 1 and the renormalisation is a no-op. It only bites when the floor
fires — and with **individual** allele frequencies it fires constantly, because
at a differentiated marker many individuals have π near 0 or 1, making `(1−π)²`
or `π²` smaller than `1e-4`.

When it fires, `p0 + p1 + p2 > 1`, so `pSum < 1` and every alternative
probability is scaled down, including the one for the genotype actually
observed — usually the common homozygote with probability near 1. The
alternative is penalised against an unpenalised null, `logAlt < logNull`, and
the statistic goes negative. The models are no longer nested once one has been
renormalised and the other has not.

**Fix:** drop the per-term floor and the renormalisation; apply a shared
`PROB_EPS = 1e-12` to both models purely to keep `log()` finite; guard the
statistic with `fmax(0.0, …)`. The same floor replaces `1e-4` in
`inbreed_coef_site()`.

---

## Fix 2 — `FileUSV::read_all()`: π is reconstructed at half the correct scale

This is the one that made the test wrong rather than merely ill-defined. Two
scale factors are missing, and they compound.

```cpp
G = U * S.asDiagonal() * V.transpose();
G(j, i) = (G(j, i) + 2.0 * F(i)) * 0.5;      // <- no rescaling
G(j, i) = fmin(fmax(G(j, i), 1e-4), 1.0 - 1e-4);
```

**(a) The reconstruction is on the standardised scale.**
`Data::standardize_E()` divides each column by `sd = sqrt(f(1-f))` and
multiplies by `sqrt(ploidy)`, so `U·S·Vᵀ` is not on the genotype scale and must
be returned to it before `2f` is added.

**(b) The stored singular values are half their true value.** On the same
matrix:

| | PC1 | PC2 | PC3 |
|---|---|---|---|
| PCAone `.sigvals` | 2390.75 | 1841.50 | 979.894 |
| true singular values | 4781.51 | 3683.00 | 1959.79 |
| ratio | 2.0000 | 2.0000 | 2.0000 |

This is internally consistent — `.eigvals` (885.469) equals `s_stored²/M` — so
the halving is upstream of both outputs, not a write bug in one file. **It may
therefore affect every other consumer of the USV files, not only `--inbreed`.**

**Fix:** multiply the reconstruction by `sd * sqrt(2)` — that is `sd/sqrt(2)`
for (a) and a further `×2` for (b) — in both `read_all()` and the block-wise
`read_block_initial()` used out of core.

### Why this produced exactly the symptoms observed

Only half the structure was being applied, so π sat too close to the population
mean and expected heterozygosity was too high. Backing the implied `expHet` out
of each program's F (`F = 1 − obsHet/expHet`, and `obsHet` is an observed count
identical in all of them):

| | expHet | implied F |
|---|---|---|
| no structure at all, 2f(1−f)N | 3193.6 | +0.121 |
| PCAone, fix 1 only | 3031.8 | +0.074 |
| PCAngsd | 2826.6 | +0.007 |

PCAone sat nearest the no-structure value because it was applying half the
structure. This also explains:

- **Rare variants were catastrophic.** The error scales as `1/sd`, so it is
  worst where `sd` is smallest. At MAF < 1% median F was **+0.83** with 96.6% of
  markers spuriously significant. PCAngsd defaults to `--maf 0.05`; PCAone
  defaults to `--maf 0`, so it will happily run where the method cannot work.
- **More components made it worse** — each added a further wrongly-scaled
  deviation.

---

## What was ruled out

Recorded so nobody repeats it: the residual was **not** caused by the number of
principal components (F is flat across K = 5, 20, 40), by missingness (0.177 vs
0.182 with and without), by the SVD solver (`--svd 0` and the default window
method give bit-identical output), by iterating the individual allele
frequencies (0.03% either way), or by the `read_usv()` row-major fix — which is
not on this path at all, since `FileUSV` uses `read_eigvecs()`.

---

## Reproducing

```bash
plink2 --bfile data --maf 0.05 --make-bed --out common   # see the MAF note above

PCAone -b common -k 3 --svd 0 --printv -o svd
PCAone -b common --inbreed 1 -P svd -k 3 --svd 0 -o hwe

awk 'NR>1 && $3 < 0' hwe.hwe | wc -l                     # negative LRTs
awk 'NR>1 && $2 < 1e-6' hwe.hwe | wc -l                  # markers removed
```

Compare against `pcangsd -p common -e 3 --inbreed-sites -o ref`.

---

## Undocumented invocation requirements

Four, three of which fail by crashing rather than erroring:

- `--inbreed` requires `--USV`; without it, a clean exit with a message.
- `-k` must match the SVD's k, or **segfault** (after only a warning).
- `-m` must not be used when the data fits in one block, or **abort**.
- The default window-based SVD needs more markers than its window; on small
  inputs, **floating point exception**. Use `--svd 0`.
