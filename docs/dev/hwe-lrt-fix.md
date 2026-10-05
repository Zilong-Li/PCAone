# Two bugs in `--inbreed`: the HWE test removed 23% of markers instead of 0.3%

`PCAone --inbreed 1` implements the structured-population Hardy–Weinberg test of
Meisner & Albrechtsen (2019, *Mol Ecol Resour* 19:1144). On our test data v0.7.2
removed **22.97%** of common markers where PCAngsd removed **0.12%** and the
truth is **0.30%**, and it returned a **negative likelihood ratio statistic for
68% of markers** (most extreme −9898).

There are two independent defects, one in each of two files. With both fixed
PCAone agrees with PCAngsd.

| | median F | removed at p<1e-6 | negative LRT | r(F, PCAngsd) |
|---|---|---|---|---|
| v0.7.2 as shipped | −0.0436 | **22.97%** | 4,399 | 0.168 |
| + fix 1 (`InbredSites.cpp`) | −0.0436 | 8.60% | 0 | 0.168 |
| + fix 2 (`FileUSV.cpp`) | **+0.0103** | **0.05%** | **0** | **0.751** |
| PCAngsd v1.36.4 | +0.0074 | 0.12% | — | — |
| stratified truth | | 0.30% | | |

Test data: chromosome 22 of a simulated cohort — 10,000 individuals, 20
populations in 4 continental groups, 6,455 markers at MAF ≥ 5%, K = 3. Ground
truth is a stratified test: HWE within each of the 20 populations, chi-squares
summed over 20 df. Same total sample size as the pooled test, so comparable
power, but no Wahlund effect. Only 0.3% of these markers genuinely deviate.

Both fixes work with the plain default SVD. Neither `--emu` nor `--pcangsd` is
needed.

---

# The model

Each individual *i* has its own allele frequency πᵢₛ at site *s*, obtained from a
rank-K decomposition of the genotypes. Under Hardy–Weinberg with a per-site
inbreeding coefficient Fₛ the genotype probabilities are

```
P(g=0) = (1−π)²   + π(1−π)F
P(g=1) = 2π(1−π)  − 2π(1−π)F  =  2π(1−π)(1−F)
P(g=2) = π²       + π(1−π)F
```

Two properties matter for what follows.

**They sum to one for every F.** The `+π(1−π)F` added to each homozygote is
exactly the `−2π(1−π)F` taken from the heterozygote:

```
(1−π)² + 2π(1−π) + π²  +  π(1−π)F + π(1−π)F − 2π(1−π)F  =  1 + 0  =  1
```

**F = 0 is nested inside F ≠ 0.** The null is the same family at F = 0, so the
likelihood ratio

```
T = 2( log L(F̂) − log L(0) )
```

is **non-negative by construction** — the maximised likelihood cannot be worse
than the likelihood at one particular point. `T < 0` is not a numerical
accident; it means the two likelihoods were not computed from the same family.

---

# Fix 1 — `src/InbredSites.cpp`: the null and the alternative were on different scales

## What was wrong, mathematically

In `calc_inbreed_site_lrt()` the alternative model floored each probability at
`1e-4` and then renormalised, while the null was left raw. Write the floored
values as p̃₀, p̃₁, p̃₂ and their sum as C. The code then evaluates

```
alternative:   log( p̃_g / C )        with C = p̃₀ + p̃₁ + p̃₂
null:          log( q_g )            with q₀+q₁+q₂ = 1, untouched
```

When no floor fires, p̃ = p and C = 1, so the renormalisation does nothing. When
a floor does fire, C > 1, and **every** alternative probability is divided by
C > 1 — including the one for the genotype actually observed. The statistic
becomes

```
T = 2 Σᵢ [ log p̃_gᵢ − log C − log q_gᵢ ]
  = 2 Σᵢ [ log p̃_gᵢ − log q_gᵢ ]  −  2 Σᵢ log C
```

with a penalty `−2 N log C` that has nothing to do with F. Once N is 10,000 that
term dominates, and T goes negative.

The floor fires constantly here because π is an **individual** allele frequency:
at a marker differentiated between groups many individuals have π near 0 or 1,
so `(1−π)²` or `π²` drops below `1e-4`. Hence 68% of markers.

## What changed in the code

The per-term floor and the renormalisation are removed from the alternative,
since the probabilities already sum to one. A single shared `PROB_EPS = 1e-12`
is applied to **both** models, purely so `log()` stays finite, and the statistic
is guarded at zero.

```diff
+// Shared floor for genotype probabilities, applied identically to the null
+// and the alternative model, solely to keep log() finite.
+static constexpr double PROB_EPS = 1e-12;
```

In `calc_inbreed_site_lrt()`:

```diff
-      p0 = fmax(1e-4, (1.0 - PI(i, j)) * (1.0 - PI(i, j)) + Fadj);
-      // p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) * (1.0 - F(jj)));
-      p1 = fmax(1e-4, 2.0 * PI(i, j) * (1.0 - PI(i, j)) - (2.0 * Fadj));
-      p2 = fmax(1e-4, PI(i, j) * PI(i, j) + Fadj);
-      // normalize
-      pSum = 1.0 / (p0 + p1 + p2);
-      p0 *= pSum;
-      p1 *= pSum;
-      p2 *= pSum;
+      p0 = fmax(PROB_EPS, (1.0 - PI(i, j)) * (1.0 - PI(i, j)) + Fadj);
+      p1 = fmax(PROB_EPS, 2.0 * PI(i, j) * (1.0 - PI(i, j)) * (1.0 - F(jj)));
+      p2 = fmax(PROB_EPS, PI(i, j) * PI(i, j) + Fadj);
```

The null gets the same floor, so both sides are treated alike:

```diff
-          logNull += log((1.0 - PI(i, j)) * (1.0 - PI(i, j)));
+          logNull += log(fmax(PROB_EPS, (1.0 - PI(i, j)) * (1.0 - PI(i, j))));
```

and the statistic cannot come back impossible:

```diff
-    T(jj) = 2.0 * (logAlt - logNull);
+    T(jj) = fmax(0.0, 2.0 * (logAlt - logNull));
```

The same `1e-4 → PROB_EPS` change is made in `inbreed_coef_site()`, so the floor
does not distort F either.

**Effect:** negative statistics 4,399 → 0, markers removed 22.97% → 8.60%. F is
unchanged, because with called genotypes the EM update reduces to
`F = 1 − obsHet/expHet` and the floor almost never touches an observed genotype.

### Negative statistics are reduced, not eliminated

On the dataset above none survive. That is not general. Re-measured on
`example/plink.chr1` (400 samples, 50,736 markers, K = 3), **5.45% of markers
still come out negative**, worst −19.52.

The reason is that `T ≥ 0` needs the alternative to be *maximised* over a family
containing the null. `log L(F)` is concave in F — every genotype probability is
linear in F — so the constrained maximiser always beats `F = 0`. But the F used
here is not that maximiser. With called genotypes the posteriors in
`inbreed_coef_site()` are deterministic, so the EM converges in one step to
`F = 1 − obsHet/expHet`, a moment estimator that matches the heterozygote cell
only. It can fit worse than `F = 0`.

Substituting the true constrained MLE removes all of them (minimum +2.5e-10) and
barely moves anything else — `cor(F, F_MLE) = 0.989`, and of 50,736 markers the
calls at p < 1e-6 go 280 → 278 — while costing several times the runtime. So the
moment estimator is kept, `fmax(0, ·)` stays, and the count is now reported:

```
2766 of 50736 sites have a negative likelihood ratio (most negative -19.5157)
and are reported as 0. ...
```

Such markers are never significant, so no call depends on this; it was worth
surfacing rather than silently zeroing.

### Relation to PCAngsd

PCAngsd's `loglike()` in `inbreed_cy.pyx` has the *same* asymmetry this fix
removes — the alternative is floored at 1e-4 and renormalised, the null is used
raw — bounds F only to [−1, 1], and applies no guard to the statistic. So this
is a fix relative to PCAngsd, not a change that reproduces it.

---

# Fix 2 — `src/FileUSV.cpp`: π was reconstructed at half the correct scale

This is the one that made the answer wrong rather than merely ill-defined.

## What was wrong, mathematically

`read_all()` built the individual allele frequencies as

```cpp
G = U * S.asDiagonal() * V.transpose();
G(j,i) = (G(j,i) + 2.0 * F(i)) * 0.5;
```

i.e. `π = (U·S·Vᵀ + 2f) / 2`. That is only correct if `U·S·Vᵀ` is a rank-K
approximation of the **centred genotype** matrix `X − 2f`. It is not, for two
separate reasons.

**(a) The decomposition is of the standardised matrix.** `Data::standardize_E()`
does

```cpp
const double sd = sqrt(f * (1.0 - f));
G.col(i) *= sqrt((double)params.ploidy) / sd;
```

so the matrix fed to the SVD is

```
Z = (X − 2f) · √2 / sd ,      sd = √(f(1−f))
```

A rank-K approximation of Z is therefore an approximation of the *standardised*
deviations. Recovering the genotype scale needs the inverse transform,
`× sd/√2`, before `2f` is added.

**(b) The stored singular values are half their true value.** On the same
matrix:

| | PC1 | PC2 | PC3 |
|---|---|---|---|
| PCAone `.sigvals` | 2390.75 | 1841.50 | 979.894 |
| true singular values of Z | 4781.51 | 3683.00 | 1959.79 |
| ratio | 2.0000 | 2.0000 | 2.0000 |

This is internally consistent — `.eigvals` (885.469) is exactly
`s_stored²/M` with M = 6455 — so the halving happens upstream of both outputs
and is not a write bug in one file. **Anything else that reconstructs from the
USV files is therefore affected too, not just `--inbreed`.**

The two compose. The reconstruction is short by `sd/√2` from (a) and by a
further factor of 2 from (b), so the correction is

```
π = ( U·S·Vᵀ · sd·√2  +  2f ) / 2
```

## Why this produced exactly the symptoms seen

Only half the structure was being applied, so π sat too close to the population
mean, and expected heterozygosity `Σᵢ 2πᵢ(1−πᵢ)` was too high. Since

```
F = 1 − obsHet / expHet
```

and `obsHet` is an observed count — identical in every implementation — an
inflated denominator drives F negative. Backing `expHet` out of each program's
F:

| | expHet | implied F |
|---|---|---|
| no structure at all, 2f(1−f)N | 3193.6 | +0.121 |
| PCAone, fix 1 only | 3031.8 | +0.074 |
| PCAngsd | 2826.6 | +0.007 |

PCAone sat nearest the *no structure* value precisely because it was applying
half of it.

It also explains two things that were otherwise puzzling:

- **Rare variants were catastrophic.** The missing factor is `sd = √(f(1−f))`,
  so the error grows as f → 0. At MAF < 1% median F was **+0.83** with 96.6% of
  markers spuriously significant. PCAngsd defaults to `--maf 0.05`; PCAone
  defaults to `--maf 0`, so it will run where the method cannot work.
- **More components made it worse** — each one adds another wrongly-scaled
  deviation.

## What changed in the code

In `read_all()`:

```diff
-        G(j, i) = (G(j, i) + 2.0 * F(i)) * 0.5;
+        {
+          const double sd = std::sqrt(F(i) * (1.0 - F(i)));
+          const double g = (sd > 1e-9) ? G(j, i) * sd * std::sqrt(2.0) : 0.0;
+          G(j, i) = (g + 2.0 * F(i)) * 0.5;
+        }
         G(j, i) = fmin(fmax(G(j, i), 1e-4), 1.0 - 1e-4);
```

and identically in the block-wise `read_block_initial()` used out of core:

```diff
-        G(j, i) = (G(j, i) + 2.0 * F(snp_idx)) * 0.5;
+        const double sd = std::sqrt(F(snp_idx) * (1.0 - F(snp_idx)));
+        const double g = (sd > 1e-9) ? G(j, i) * sd * std::sqrt(2.0) : 0.0;
+        G(j, i) = (g + 2.0 * F(snp_idx)) * 0.5;
```

`sd > 1e-9` guards a monomorphic site, where the deviation is zero anyway.

**Effect:** markers removed 8.60% → **0.05%**, median F −0.0436 → **+0.0103**,
correlation of F with PCAngsd 0.168 → **0.751**.

---

# What was ruled out

Recorded so nobody repeats the dead ends. None of these was the cause:

| candidate | evidence |
|---|---|
| number of components | F flat across K = 5, 20, 40 (+0.167, +0.177, +0.178) |
| missingness | median F 0.177 with 3% missing, 0.182 with none |
| SVD solver | `--svd 0` and the default window method give bit-identical output |
| iterating π | 0.03% removed either way |
| `read_usv()` transposition | not on this path — `FileUSV` uses `read_eigvecs()`, which reads correctly |

---

# Reproducing

```bash
plink2 --bfile data --maf 0.05 --make-bed --out common   # see the MAF note above

PCAone -b common -k 3 --svd 0 --printv -o svd
PCAone -b common --inbreed 1 -P svd -k 3 --svd 0 -o hwe

awk 'NR>1 && $3 < 0'    hwe.hwe | wc -l   # negative LRTs: 4399 before, 0 after
awk 'NR>1 && $2 < 1e-6' hwe.hwe | wc -l   # markers removed: 22.97% before, 0.05% after
```

Reference: `pcangsd -p common -e 3 --inbreed-sites -o ref`.

---

# Undocumented invocation requirements

Four, three of which fail by crashing rather than erroring:

- `--inbreed` requires `--USV`; without it, a clean exit with a message.
- `-k` must match the SVD's k, or **segfault** (after only a warning).
- `-m` must not be used when the data fits in one block, or **abort**.
- The default window-based SVD needs more markers than its window; on small
  inputs, **floating point exception**. Use `--svd 0`.
