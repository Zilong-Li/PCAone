# `--evaladmix`: correlation of residuals from the top PCs

`--evaladmix` computes the evalAdmix statistic (Garcia-Erill & Albrechtsen 2020,
*Mol Ecol Resour* 20:936) using the analytic projection estimator of van Waaij
et al. (2023, *Genetics* 225:iyad157), with **PCA rather than an admixture model**
as the front-end.

```bash
PCAone -b plink -k <K-1> --evaladmix --maf 0.05 -o out
```

writes

- `out.corres` — the N x N correlation of residuals, and
- `out.kinship` — the same divided by 2, which is the kinship scale,

both with sample IDs from the `.fam`. `--evaladmix-k` selects how many of the
computed PCs enter the projection, so one run can produce many PCs for other
purposes while the statistic uses `K-1`.

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

**Use `K-1` PCs.** `method="standard"` appends an intercept, so `K-1` PCs span
the same dimension as `Q` for a `K`-population admixture model. One PC too many
costs roughly a factor of seven in RMSE, because with `K` real populations the
later PCs fit relatedness rather than ancestry.

## Why no residual matrix is needed

Assuming no missing genotypes, `Rtilde' Rtilde` can be written without ever
forming `R`:

```
Rtilde' Rtilde = (I-P) [ G'G - M gbar gbar' ] (I-P)
```

so the statistic depends on the data only through three quantities:

- `G'G`, the N x N genotype Gram matrix,
- `gbar`, the per-sample mean genotype,
- `Dhat`, the per-sample mean heterozygosity.

All three are accumulations over sites, so **one streaming pass** suffices.
Memory is `O(N^2)`, independent of the number of sites, which is what makes the
approach usable out-of-core. Cost is one `O(MN^2)` accumulation — the same order
as the covariance matrix a PCA already forms — plus an `O(N^3)` tail. On 126
samples x 102k sites the statistic takes **0.12 s** on top of a 0.54 s PCA.

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

Compared against the other ways of reaching the same statistic, on the same
data (RMSE over the 60 related pairs):

| route | RMSE | needs ADMIXTURE | genotypes in RAM |
|---|---|---|---|
| ADMIXTURE + evalAdmix EM | 0.00292 | yes | yes |
| ADMIXTURE + projection | 0.00429 | yes | yes |
| PCA (K−1) + projection, R | 0.00443 | no | yes |
| **`PCAone --evaladmix`** | **0.00251** | **no** | **no** |

Accuracy is essentially identical across all four — correlations with the truth
differ in the fourth decimal. Most of the RMSE spread is the `[-1,1]` clip: apply
it uniformly and the three projection routes land on 0.00251–0.00255 against the
EM route's 0.00292. The real differences are that the EM route refits `F` with
each individual excluded (so its work grows with sample size), the two PCA routes
need no admixture run at all, and only this one avoids holding the genotype
matrix in RAM.

All four overshoot equally at second cousins (~1.13x), which is a property of the
statistic near the noise floor rather than of any front-end.

**Validation.** Against `evalPCA()` in
[popgenDK/evalPopStructure](https://github.com/popgenDK/evalPopStructure):
r = 0.999955 over all 7875 pairs, identical RMSE on related non-duplicate pairs.
The only entries differing by more than 1.6e-04 are the ten duplicate pairs,
where this implementation applies the `[-1,1]` clip and the R reference does not.
In-core and out-of-core runs, and `--evaladmix-k 1` against a `-k 1` run, agree
exactly.

## Implementation notes

- **Genotype scale.** PCAone codes genotypes as `{0, 0.5, 1}` (`BED2GENO`), i.e.
  `g/2`. Since `bhat` and `chat` are each converted to correlations, the constant
  factor cancels and no rescaling is needed.

- **Centring.** The out-of-core block reader always returns *centred* genotypes,
  ignoring `params.center`, and estimates `F` on the fly — so `F` and
  `centered_geno_lookup` are allocated before the loop and the centring is undone
  to recover heterozygosity. `A` and `gbar` need no correction: per-site centring
  subtracts `c_s * 1` from column `s`, and every resulting term carries a factor
  `(I-P)1 = 0` because the projection contains the intercept. The in-core path
  runs with `center = false` and is already raw.

- **`read_usv()`.** The PC scores are read with `Utils::read_usv()`, which this
  branch also fixes — see [read-usv-fix.md](read-usv-fix.md). `--evaladmix-k` is
  a direct regression test for that bug.

## Limitations

- PLINK bed / PLINK2 pgen input only.
- **Complete data is assumed.** With missingness, evalAdmix skips sites absent in
  *either* member of a pair, so each pair sees a different site set, `G'G` is no
  longer a sufficient summary, and pairwise site counts would have to be
  accumulated alongside it.
- `Dhat` assumes Hardy–Weinberg within individuals, so inbred samples need care.
- `--maf` with out-of-core is rejected by PCAone itself, unrelated to this
  feature.
- Tested on one dataset: 126 samples, complete genotypes, K=2, PLINK input,
  in-core and out-of-core. Not tested on pgen, or at a scale where the `O(N^2)`
  memory matters.
