# Ancestry-adjusted LD without a residual matrix

Since v0.8.0 the LD analyses (`-R/--print-r2`, `--ld-r2`, `--clump`) read the
genotypes and remove the PCs of a previous run as they are read:

```bash
PCAone -b plink -k 3 -o pcs                                  # the PCs
PCAone -b plink -P pcs --ld-r2 0.8 --ld-bp 1000000 -o adj    # LD on the residuals
```

This replaces the old two-step route, where `-D/--ld` wrote a `.residuals`
matrix `G - U S V'` and `-B/--binary` read it back. This note explains why the
residuals `(I - Q Q') G` are the same thing, where the two differ, and how the
new path is computed.

## The residuals are a projection

`G` is the N x M matrix of centred genotypes, one column per site. Missing calls
are imputed to the site mean, i.e. 0. `U` (N x k) holds the top k PCs, with
orthonormal columns: `U'U = I`.

Regress one site `g` on the PCs. The coefficients are `b = (U'U)^-1 U' g = U' g`,
and the residual is

```
r = g - U b = (I - U U') g
```

Over all sites, `R = (I - U U') G`. This is what "ancestry adjusted" means: what
is left of each site once the PCs are regressed out.

## Why it equals `G - U S V'`

Write the full SVD of `G` as the top k components plus the rest:

```
G = U S V' + U2 S2 V2'        with   U' U2 = 0
```

Multiply on the left by `U'`:

```
U' G = (U'U) S V' + (U'U2) S2 V2' = S V'
```

so

```
G - U S V' = G - U (U' G) = (I - U U') G
```

Row by row, row `j` of `V` is `S^-1 U' g_j`. `V` never adds information:
it can always be recomputed from `U` and the genotypes. That is why the new path
needs only the `.eigvecs`, and neither `.sigvals`, `.loadings` nor `.mbim`.

## Where the two can differ

- **The SVD method.** The identity needs `S V' = U' G` exactly. It holds for an
  exact truncated SVD, and for any method whose final `V` is the projection of
  `G` onto `U`: IRAM (`V = G' U S^-1`) and a Halko-style last step (the SVD of
  `B = Q'G`). If an approximate method's `V` does not satisfy it, `G - U S V'`
  keeps part of the PCs, while `(I - U U') G` removes them exactly, whatever
  produced `U`. The old `--ld` recommended `--svd 1` for this reason.
- **Standardisation.** A default PCA decomposes `G D`, where `D` is a diagonal
  matrix of per-site scales `sqrt(ploidy) / sd`. Then `U S V'` approximates
  `G D`, not `G`, and `G - U S V'` mixes two scales. This is why `-D/--ld` ran
  its PCA unstandardised. The projection does not care:

  ```
  (I - U U') (G D) = ((I - U U') G) D
  ```

  The two differ only by a per-site scale, and a correlation ignores that. The
  PCs themselves do change with standardisation. Add `--scale 0` to the PCA to
  get the PCs `-D/--ld` used, and hence the same R2.
- **Centring.** Every column of `G` sums to 0, so the PCs are orthogonal to the
  intercept and the residuals stay centred. They are re-centred anyway, which
  only removes rounding.
- **Precision of `.eigvecs`.** The file keeps 6 significant digits, so `U'U = I`
  only to about 1e-6. PCAone uses `Q` from a QR decomposition of `U`. `Q` spans
  the same PCs and `Q Q'` is an exact projector; for an exactly orthonormal `U`,
  `Q Q' = U U'`.
- **`--emu`.** The old in-core path subtracted the PCs from the EM-imputed
  genotypes. The new path uses mean imputation, as every other LD path does,
  including the old out-of-core one.

## How it is computed

`LDColumns` (`src/LD.cpp`) holds the residuals of the sites in memory: the whole
matrix in-core, or with `-m` the two consecutive blocks the current windows
need. Each column is scaled to unit norm once, so the correlation of two sites is
the dot product of their columns. A site with no variance (monomorphic, or fully
explained by the PCs) becomes a zero column, and its R2 is reported as 0.

The correlations of a batch of windows come from one matrix product,
`C = X_rows' X_cols`, rather than one dot product per pair. That reuses every
column while it is in cache.

- **Pruning** takes the next windows whose lead site is still kept. It
  multiplies them against the kept sites of their span, then applies the greedy
  rule (remove the lower-MAF site of each pair above `--ld-r2`) in window order.
  This makes the same decisions as one window at a time. The batch grows while
  leads survive and shrinks while they are pruned away.
- **`-R`** multiplies consecutive windows over their span. It then formats and
  gzip-compresses the lines on all threads, one gzip member per thread. The
  `.ld.gz` is a series of members, which `zcat`, `gzip -d`, R and Python read as
  one file.
- **Clumping** reads its target sites and takes dot products.

In-core and `-m` go through the same code. A window must fit in two blocks of
`-m`; otherwise PCAone asks for a larger `-m` or a smaller `--ld-bp`.

## Checks

- **Old vs new, the same PCs.** Both paths used the IRAM `.eigvecs` of
  `example/plink.chr1` (400 x 50,736).
  - R2 is equal to the printed precision (max difference 1e-6).
  - Pruning at r2 0.2/0.8/0.99 and clumping are identical, in-core and `-m`.
  - With `--ld-stats 1` and with `--maf 0.1`, R2 is also equal to the printed
    precision.
- **Independent reference.** numpy decodes the `.bed`, centres it, and projects
  with `Q` from `numpy.linalg.qr`. It matches PCAone to 5e-7.
- **PCs.** `--scale 0` reproduces the PCs of the old `-D/--ld` exactly.
- **Speed.** On 4 threads and simulated data with N x M = 1e8 (N = 2,000 to
  20,000, 1 Mb windows, up to 50M pairs), the new path reads 16x less data than
  a `.residuals` file (2 bits against 4 bytes per genotype). The batched
  products make it faster still:
  - pruning is 1.4-3.4x faster in-core and 5-14x faster with `-m`;
  - `-R` is 2.5-10x faster in-core and up to 38x faster with `-m`;
  - the output is bit-identical before and after.
