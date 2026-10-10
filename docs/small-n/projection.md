# Projection

Projection places new samples on the PCs of a reference, without recomputing
them: ancient or low-coverage samples onto a panel of modern ones, or a study
cohort onto 1000 Genomes.

First run PCAone on the reference samples. It writes the loadings and the
`.mbim` (variants and allele frequencies) by default:

```shell
PCAone -b example/ref -k 10 -o ref
```

Pass the reference prefix with `-P/--USV` to read its loadings, singular
values and allele frequencies. PCAone matches the overlapping markers against
`.mbim` and corrects flipped alleles. Both datasets must use the same genome
build and compatible alleles.

```shell
## --USV: prefix to .eigvecs, .sigvals, .loadings, .mbim
## --project: the projection method, see below
PCAone -b example/new \
       -P ref \
       --project 2 \
       -o new
```

The target's PC scores are written to `new.eigvecs`, in the order of its
`.fam`. Plot them over the reference's `ref.eigvecs` to see where the new
samples fall.

The projection uses every PC in the reference; add `-k` to use only the leading
ones (e.g. `-k 3`). This holds for every analysis that reads a reference with
`-P/--USV`: `--project`, `--selection`, `--inbreed`, `--evaladmix` and the
[LD analyses](ld.md), so one reference run with a generous `-k` serves them all.

## Choosing a method

| `--project` | input | missing calls |
|---|---|---|
| `1` | PLINK, PGEN | replaced by the reference allele frequency |
| `2` | PLINK, PGEN | left out: each sample is fitted on its own called sites |
| `3` | BEAGLE | genotype likelihoods, by EM |

Use `--project 2` for called genotypes with missing data, such as ancient DNA
or a different array. `--project 1` is the classic projection and suits a
target with few missing calls; it warns when there are some. A sample with
fewer called sites than PCs gets NA with `--project 2`, with a warning.

## Uncertainty of projected samples

A sample with few called sites has an uncertain position. `--project-bootstrap`
resamples the matched SNPs with replacement after the baseline `--project 2`
projection and projects each sample again, once per replicate:

```shell
PCAone -b example/new \
       -P ref \
       --project 2 \
       --project-bootstrap 100 \
       -o new
```

PCAone writes per-sample and per-PC summaries to `.proj.bootstrap.tsv`, and
the covariance of each pair of PCs to `.proj.bootstrap.cov.tsv`, which can be
used to draw uncertainty ellipses around the baseline estimates. Add
`--project-bootstrap-save` to write the coordinates of every replicate to
`.proj.bootstrap.eigvecs`.

| file | columns |
|---|---|
| `.proj.bootstrap.tsv` | `sample pc baseline mean bootstrap_se min max rmsd` |
| `.proj.bootstrap.cov.tsv` | `sample pc_x pc_y baseline_x baseline_y mean_x mean_y var_x var_y cov corr` |
| `.proj.bootstrap.eigvecs` | `replicate sample PC1 PC2 ...` |

`sample` is the 1-based position of the sample in the target, and `rmsd` is
the root mean square distance of the replicates from the baseline.

## Genotype likelihoods

For BEAGLE genotype likelihood input, use `--project 3` to perform
genotype-likelihood aware projection with an EM procedure:

```shell
PCAone -G example/target.beagle.gz \
       -P ref \
       --project 3 \
       -o target
```

The target is scaled as the reference PCA was (recorded in its `.sigvals`), so
`--scale` is not needed here. The reference can be a PCA of called genotypes or
a PCAngsd run.

**NB:** The two BEAGLE alleles must be the two alleles of the reference `.mbim`,
in either order. BEAGLE likelihoods count `allele2`, which PCAone matches to the
reference's counted allele `A1` (5th column of `.mbim`), and sites where they
are swapped are flipped automatically; sites with other alleles are dropped.
Use angsd `-doMajorMinor 3` with a `-sites` file taken from the reference
`.mbim` to fix the alleles.

## Method details

The reference PCA decomposed `X ≈ U S V'`, with one row of `X` per sample and
one column per site. Projection keeps `V` and `S` fixed and solves for the new
rows of `U`.

**Matching.** The target's sites are matched to the reference `.mbim` by
chromosome, position and alleles. Sites present in only one are dropped. At a
site where the target counts the other allele, the row of `V` changes sign.
The reference allele frequencies `f` from the `.mbim` centre the target, and
its genotypes are put on the scale the reference PCA used, which `.sigvals`
records (`scale`, `ploidy` and `gscale` on its header line). So each target
sample is projected on its own, independently of the other target samples.

**`--project 1`** multiplies by the loadings, `u = g V S^-1`, with the missing
calls of `g` at the reference frequency, i.e. 0 once centred. A missing call
pulls the sample towards the origin, which is why `--project 2` is preferred
when there are many.

**`--project 2`** solves the least squares problem `(V S) u ≈ g` for each
sample over its called sites only. When every site of the reference is matched
and called, the columns of `V` are orthonormal and it equals `--project 1`.

**`--project 3`** alternates two steps, starting from the least-squares
projection of the expected genotypes:

- E-step: the individual allele frequency of the sample at each site is
  `pi = f + (u S V')_j / a_j`, bounded to `[1e-4, 1 - 1e-4]`. With it as the
  prior, the genotype likelihoods give the posterior expected genotype.
  `a_j` is the per-site factor of the reference's scaling: `sqrt(ploidy)/sd`
  for a standardized reference, 2 for a PCAngsd reference (centred dosages),
  1 for an unscaled one. PCAngsd's own algorithm is the `a_j = 2` case.
- M-step: least squares of the expected genotypes on `V S`, as in
  `--project 2`.

It stops when the relative change of the scores is below `--tol-em`, or after
`--maxiter` iterations. Equivalent 0..1 and 0..2 references give identical
coordinates.

**Bootstrap.** Each replicate draws as many sites as were matched, with
replacement, and solves the least squares problem of `--project 2` on them,
weighted by how often each site was drawn. The design `V S` is fixed, so a
replicate has no arbitrary sign and is not flipped towards the baseline. The
draws follow `--seed` and are the same on every platform.
