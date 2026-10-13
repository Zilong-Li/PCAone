# Projection at biobank scale

In a biobank, projection puts every sample on PCs computed from a subset of
it. Two common designs:

- **PCs of the unrelated samples.** Close relatives pull PCs towards their
  families, so compute the PCs from a maximal unrelated set and project the
  relatives onto them, so that every sample has comparable coordinates.
- **PCs of a reference panel.** Project the cohort onto the PCs of a panel
  with known ancestry, such as 1000 Genomes or HGDP, to place it on the
  panel's axes and label its ancestry.

```shell
# the unrelated set, e.g. from --evaladmix-kin (see Relatedness)
plink2 --bfile cohort --keep rel.unrelated --make-bed --out unrel
plink2 --bfile cohort --remove rel.unrelated --make-bed --out related
# PCs of the unrelated samples, out-of-core; writes the loadings by default
PCAone -b unrel -k 20 -m 16 -o ref
# project the relatives, out-of-core
PCAone -b related -P ref --project 2 -m 16 -o proj
```

`proj.eigvecs` then holds the coordinates of the relatives on the axes of
`ref.eigvecs`. The projection methods, the matching of variants and the
bootstrap are described in [Projection](../small-n/projection.md); everything
there applies here.

## Memory: stream the target in blocks

In-core, projection reads the genotypes of the matched variants into memory:
`9 x N x M` bytes (a value and a missing flag per genotype), e.g. 45 GB for
50,000 samples and 100,000 variants. With `-m`, `--project 1` and
`--project 2` read the target from the PLINK or PGEN files a block of variants
at a time and keep only sums per sample, so a target of any size is projected
in one run, within `-m` GB. There is no need to split it into batches of
samples.

| | held per sample, `K` PCs | passes over the genotypes |
|---|---|---|
| `--project 1` | the `K` scores | one |
| `--project 2` | `(K+1)(K+2)/2` numbers: 66 for `K = 10`, 0.5 GB per 1M samples | one, or one per group of samples |
| `--project-bootstrap R` | about `R` times as many: 5.3 GB per 100,000 samples for `R = 100`, `K = 10` | one, or one per group of samples |

When the sums of all the samples and a block of variants do not fit in `-m`,
the samples are taken in groups, with one pass over the genotypes each. The log
gives the plan, e.g.

```text
out-of-core projection by -m/--memory: blocksize = 3693 , nblocks = 6 , samples per pass = 60000 , passes = 1
```

The scores, the bootstrap files and the order of their rows are those of the
in-core run, up to rounding (in the tests, every digit written agrees).
`--project 3` reads genotype likelihoods and runs in-core only. LD-pruned
variants, which the PCs need anyway, keep `M` small and each pass short.

## Shrinkage of projected samples

When the variants far outnumber the reference samples, projected samples are
pulled towards the origin compared with the reference samples on the same PCs
(Lee et al. 2010, *Ann Stat* 38:3605). The effect is strongest on the later
PCs, with small eigenvalues. A large reference, such as the unrelated part of
a biobank, makes it small for the leading PCs; with a small reference panel,
compare projected and reference samples of the same ancestry before reading
distances between them.

## Method details

See [Projection: method details](../small-n/projection.md#method-details) for
how each method solves for the scores. Three points matter at this scale:

- **Independence of samples.** The target is centred by the allele
  frequencies of the reference `.mbim` and scaled as the reference PCA was
  (recorded in `.sigvals`), never by statistics of the target itself. Each
  sample's scores are a function of its own genotypes and the reference only,
  which is why streaming the target with `-m`, in groups of samples or not,
  gives the same coordinates as one in-core run.
- **Matched variants.** Only the variants in both the target and the
  reference `.mbim` are used, matched by chromosome, position and alleles. With
  `--project 1` the scores are `g V S^-1` over the matched variants, which
  shrinks every sample when many reference variants are missing from the
  target. `--project 2` fits each sample on its own called variants, so it is
  the safer default when the two were genotyped or imputed differently.
- **Out-of-core least squares.** In-core, `--project 2` takes the pivoted QR
  of the rows of `V S` at the called sites of a sample. With `-m` it solves the
  same least squares from its normal equations, `A z = b`, with
  `A = sum v_j v_j'` and `b = sum g_j v_j` over the called sites `j` (`v_j` the
  row of `V` of site `j`), and `u = S^-1 z`. `b` is a product of each block
  with `V`. `A` is the sum over all the sites of a block less the sites the
  sample misses, or the sum over the called ones when it misses most of the
  block, so the work grows with the missing calls rather than with
  `N x M x K^2`. The bootstrap accumulates the same sums with each site
  weighted by the number of times a replicate drew it, all replicates in one
  pass, with the draws of the in-core run (`--seed`).
