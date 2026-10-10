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
# project the relatives
PCAone -b related -P ref --project 2 -o proj
```

`proj.eigvecs` then holds the coordinates of the relatives on the axes of
`ref.eigvecs`. The projection methods, the matching of variants and the
bootstrap are described in [Projection](../small-n/projection.md); everything
there applies here.

## Memory: project in batches of samples

Projection reads the target in-core: the genotypes of the matched variants
take `8 x N x M` bytes, e.g. 40 GB for 50,000 samples and 100,000 variants.
`-m` is ignored. Each sample is projected on its own, from the reference's
allele frequencies, loadings and singular values, so a large target can be
split into batches of samples and projected batch by batch:

```shell
# split the target into batches of 20,000 samples
split -l 20000 -d related.fam batch_
for f in batch_*; do
  plink2 --bfile related --keep $f --make-bed --out related_$f
  PCAone -b related_$f -P ref --project 2 -o proj_$f
done
cat proj_batch_*.eigvecs > proj.eigvecs
```

Keep the batches in the order of `split`, so the rows of the concatenated
`.eigvecs` follow `related.fam`. LD-pruned variants, which the PCs need anyway,
also keep `M` small.

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
how each method solves for the scores. Two points matter at this scale:

- **Independence of samples.** The target is centred by the allele
  frequencies of the reference `.mbim` and scaled as the reference PCA was
  (recorded in `.sigvals`), never by statistics of the target itself. Each
  sample's scores are a function of its own genotypes and the reference only,
  which is why batches give the same coordinates as one run.
- **Matched variants.** Only the variants in both the target and the
  reference `.mbim` are used, matched by chromosome, position and alleles. With
  `--project 1` the scores are `g V S^-1` over the matched variants, which
  shrinks every sample when many reference variants are missing from the
  target. `--project 2` fits each sample on its own called variants, so it is
  the safer default when the two were genotyped or imputed differently.
