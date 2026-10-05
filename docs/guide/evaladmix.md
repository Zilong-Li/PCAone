# Evaluating the PCA fit (evalAdmix)

`--evaladmix` computes the N x N correlation of residuals after predicting
each genotype from the PCs of a previous run of the same samples, in the same
order. Near-zero entries mean the PCs capture the structure; for related
samples the entry estimates twice their kinship. PCAone writes `.corres` and
`.kinship` (`.corres` divided by 2). Use K-1 PCs for K ancestral groups and apply
the same `--maf` in both runs.

```shell
PCAone -b example/plink -k 3 --maf 0.05 -o pcs
PCAone -b example/plink -P pcs --evaladmix --maf 0.05 -o eval
```

For biobank-scale cohorts, `--evaladmix-kin <cutoff>` computes the same
statistic in stripes that fit in `-m` and writes only the pairs with kinship at
or above the cutoff (`.kin0`), plus a maximal unrelated set (`.unrelated`, for
`--keep`), instead of the two N x N files:

```shell
PCAone -b cohort -P pcs --evaladmix --evaladmix-kin 0.0442 -m 32 -o rel
```

Kinship alone cannot tell parent–offspring from full sibs (both 1/4).
`--evaladmix-ibd` adds the probabilities of sharing 0, 1 and 2 alleles IBD:
`.k0` and `.k2` beside `.kinship`, or `K0 K1 K2` columns in `.kin0` with
`--evaladmix-kin`. Parent–offspring have `k0` near 0, full sibs near 1/4.

```shell
PCAone -b example/plink -P pcs --evaladmix --evaladmix-ibd --maf 0.05 -o eval
```

See [the evalAdmix guide](../methods/evaladmix.md) for the method and benchmarks.
