# Biobank analysis

These guides are for cohorts of tens of thousands to millions of samples, such
as UK Biobank or All of Us, where the genotypes may not fit in memory and an
`N x N` matrix certainly does not. PCAone reads the data block by block with
`-m/--memory` and keeps the memory bounded by it.

A typical workflow, from QC'd genotypes to PCs of the whole cohort:

```shell
# 1. PCs of all samples, out-of-core
PCAone -b cohort -k 20 -m 16 -o pcs
# 2. prune on the ancestry-adjusted LD
PCAone -b cohort -P pcs -k 10 --ld-r2 0.1 -m 16 -o pruned
plink2 --bfile cohort --extract pruned.ld.prune.in --make-bed --out cohort.pruned
# 3. relatives and a maximal unrelated set
PCAone -b cohort.pruned -P pcs -k 10 --evaladmix --evaladmix-kin 0.0442 -m 64 -o rel
# 4. PCs of the unrelated samples, then project the relatives onto them
plink2 --bfile cohort.pruned --keep rel.unrelated --make-bed --out unrel
plink2 --bfile cohort.pruned --remove rel.unrelated --make-bed --out related
PCAone -b unrel -k 20 -m 16 -o ref
PCAone -b related -P ref --project 2 -o proj
```

The numbers of PCs above are placeholders: choose them from the data, e.g. from
the [SNP loadings](../guide/plotting.md#snp-loadings) and the scatter check of
[`--evaladmix-kin`](relatedness.md#checks-in-the-log).

| Guide | What it does |
|---|---|
| [Out-of-core PCA](pca.md) | `-m` with BED, PGEN and BGEN; disk, threads and prefetching |
| [Projection](projection.md) | PCs of a subset, such as the unrelated samples, applied to everyone |
| [LD pruning and clumping](ld.md) | ancestry-adjusted LD out-of-core |
| [Relatedness](relatedness.md) | `--evaladmix-kin`: related pairs and an unrelated set within `-m` |
| [Plotting large cohorts](plotting.md) | density plots of the PCs of hundreds of thousands of samples |

The [selection scan](../small-n/selection.md) and the per-site and per-sample
[inbreeding](../small-n/hwe.md) analyses also run with `-m`. Each guide ends
with **Method details**.
